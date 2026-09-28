#!/usr/bin/env python3
"""
Validate the MOE score's MOE>=4 prioritization threshold against ClinVar
pathogenic/likely-pathogenic (P/LP) variant status, as an independent
genetics-derived label -- addressing reviewer request R1.3 ("provide
stronger independent validation ... reporting precision, recall, odds
ratios and confidence intervals").

ClinVar P/LP status is independent of how the MOE score itself is built
(GO/KEGG/Reactome pathway enrichment vs. clinical variant curation), so it
is a genuine external check rather than testing the score against data it
was derived from.

For every MOE score cutoff (>=0 through >=5), a gene is called "predicted
positive" if MOE_score >= cutoff, and "true positive" if it carries >=1
ClinVar P/LP variant (or, with --strict, >=1 P/LP variant with >=1-star
review status). Reports precision, recall (sensitivity), specificity,
accuracy, odds ratio with 95% CI (Haldane-Anscombe corrected when any
confusion-matrix cell is zero) and a Fisher's exact test p-value per
cutoff, plus ROC and precision-recall curves treating MOE_score as a
6-level ordinal classifier score.

Input: the merged gene master table (ndd_master_table.tsv), which must
have MOE_score, n_pathogenic_likely_pathogenic and n_plp_ge1star columns.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def confusion_counts(pred_pos, true_pos):
    tp = int(np.sum(pred_pos & true_pos))
    fp = int(np.sum(pred_pos & ~true_pos))
    fn = int(np.sum(~pred_pos & true_pos))
    tn = int(np.sum(~pred_pos & ~true_pos))
    return tp, fp, fn, tn


def odds_ratio_ci(tp, fp, fn, tn):
    """Odds ratio with 95% CI on the log scale; Haldane-Anscombe +0.5
    correction applied uniformly when any cell is zero, to avoid an
    undefined or infinite OR."""
    a, b, c, d = tp, fp, fn, tn
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_est = (a * d) / (b * c)
    log_or = np.log(or_est)
    se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    ci_low = np.exp(log_or - 1.96 * se_log_or)
    ci_high = np.exp(log_or + 1.96 * se_log_or)
    return or_est, ci_low, ci_high


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv")
    ap.add_argument("--strict", action="store_true",
                     help="Use n_plp_ge1star (>=1-star reviewed P/LP variants) as the "
                          "ClinVar positive label instead of n_pathogenic_likely_pathogenic")
    ap.add_argument("--class-filter", default=None,
                     help="Restrict to one Class value (curated, candidate, no_evidence). "
                          "Default: all genes.")
    ap.add_argument("--output-dir", default=".", help="Directory for output table and figure")
    args = ap.parse_args()

    if not os.path.exists(args.master_table):
        print(f"Error: File not found: {args.master_table}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.master_table, sep="\t", dtype=str)
    required = ["MOE_score", "n_pathogenic_likely_pathogenic", "n_plp_ge1star"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Error: master table is missing required column(s): {missing}", file=sys.stderr)
        sys.exit(1)

    df["MOE_score"] = pd.to_numeric(df["MOE_score"], errors="coerce")
    df["n_pathogenic_likely_pathogenic"] = pd.to_numeric(df["n_pathogenic_likely_pathogenic"], errors="coerce").fillna(0)
    df["n_plp_ge1star"] = pd.to_numeric(df["n_plp_ge1star"], errors="coerce").fillna(0)
    df = df.dropna(subset=["MOE_score"])

    if args.class_filter:
        if "Class" not in df.columns:
            print("Error: --class-filter given but master table has no 'Class' column", file=sys.stderr)
            sys.exit(1)
        df = df[df["Class"] == args.class_filter]
        if df.empty:
            print(f"Error: no genes found for Class == '{args.class_filter}'", file=sys.stderr)
            sys.exit(1)

    clinvar_col = "n_plp_ge1star" if args.strict else "n_pathogenic_likely_pathogenic"
    true_pos_all = (df[clinvar_col] > 0).to_numpy()
    moe = df["MOE_score"].to_numpy()

    n = len(df)
    prevalence = true_pos_all.mean()
    print(f"Genes analyzed: {n:,}{' (Class = ' + args.class_filter + ')' if args.class_filter else ''}")
    print(f"ClinVar-positive label: {clinvar_col} > 0 ({'>=1-star reviewed' if args.strict else 'any review status'})")
    print(f"ClinVar-positive prevalence: {true_pos_all.sum():,}/{n:,} ({prevalence * 100:.1f}%)")
    print()

    rows = []
    for cutoff in range(6):
        pred_pos = moe >= cutoff
        tp, fp, fn, tn = confusion_counts(pred_pos, true_pos_all)
        precision = tp / (tp + fp) if (tp + fp) else float("nan")
        recall = tp / (tp + fn) if (tp + fn) else float("nan")
        specificity = tn / (tn + fp) if (tn + fp) else float("nan")
        accuracy = (tp + tn) / n
        or_est, ci_low, ci_high = odds_ratio_ci(tp, fp, fn, tn)
        _, p_value = fisher_exact([[tp, fp], [fn, tn]])
        rows.append({
            "moe_cutoff": f">={cutoff}", "n_predicted_positive": tp + fp,
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": precision, "recall": recall, "specificity": specificity,
            "accuracy": accuracy, "odds_ratio": or_est,
            "or_ci_low": ci_low, "or_ci_high": ci_high, "fisher_p": p_value,
        })

    table = pd.DataFrame(rows)
    out_table = os.path.join(args.output_dir, "moe_clinvar_operating_characteristics.tsv")
    table.to_csv(out_table, sep="\t", index=False, float_format="%.4g")
    print(f"Operating-characteristics table saved to: {out_table}\n")
    print(table.to_string(index=False))

    # ROC and precision-recall curves. MOE_score is a 6-level (0-5) ordinal score, so the
    # natural set of thresholds is >=0 (everyone positive) through >=6 (nobody positive).
    fpr_pts, tpr_pts, prec_pts, rec_pts = [], [], [], []
    for cutoff in range(7):
        pred_pos = moe >= cutoff
        tp, fp, fn, tn = confusion_counts(pred_pos, true_pos_all)
        tpr = tp / (tp + fn) if (tp + fn) else 0.0
        fpr = fp / (fp + tn) if (fp + tn) else 0.0
        prec = tp / (tp + fp) if (tp + fp) else 1.0  # by convention at the empty-positive-set end
        fpr_pts.append(fpr)
        tpr_pts.append(tpr)
        prec_pts.append(prec)
        rec_pts.append(tpr)

    order = np.argsort(fpr_pts)
    fpr_sorted = np.array(fpr_pts)[order]
    tpr_sorted = np.array(tpr_pts)[order]
    roc_auc = np.trapezoid(tpr_sorted, fpr_sorted)

    order_pr = np.argsort(rec_pts)
    rec_sorted = np.array(rec_pts)[order_pr]
    prec_sorted = np.array(prec_pts)[order_pr]
    pr_auc = np.trapezoid(prec_sorted, rec_sorted)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    ax1.plot(fpr_pts, tpr_pts, "o-", color="#4689a3", lw=2, markersize=5)
    ax1.plot([0, 1], [0, 1], "--", color="#999999", lw=1)
    for cutoff, x, y in zip(range(7), fpr_pts, tpr_pts):
        if cutoff <= 5:
            ax1.annotate(f">={cutoff}", (x, y), textcoords="offset points", xytext=(4, -8), fontsize=8)
    ax1.set_xlabel("False positive rate")
    ax1.set_ylabel("True positive rate")
    ax1.set_title(f"ROC curve (AUC = {roc_auc:.3f})")
    ax1.set_xlim(-0.02, 1.02)
    ax1.set_ylim(-0.02, 1.02)
    for spine in ax1.spines.values():
        spine.set_visible(True)

    ax2.plot(rec_pts, prec_pts, "o-", color="#4689a3", lw=2, markersize=5)
    ax2.axhline(prevalence, ls="--", color="#999999", lw=1, label=f"baseline ({prevalence:.2f})")
    for cutoff, x, y in zip(range(7), rec_pts, prec_pts):
        if cutoff <= 5:
            ax2.annotate(f">={cutoff}", (x, y), textcoords="offset points", xytext=(4, -8), fontsize=8)
    ax2.set_xlabel("Recall")
    ax2.set_ylabel("Precision")
    ax2.set_title(f"Precision-recall curve (AUC = {pr_auc:.3f})")
    ax2.set_xlim(-0.02, 1.02)
    ax2.set_ylim(-0.02, 1.02)
    ax2.legend(loc="lower left", fontsize=8)
    for spine in ax2.spines.values():
        spine.set_visible(True)

    fig.suptitle(
        f"MOE score vs. ClinVar {'>=1-star ' if args.strict else ''}P/LP status"
        f"{' -- ' + args.class_filter if args.class_filter else ''} (n={n:,})",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_fig = os.path.join(args.output_dir, "moe_clinvar_roc_pr.png")
    fig.savefig(out_fig, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nROC/PR figure saved to: {out_fig}")
    print(f"ROC AUC: {roc_auc:.4f}  |  PR AUC: {pr_auc:.4f}")


if __name__ == "__main__":
    main()
