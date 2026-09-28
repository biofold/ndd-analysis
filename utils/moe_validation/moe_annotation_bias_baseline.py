#!/usr/bin/env python3
"""
Test the reviewer's literal concern (R1.3, R2.9): "[MOE] may favor genes
that are simply better annotated." Compares the MOE score against three
baseline predictors of the same outcome (ClinVar P/LP status) that
measure how much attention a gene has received in the literature/GO
curation rather than any specific NDD evidence:

  - GO annotation count  (n_go_annotations: number of distinct GO terms,
    across BP/CC/MF, a gene is annotated to -- from the same GO libraries
    MOE's own enrichment is computed against)
  - PubMed citation count (n_pubmed: from NCBI's gene2pubmed)
  - random ranking (a permuted copy of MOE_score, as a null baseline)

For each predictor, ROC-AUC / PR-AUC is computed by sweeping the raw
predictor's own value as the classification threshold directly (AUC is
invariant to monotonic score transforms, so this gives the same AUC a
single-variable logistic regression on that predictor would) on the same
candidate-only outcome MOE is validated against elsewhere (see
moe_clinvar_validation.py), so ROC curves can be compared directly on one
plot: if MOE's curve does not clearly exceed the annotation-count and
PubMed-count curves, the reviewer's concern stands.

Then fits the WS4.7 adjusted model:
  ClinVar_P/LP ~ MOE_score + log1p(n_go_annotations) + log1p(n_pubmed)
and reports the MOE coefficient/OR with 95% CI -- if MOE remains
significant here, its association with ClinVar evidence is not merely a
proxy for annotation/study volume.

Restricted to the candidate set by default (--all-classes to override),
for the same non-circularity reason as moe_clinvar_validation.py /
moe_pli_validation.py.

Inputs: the merged gene master table (ndd_master_table.tsv), plus
data/gene_go_annotation_counts.tsv and data/gene_pubmed_counts.tsv
(built from the existing GO_*_2026.gmt libraries and NCBI's
gene2pubmed.gz + Homo_sapiens.gene_info.gz respectively).
"""

import argparse
import os
import sys

import numpy as np

# numpy.trapezoid is the renamed numpy.trapz, available only from NumPy 2.0
# onward; numpy.trapz itself is deprecated in 2.x but still present, and is
# the only name available on older NumPy. Resolve whichever exists so this
# runs on either.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr


def roc_pr_points(score, true_pos):
    """ROC/PR points by sweeping unique score values as thresholds
    (works for both continuous and discrete/ordinal scores)."""
    thresholds = np.unique(score)
    thresholds = np.concatenate([[thresholds.min() - 1], thresholds, [thresholds.max() + 1]])
    fpr, tpr, prec, rec = [], [], [], []
    n_pos = true_pos.sum()
    n_neg = len(true_pos) - n_pos
    for t in thresholds:
        pred_pos = score >= t
        tp = int(np.sum(pred_pos & true_pos))
        fp = int(np.sum(pred_pos & ~true_pos))
        fn = int(np.sum(~pred_pos & true_pos))
        tpr.append(tp / n_pos if n_pos else 0.0)
        fpr.append(fp / n_neg if n_neg else 0.0)
        prec.append(tp / (tp + fp) if (tp + fp) else 1.0)
        rec.append(tp / n_pos if n_pos else 0.0)
    order = np.argsort(fpr)
    roc_auc = _trapezoid(np.array(tpr)[order], np.array(fpr)[order])
    order_pr = np.argsort(rec)
    pr_auc = _trapezoid(np.array(prec)[order_pr], np.array(rec)[order_pr])
    return (np.array(fpr)[order], np.array(tpr)[order], roc_auc), (np.array(rec)[order_pr], np.array(prec)[order_pr], pr_auc)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv")
    ap.add_argument("--go-counts", default="data/gene_go_annotation_counts.tsv",
                     help="Path to gene_go_annotation_counts.tsv (Gene, n_go_annotations)")
    ap.add_argument("--pubmed-counts", default="data/gene_pubmed_counts.tsv",
                     help="Path to gene_pubmed_counts.tsv (Symbol, n_pubmed)")
    ap.add_argument("--strict", action="store_true",
                     help="Use n_plp_ge1star instead of n_pathogenic_likely_pathogenic as the "
                          "ClinVar positive label")
    ap.add_argument("--class-filter", default="candidate",
                     help="Restrict to one Class value. Default: candidate (see docstring).")
    ap.add_argument("--all-classes", action="store_true", help="Override --class-filter")
    ap.add_argument("--seed", type=int, default=0, help="Random seed for the permuted-MOE null baseline")
    ap.add_argument("--output-dir", default=".", help="Directory for output tables and figure")
    args = ap.parse_args()
    if args.all_classes:
        args.class_filter = None

    for path in [args.master_table, args.go_counts, args.pubmed_counts]:
        if not os.path.exists(path):
            print(f"Error: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    df = pd.read_csv(args.master_table, sep="\t", dtype=str)
    required = ["Gene", "MOE_score", "n_pathogenic_likely_pathogenic", "n_plp_ge1star"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Error: master table is missing required column(s): {missing}", file=sys.stderr)
        sys.exit(1)

    go_counts = pd.read_csv(args.go_counts, sep="\t").rename(columns={"n_go_annotations": "n_go_annotations"})
    pubmed_counts = pd.read_csv(args.pubmed_counts, sep="\t")[["Symbol", "n_pubmed"]].rename(columns={"Symbol": "Gene"})

    # NCBI's gene2pubmed keys on GeneID, so a handful of symbols (e.g. TEC,
    # MMD2) appear more than once. Left-merging those unaggregated would
    # duplicate candidate genes and inflate n relative to the other validation
    # scripts. Collapse to one row per symbol, keeping the largest count --
    # the most favourable value for the annotation-volume baseline, so the
    # comparison against MOE is not made artificially easy to win.
    go_counts = go_counts.groupby("Gene", as_index=False).max(numeric_only=True)
    pubmed_counts = pubmed_counts.groupby("Gene", as_index=False).max(numeric_only=True)

    # A master table built by utils/build_master_table.py already carries
    # n_go_annotations / n_pubmed columns. Drop them before merging so pandas
    # does not disambiguate the duplicates into _x/_y suffixes; the standalone
    # count files passed on the command line are taken as authoritative.
    df = df.drop(columns=[c for c in ("n_go_annotations", "n_pubmed") if c in df.columns])

    df = df.merge(go_counts, on="Gene", how="left")
    df = df.merge(pubmed_counts, on="Gene", how="left")
    df["n_go_annotations"] = df["n_go_annotations"].fillna(0)
    df["n_pubmed"] = df["n_pubmed"].fillna(0)

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
    true_pos = (df[clinvar_col] > 0).to_numpy()
    n = len(df)
    print(f"Genes analyzed: {n:,}{' (Class = ' + args.class_filter + ')' if args.class_filter else ''}")
    print(f"ClinVar-positive prevalence: {true_pos.sum():,}/{n:,} ({true_pos.mean() * 100:.1f}%)\n")

    # Reviewer R2.9 asks literally whether the MOE score correlates with
    # annotation density per gene. Report that directly (Spearman, since MOE is
    # a 0-5 ordinal and the count variables are heavily right-skewed) before
    # moving on to the harder question of whether MOE beats those counts as a
    # predictor. A 95% CI is obtained from the Fisher z-transform of rho.
    corr_rows = []
    for label in ("n_go_annotations", "n_pubmed"):
        rho, p_rho = spearmanr(df["MOE_score"], df[label])
        se = 1.0 / np.sqrt(n - 3)
        z = np.arctanh(rho)
        lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
        # scipy returns exactly 0.0 when the p-value underflows double
        # precision; reporting "p = 0" in a supplementary table is misleading,
        # so record the representable bound instead.
        p_display = f"{p_rho:.4g}" if p_rho > 0 else "<2.2e-308"
        corr_rows.append({"variable": label, "n": n, "spearman_rho": rho,
                          "rho_ci_low": lo, "rho_ci_high": hi, "p_value": p_display})
    corr_table = pd.DataFrame(corr_rows)
    print("Spearman correlation of MOE score with annotation density:")
    print(corr_table.to_string(index=False))
    out_corr = os.path.join(args.output_dir, "moe_annotation_density_correlations.tsv")
    corr_table.to_csv(out_corr, sep="\t", index=False, float_format="%.4g")
    print(f"\nCorrelation table saved to: {out_corr}\n")

    rng = np.random.default_rng(args.seed)
    permuted_moe = rng.permutation(df["MOE_score"].to_numpy())

    predictors = {
        "MOE_score": df["MOE_score"].to_numpy(),
        "n_go_annotations": df["n_go_annotations"].to_numpy(),
        "n_pubmed": df["n_pubmed"].to_numpy(),
        "random_ranking (permuted MOE)": permuted_moe,
    }

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
    colors = {"MOE_score": "#4689a3", "n_go_annotations": "#c05746", "n_pubmed": "#e0a03c",
              "random_ranking (permuted MOE)": "#999999"}
    summary_rows = []
    for name, score in predictors.items():
        (fpr, tpr, roc_auc), (rec, prec, pr_auc) = roc_pr_points(score, true_pos)
        style = "-" if name == "MOE_score" else "--"
        ax1.plot(fpr, tpr, style, color=colors[name], lw=2 if name == "MOE_score" else 1.3,
                 label=f"{name} (AUC={roc_auc:.3f})")
        ax2.plot(rec, prec, style, color=colors[name], lw=2 if name == "MOE_score" else 1.3,
                 label=f"{name} (AUC={pr_auc:.3f})")
        summary_rows.append({"predictor": name, "roc_auc": roc_auc, "pr_auc": pr_auc})

    ax1.plot([0, 1], [0, 1], ":", color="#cccccc", lw=1)
    ax1.set_xlabel("False positive rate")
    ax1.set_ylabel("True positive rate")
    ax1.set_title("ROC: MOE vs. annotation-volume baselines")
    ax1.legend(loc="lower right", fontsize=8)
    for spine in ax1.spines.values():
        spine.set_visible(True)

    ax2.axhline(true_pos.mean(), ls=":", color="#cccccc", lw=1)
    ax2.set_xlabel("Recall")
    ax2.set_ylabel("Precision")
    ax2.set_title("Precision-recall: MOE vs. annotation-volume baselines")
    ax2.legend(loc="upper right", fontsize=8)
    for spine in ax2.spines.values():
        spine.set_visible(True)

    fig.suptitle(
        f"Does MOE beat simple annotation-volume proxies?"
        f"{' -- ' + args.class_filter if args.class_filter else ''} (n={n:,})",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out_fig = os.path.join(args.output_dir, "moe_annotation_bias_baseline_roc_pr.png")
    fig.savefig(out_fig, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    summary = pd.DataFrame(summary_rows)
    out_summary = os.path.join(args.output_dir, "moe_annotation_bias_baseline_auc.tsv")
    summary.to_csv(out_summary, sep="\t", index=False, float_format="%.4g")
    print("Single-predictor AUC comparison:")
    print(summary.to_string(index=False))
    print(f"\nBaseline-comparison figure saved to: {out_fig}")
    print(f"AUC summary saved to: {out_summary}\n")

    moe_beats_go = summary.loc[summary["predictor"] == "MOE_score", "roc_auc"].iloc[0] > \
                   summary.loc[summary["predictor"] == "n_go_annotations", "roc_auc"].iloc[0]
    moe_beats_pubmed = summary.loc[summary["predictor"] == "MOE_score", "roc_auc"].iloc[0] > \
                       summary.loc[summary["predictor"] == "n_pubmed", "roc_auc"].iloc[0]
    print(f"MOE_score ROC-AUC exceeds n_go_annotations: {moe_beats_go}")
    print(f"MOE_score ROC-AUC exceeds n_pubmed: {moe_beats_pubmed}\n")

    # WS4.7 adjusted model
    X = pd.DataFrame({
        "MOE_score": df["MOE_score"].astype(float),
        "log1p_n_go_annotations": np.log1p(df["n_go_annotations"].astype(float)),
        "log1p_n_pubmed": np.log1p(df["n_pubmed"].astype(float)),
    })
    X = sm.add_constant(X)
    model = sm.Logit(true_pos.astype(int), X).fit(disp=False)

    coef_table = pd.DataFrame({
        "term": model.params.index, "coef": model.params.values,
        "std_err": model.bse.values, "z": model.tvalues.values, "p_value": model.pvalues.values,
    })
    coef_table["odds_ratio"] = np.exp(coef_table["coef"])
    ci = model.conf_int(alpha=0.05)
    coef_table["or_ci_low"] = np.exp(ci[0].values)
    coef_table["or_ci_high"] = np.exp(ci[1].values)

    out_model = os.path.join(args.output_dir, "moe_annotation_bias_adjusted_model.tsv")
    coef_table.to_csv(out_model, sep="\t", index=False, float_format="%.4g")
    print("Adjusted model: P(ClinVar P/LP) ~ MOE_score + log1p(n_go_annotations) + log1p(n_pubmed)")
    print(coef_table.to_string(index=False))
    print(f"\nAdjusted model table saved to: {out_model}")

    moe_p = coef_table.loc[coef_table["term"] == "MOE_score", "p_value"].iloc[0]
    print(f"\nMOE_score remains significant after adjusting for annotation/PubMed volume: "
          f"{moe_p < 0.05} (p = {moe_p:.3g})")


if __name__ == "__main__":
    main()
