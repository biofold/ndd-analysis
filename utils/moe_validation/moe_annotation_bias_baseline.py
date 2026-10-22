#!/usr/bin/env python3
"""
MOE score versus GO annotation density (R1.3, R2.9).

Reviewer 2 (comment 9) asks whether the MOE score "correlates with annotation
density per gene", since a score built from GO/KEGG/Reactome enrichment may
favour well-annotated genes. Annotation density is defined here as

  n_go_annotations: the number of distinct GO terms a gene is annotated to in
  the GO Biological Process, Cellular Component and Molecular Function 2026
  libraries (libs/GO_*_2026.gmt) -- the same libraries the MOE score is
  computed from (data/gene_go_annotation_counts.tsv; rebuilt and checked
  against the .gmt files by --check-libs).

For the candidate set (MOE is defined only for candidates) the script reports
  1. GO-term count per MOE tier (n, mean, median, quartiles, genes with no
     GO term)                                   -> moe_annotation_density_by_tier.tsv
  2. Spearman correlation of MOE with the GO-term count, 95% CI from the
     Fisher z-transform                         -> moe_annotation_density_correlations.tsv
  3. MOE, GO-term count and a permuted MOE score as single predictors of each
     outcome: ROC-AUC and average precision (step-wise sum of precision over
     recall increments)                         -> moe_annotation_bias_baseline_auc.tsv
  4. Logistic regression of each outcome on MOE, unadjusted and adjusted for
     log(1 + GO-term count)                     -> moe_annotation_bias_adjusted_model.tsv
  5. ROC curves of the three predictors, one panel per outcome
                                                -> moe_annotation_bias_baseline_roc_pr.png

Outcomes: ClinVar P/LP status (from the master table) and, if --tada-input is
given (utils/moe_validation/build_tada_input.py), TADA-NDD FDR <= 0.1 from
Fu et al. 2022.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def roc_pr(score, true_pos):
    """ROC points and AUC (trapezoidal), and step-wise average precision,
    sweeping every unique score value as a '>=' threshold."""
    thresholds = np.unique(score)
    thresholds = np.concatenate([[thresholds.min() - 1], thresholds, [thresholds.max() + 1]])
    n_pos = true_pos.sum()
    n_neg = len(true_pos) - n_pos
    fpr, tpr, prec = [], [], []
    for t in thresholds:
        pred = score >= t
        tp = int(np.sum(pred & true_pos))
        fp = int(np.sum(pred & ~true_pos))
        tpr.append(tp / n_pos if n_pos else 0.0)
        fpr.append(fp / n_neg if n_neg else 0.0)
        prec.append(tp / (tp + fp) if (tp + fp) else 1.0)
    fpr, tpr, prec = map(np.array, (fpr, tpr, prec))
    order = np.lexsort((tpr, fpr))
    roc_auc = float(_trapezoid(tpr[order], fpr[order]))
    o = np.argsort(tpr, kind="stable")
    rec_s, prec_s = tpr[o], prec[o]
    ap = float(np.sum(np.diff(rec_s) * prec_s[1:]))
    return fpr[order], tpr[order], roc_auc, ap


def fmt_p(p):
    return f"{p:.4g}" if p > 0 else "<2.2e-308"


def logit_rows(y, df, outcome):
    rows = []
    specs = {
        "unadjusted": ["MOE_score"],
        "adjusted for GO-term count": ["MOE_score", "log1p_n_go_annotations"],
    }
    for model_name, cols in specs.items():
        X = sm.add_constant(df[cols].astype(float))
        fit = sm.Logit(y.astype(int), X).fit(disp=False)
        ci = fit.conf_int(alpha=0.05)
        for term in cols:
            rows.append({
                "outcome": outcome, "model": model_name, "term": term,
                "n": len(y), "n_positive": int(y.sum()),
                "coef": fit.params[term], "std_err": fit.bse[term],
                "odds_ratio": np.exp(fit.params[term]),
                "or_ci_low": np.exp(ci.loc[term, 0]), "or_ci_high": np.exp(ci.loc[term, 1]),
                "p_value": fmt_p(fit.pvalues[term]),
            })
    return rows


def check_against_libs(go_counts, lib_dir):
    from collections import Counter
    cnt = Counter()
    for lib in ("GO_Biological_Process_2026", "GO_Cellular_Component_2026", "GO_Molecular_Function_2026"):
        path = os.path.join(lib_dir, f"{lib}.gmt")
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                for g in {x for x in f[2:] if x}:
                    cnt[g] += 1
    rebuilt = go_counts["Gene"].map(cnt).fillna(0).astype(int)
    n_diff = int((rebuilt != go_counts["n_go_annotations"].astype(int)).sum())
    print(f"GO-term counts rebuilt from {lib_dir}: {len(cnt):,} genes; "
          f"mismatches with the count file: {n_diff}")
    if n_diff:
        sys.exit("Error: gene_go_annotation_counts.tsv does not match the GO 2026 libraries")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv")
    ap.add_argument("--go-counts", default="data/gene_go_annotation_counts.tsv",
                    help="Gene, n_go_annotations")
    ap.add_argument("--tada-input", default=None,
                    help="Optional output of build_tada_input.py (Gene, TADA_NDD)")
    ap.add_argument("--check-libs", default=None,
                    help="Optional libs/ directory: verify the GO-term counts against the .gmt files")
    ap.add_argument("--strict", action="store_true",
                    help="Use n_plp_ge1star instead of n_pathogenic_likely_pathogenic for ClinVar")
    ap.add_argument("--class-filter", default="candidate")
    ap.add_argument("--seed", type=int, default=0, help="Seed for the permuted-MOE null")
    ap.add_argument("--output-dir", default=".")
    args = ap.parse_args()

    for path in [args.master_table, args.go_counts] + ([args.tada_input] if args.tada_input else []):
        if not os.path.exists(path):
            print(f"Error: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    go_counts = pd.read_csv(args.go_counts, sep="\t")
    go_counts = go_counts.groupby("Gene", as_index=False)["n_go_annotations"].max()
    if args.check_libs:
        check_against_libs(go_counts, args.check_libs)

    df = pd.read_csv(args.master_table, sep="\t", dtype=str)
    df = df.drop(columns=[c for c in ("n_go_annotations",) if c in df.columns])
    df = df.merge(go_counts, on="Gene", how="left")
    df["n_go_annotations"] = df["n_go_annotations"].fillna(0).astype(int)
    df["MOE_score"] = pd.to_numeric(df["MOE_score"], errors="coerce")
    for c in ("n_pathogenic_likely_pathogenic", "n_plp_ge1star"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    df = df.dropna(subset=["MOE_score"])
    if args.class_filter:
        df = df[df["Class"] == args.class_filter]
    df["log1p_n_go_annotations"] = np.log1p(df["n_go_annotations"].astype(float))
    n = len(df)
    print(f"Genes analyzed: {n:,} (Class = {args.class_filter})")

    # 1. GO-term count per MOE tier
    g = df.groupby(df["MOE_score"].astype(int))["n_go_annotations"]
    tier = pd.DataFrame({
        "n_genes": g.size(), "mean": g.mean(), "median": g.median(),
        "q1": g.quantile(0.25), "q3": g.quantile(0.75),
        "n_no_go_term": g.apply(lambda s: int((s == 0).sum())),
    })
    tier["fraction_no_go_term"] = tier["n_no_go_term"] / tier["n_genes"]
    tier.index.name = "MOE_score"
    out_tier = os.path.join(args.output_dir, "moe_annotation_density_by_tier.tsv")
    tier.reset_index().to_csv(out_tier, sep="\t", index=False, float_format="%.4g")
    print("\nGO-term count by MOE tier:")
    print(tier.to_string())

    # 2. Spearman correlation
    rho, p_rho = spearmanr(df["MOE_score"], df["n_go_annotations"])
    se = 1.0 / np.sqrt(n - 3)
    z = np.arctanh(rho)
    corr = pd.DataFrame([{
        "variable": "n_go_annotations", "n": n, "spearman_rho": rho,
        "rho_ci_low": np.tanh(z - 1.96 * se), "rho_ci_high": np.tanh(z + 1.96 * se),
        "p_value": fmt_p(p_rho)}])
    out_corr = os.path.join(args.output_dir, "moe_annotation_density_correlations.tsv")
    corr.to_csv(out_corr, sep="\t", index=False, float_format="%.4g")
    print("\nSpearman correlation of MOE with GO-term count:")
    print(corr.to_string(index=False))

    # Outcomes
    clinvar_col = "n_plp_ge1star" if args.strict else "n_pathogenic_likely_pathogenic"
    outcomes = {"ClinVar P/LP": df.assign(y=(df[clinvar_col] > 0).astype(int))}
    if args.tada_input:
        tada = pd.read_csv(args.tada_input, sep="\t")[["Gene", "TADA_NDD"]]
        outcomes["TADA-NDD FDR <= 0.1"] = df.merge(tada, on="Gene", how="inner").rename(
            columns={"TADA_NDD": "y"})

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(args.seed)
    colors = {"MOE_score": "#4689a3", "n_go_annotations": "#c05746",
              "permuted MOE score": "#999999"}
    fig, axes = plt.subplots(1, len(outcomes), figsize=(6 * len(outcomes), 5.5), squeeze=False)
    auc_rows, model_rows = [], []
    for ax, (name, d) in zip(axes[0], outcomes.items()):
        y = d["y"].to_numpy().astype(bool)
        preds = {"MOE_score": d["MOE_score"].to_numpy(),
                 "n_go_annotations": d["n_go_annotations"].to_numpy(),
                 "permuted MOE score": rng.permutation(d["MOE_score"].to_numpy())}
        for pname, s in preds.items():
            fpr, tpr, auc, apr = roc_pr(s, y)
            auc_rows.append({"outcome": name, "n": len(y), "n_positive": int(y.sum()),
                             "predictor": pname, "roc_auc": auc, "average_precision": apr})
            ax.plot(fpr, tpr, "-" if pname == "MOE_score" else "--", color=colors[pname],
                    lw=2 if pname == "MOE_score" else 1.3, label=f"{pname} (AUC = {auc:.3f})")
        ax.plot([0, 1], [0, 1], ":", color="#cccccc", lw=1)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title(f"{name} (n = {len(y):,}; {int(y.sum()):,} positive)")
        ax.legend(loc="lower right", fontsize=8)
        model_rows += logit_rows(d["y"].to_numpy(), d, name)
    fig.tight_layout()
    out_fig = os.path.join(args.output_dir, "moe_annotation_bias_baseline_roc_pr.png")
    fig.savefig(out_fig, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    auc_table = pd.DataFrame(auc_rows)
    out_auc = os.path.join(args.output_dir, "moe_annotation_bias_baseline_auc.tsv")
    auc_table.to_csv(out_auc, sep="\t", index=False, float_format="%.4g")
    print("\nSingle-predictor comparison:")
    print(auc_table.to_string(index=False))

    model_table = pd.DataFrame(model_rows)
    out_model = os.path.join(args.output_dir, "moe_annotation_bias_adjusted_model.tsv")
    model_table.to_csv(out_model, sep="\t", index=False, float_format="%.4g")
    print("\nLogistic regression (odds ratio per MOE point; per unit of log(1 + GO terms)):")
    print(model_table[["outcome", "model", "term", "odds_ratio", "or_ci_low", "or_ci_high",
                       "p_value"]].to_string(index=False))
    for path in (out_tier, out_corr, out_auc, out_model, out_fig):
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
