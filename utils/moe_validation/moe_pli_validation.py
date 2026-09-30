#!/usr/bin/env python3
"""
Test whether the MOE score's association with ClinVar pathogenic/likely
pathogenic (P/LP) variant status holds independently of gnomAD pLI (the
probability of loss-of-function intolerance) -- a robustness check for
reviewer request R1.3's annotation-bias concern ("may favor genes that
are simply better annotated").

gnomAD pLI is derived purely from population allele-frequency data across
~800,000 individuals and carries no dependence on disease-literature
curation, GO/pathway annotation density, or publication count. If MOE
score remains a significant predictor of ClinVar P/LP status after
adjusting for pLI in a joint logistic regression, that is evidence the
MOE score is not simply tracking how well-studied a gene is.

Fits: P(ClinVar P/LP positive) ~ MOE_score + gnomAD_pLI  (logistic
regression via statsmodels), reporting coefficients, odds ratios, 95% CI
and p-values for both terms, plus the Spearman correlation between
MOE_score and gnomAD_pLI, and per-MOE-tier pLI summary statistics.

Restricted to the candidate set by default (--all-classes to override):
including high-confidence genes would be circular, since the criteria used
to call a gene high-confidence (GeneTrek-HC, Orphanet-both, SFARI+,
Sanchis-Juan) already independently correlate with ClinVar evidence.

Input: the merged gene master table (ndd_master_table.tsv), which must
have MOE_score, Class, gnomAD_pLI, n_pathogenic_likely_pathogenic and
n_plp_ge1star columns.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.api as sm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv")
    ap.add_argument("--strict", action="store_true",
                     help="Use n_plp_ge1star (>=1-star reviewed P/LP variants) as the "
                          "ClinVar positive label instead of n_pathogenic_likely_pathogenic")
    ap.add_argument("--class-filter", default="candidate",
                     help="Restrict to one Class value (curated, candidate, no_evidence). "
                          "Default: candidate (see module docstring for why).")
    ap.add_argument("--all-classes", action="store_true",
                     help="Override --class-filter and run on all genes regardless of Class.")
    ap.add_argument("--output-dir", default=".", help="Directory for output tables and figure")
    ap.add_argument("--metric", choices=["pLI", "LOEUF"], default="pLI",
                     help="gnomAD v2.1.1 constraint metric (default: pLI). LOEUF is the "
                          "upper bound of the o/e LoF 90%% CI; LOWER LOEUF = MORE constrained, "
                          "the opposite direction to pLI.")
    args = ap.parse_args()
    if args.all_classes:
        args.class_filter = None
    col = f"gnomAD_{args.metric}"
    tag = args.metric.lower()
    label = f"gnomAD v2.1.1 {args.metric}"

    if not os.path.exists(args.master_table):
        print(f"Error: File not found: {args.master_table}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.master_table, sep="\t", dtype=str)
    required = ["MOE_score", col, "n_pathogenic_likely_pathogenic", "n_plp_ge1star"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"Error: master table is missing required column(s): {missing}", file=sys.stderr)
        sys.exit(1)

    df["MOE_score"] = pd.to_numeric(df["MOE_score"], errors="coerce")
    df[col] = pd.to_numeric(df[col], errors="coerce")
    df["n_pathogenic_likely_pathogenic"] = pd.to_numeric(df["n_pathogenic_likely_pathogenic"], errors="coerce").fillna(0)
    df["n_plp_ge1star"] = pd.to_numeric(df["n_plp_ge1star"], errors="coerce").fillna(0)

    if args.class_filter:
        if "Class" not in df.columns:
            print("Error: --class-filter given but master table has no 'Class' column", file=sys.stderr)
            sys.exit(1)
        df = df[df["Class"] == args.class_filter]

    n_total = len(df)
    df_complete = df.dropna(subset=["MOE_score", col])
    n_dropped = n_total - len(df_complete)
    print(f"Genes analyzed: {len(df_complete):,} (dropped {n_dropped:,} with missing {col}"
          f"{' / MOE_score' if df['MOE_score'].isna().any() else ''}, out of {n_total:,} total"
          f"{' in Class == ' + args.class_filter if args.class_filter else ''})\n")

    clinvar_col = "n_plp_ge1star" if args.strict else "n_pathogenic_likely_pathogenic"
    y = (df_complete[clinvar_col] > 0).astype(int)
    X = df_complete[["MOE_score", col]].astype(float)
    X = sm.add_constant(X)

    model = sm.Logit(y, X).fit(disp=False)

    coef_table = pd.DataFrame({
        "term": model.params.index,
        "coef": model.params.values,
        "std_err": model.bse.values,
        "z": model.tvalues.values,
        "p_value": model.pvalues.values,
    })
    coef_table["odds_ratio"] = np.exp(coef_table["coef"])
    ci = model.conf_int(alpha=0.05)
    coef_table["or_ci_low"] = np.exp(ci[0].values)
    coef_table["or_ci_high"] = np.exp(ci[1].values)

    out_table = os.path.join(args.output_dir, f"moe_{tag}_logistic_regression.tsv")
    coef_table.to_csv(out_table, sep="\t", index=False, float_format="%.4g")
    print(f"Joint logistic regression: P({clinvar_col} > 0) ~ MOE_score + {col}")
    print(coef_table.to_string(index=False))
    print(f"\nModel table saved to: {out_table}")

    out_summary = os.path.join(args.output_dir, f"moe_{tag}_logistic_regression_summary.txt")
    with open(out_summary, "w") as f:
        f.write(str(model.summary()))
    print(f"Full statsmodels summary saved to: {out_summary}\n")

    moe_significant = coef_table.loc[coef_table["term"] == "MOE_score", "p_value"].iloc[0] < 0.05
    print(f"MOE_score remains significant after adjusting for {col}: {moe_significant} "
          f"(p = {coef_table.loc[coef_table['term'] == 'MOE_score', 'p_value'].iloc[0]:.3g})")

    # Diagnostic: does MOE score itself correlate with pLI? (If strongly, that alone would
    # be a confound worth flagging; a weak/moderate correlation with MOE staying significant
    # above is the reassuring combination.)
    rho, rho_p = spearmanr(df_complete["MOE_score"], df_complete[col])
    print(f"\nSpearman correlation MOE_score vs {col}: rho = {rho:.3f} (p = {rho_p:.3g})")

    tier_stats = (
        df_complete.groupby("MOE_score")[col]
        .agg(n="size", mean="mean", median="median", std="std")
        .reset_index()
    )
    out_tiers = os.path.join(args.output_dir, f"moe_{tag}_by_tier.tsv")
    tier_stats.to_csv(out_tiers, sep="\t", index=False, float_format="%.4g")
    print(f"\n{col} by MOE score tier:\n{tier_stats.to_string(index=False)}")
    print(f"\nPer-tier summary saved to: {out_tiers}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    tiers = sorted(df_complete["MOE_score"].unique())
    data = [df_complete.loc[df_complete["MOE_score"] == t, col].values for t in tiers]
    bp = ax.boxplot(data, positions=tiers, widths=0.6, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor("#4689a3")
        patch.set_alpha(0.7)
    for median in bp["medians"]:
        median.set_color("#333333")
    ax.set_xlabel("MOE score")
    ax.set_ylabel(label)
    ax.set_title(
        f"{label} by MOE score tier (Spearman rho = {rho:.2f}, p = {rho_p:.2g})"
        f"{chr(10) + args.class_filter if args.class_filter else ''}",
        fontsize=11,
    )
    ax.set_xticks(tiers)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
    fig.tight_layout()

    out_fig = os.path.join(args.output_dir, f"moe_{tag}_by_tier_boxplot.png")
    fig.savefig(out_fig, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nBoxplot figure saved to: {out_fig}")


if __name__ == "__main__":
    main()
