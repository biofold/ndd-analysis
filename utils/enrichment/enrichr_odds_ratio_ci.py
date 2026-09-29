#!/usr/bin/env python3
"""
Add a 95% confidence interval to the Odds Ratio column of an Enrichr/gseapy
enrichment output file.

WHY THIS IS NEEDED
------------------
scripts/1_enrichr_all.py's output has Term, Overlap, P-value, Adjusted
P-value, Odds Ratio, Combined Score, Genes -- an odds ratio point estimate
with no interval. The odds ratio is that of the standard 2x2 enrichment
table:

                    in term        not in term
    in gene list       a                b
    not in gene list   c                d

`Overlap` is written as "k/n": k is the number of genes from the query list
annotated to the term (cell a), n is the total number of genes annotated to
the term in the background (a + c). Given the query list size and the
background size, the whole table is determined:

    a = k
    b = list_size - k
    c = n - k
    d = background_size - list_size - (n - k)

which reproduces GSEApy's own Odds Ratio column — not Enrichr's, since the
pipeline's calls to gp.enrichr() always pass a local .gmt file plus a custom
background, routing them to GSEApy's offline enrich_local() rather than the real
Enrichr web service — once the same Haldane-Anscombe +0.5 correction gseapy itself
applies to every cell, on every row, is included (verified across all 3,352 rows of
results/main/gene_set1_GO_Biological_Process_Cancer_2026.tsv (list_size=6,856), where
the corrected formula matches the file's Odds Ratio column exactly).

The 95% CI is the standard Wald interval on the log odds ratio (mirroring
utils/moe_clinvar_validation.py's odds_ratio_ci, which does the same thing for
the MOE-threshold operating-characteristics table): SE(log OR) =
sqrt(1/a + 1/b + 1/c + 1/d),computed on the same Haldane-Anscombe-corrected cells
(+0.5 to a, b, c, d unconditionally, on every row) used for the Odds Ratio itself
above.

Usage:
  python3 utils/enrichr_odds_ratio_ci.py \
      results/main/gene_set1_GO_Biological_Process_2026.tsv \
      --list-size 6856 --background-size 19354 \
      --output results/main/gene_set1_GO_Biological_Process_2026_with_ci.tsv
"""

import argparse
import sys

import numpy as np
import pandas as pd


def odds_ratio_ci(a, b, c, d):
    """Odds ratio with 95% CI on the log scale. Haldane-Anscombe +0.5 is
    applied to every cell unconditionally, matching gseapy's own Odds Ratio
    column (gseapy.stats.calc_pvalues uses bu = 0.5 added to every cell on
    every row, not only when a cell is zero -- see GSEApy issue #132), so
    the point estimate computed here reproduces the Odds Ratio already in
    the file, and the CI is centered on that same value."""
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_est = (a * d) / (b * c)
    log_or = np.log(or_est)
    se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return np.exp(log_or - 1.96 * se_log_or), np.exp(log_or + 1.96 * se_log_or)


def _reorder_ci_columns(df):
    """Move OR_CI_low/OR_CI_high to sit immediately after 'Odds Ratio'
    (falls back to leaving them at the end if that column is absent)."""
    cols = [c for c in df.columns if c not in ("OR_CI_low", "OR_CI_high")]
    if "Odds Ratio" in cols:
        insert_at = cols.index("Odds Ratio") + 1
        cols = cols[:insert_at] + ["OR_CI_low", "OR_CI_high"] + cols[insert_at:]
    else:
        cols = cols + ["OR_CI_low", "OR_CI_high"]
    return df[cols]


def add_ci(df, list_size, background_size, on_bad_row="raise"):
    """Append OR_CI_low/OR_CI_high to an Enrichr/gseapy results dataframe.

    on_bad_row: "raise" exits with an error message naming the offending
    term(s) (used by the CLI below, where a bad --list-size/--background-size
    is a usage mistake worth stopping on). "nan" leaves OR_CI_low/OR_CI_high
    as NaN for any row whose reconstructed 2x2 table would have a negative
    cell, instead of aborting -- used when this is called inline from
    scripts/1_enrichr_all.py, where one malformed row must not take down an
    enrichment run covering many libraries and gene sets.
    """
    if len(df) == 0:
        df = df.copy()
        df["OR_CI_low"] = pd.Series(dtype=float)
        df["OR_CI_high"] = pd.Series(dtype=float)
        return _reorder_ci_columns(df)

    k = df["Overlap"].str.split("/").str[0].astype(int)
    n = df["Overlap"].str.split("/").str[1].astype(int)

    a = k
    b = list_size - k
    c = n - k
    d = background_size - list_size - c

    bad = (b < 0) | (c < 0) | (d < 0)
    if bad.any():
        bad_terms = df.loc[bad, "Term"].tolist()
        if on_bad_row == "raise":
            sys.exit(
                "Error: --list-size/--background-size do not fit these Overlap "
                f"values (negative table cell) for term(s): {bad_terms[:3]}. Check "
                "that both sizes are the ones scripts/1_enrichr_all.py was run with."
            )
        sys.stderr.write(
            f"Warning: {bad.sum()} term(s) have an Overlap inconsistent with "
            f"list_size={list_size}/background_size={background_size}; leaving "
            f"OR_CI blank for them, e.g. {bad_terms[:3]}\n"
        )

    df = df.copy()
    ci_low = pd.Series(np.nan, index=df.index, dtype=float)
    ci_high = pd.Series(np.nan, index=df.index, dtype=float)
    ok = ~bad
    for idx, ai, bi, ci_, di in zip(df.index[ok], a[ok], b[ok], c[ok], d[ok]):
        ci_low.loc[idx], ci_high.loc[idx] = odds_ratio_ci(ai, bi, ci_, di)
    df["OR_CI_low"] = ci_low
    df["OR_CI_high"] = ci_high
    return _reorder_ci_columns(df)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("enrichment_file", help="Enrichr/gseapy output TSV")
    ap.add_argument("--list-size", type=int, required=True,
                    help="Number of genes in the query gene list")
    ap.add_argument("--background-size", type=int, required=True,
                    help="Number of genes in the background/universe")
    ap.add_argument("--output", help="Output path (default: overwrite input)")
    args = ap.parse_args()

    df = pd.read_csv(args.enrichment_file, sep="\t")
    if "Overlap" not in df.columns:
        sys.exit("Error: expected an 'Overlap' column (Enrichr/gseapy output)")

    out = add_ci(df, args.list_size, args.background_size)
    out_path = args.output or args.enrichment_file
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(out)} rows with OR_CI_low/OR_CI_high to: {out_path}")
    print(out[["Term", "Overlap", "Odds Ratio", "OR_CI_low", "OR_CI_high"]].head(5).to_string(index=False))
