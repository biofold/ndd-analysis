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

which reproduces Enrichr's own Odds Ratio = a*d / (b*c) (matches to within
floating-point/universe-filtering differences -- verified against
results/main/gene_set1_GO_Biological_Process_2026.tsv).

The 95% CI is the standard Wald interval on the log odds ratio (mirroring
utils/moe_clinvar_validation.py's odds_ratio_ci, which does the same thing for
the MOE-threshold operating-characteristics table): SE(log OR) =
sqrt(1/a + 1/b + 1/c + 1/d), with a Haldane-Anscombe +0.5 correction to every
cell when any cell is zero (only possible here if a term's Overlap equals the
whole query list or the whole background, which does not happen in practice
but is guarded against).

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
    """Odds ratio with 95% CI on the log scale; Haldane-Anscombe +0.5
    correction applied uniformly when any cell is zero, to avoid an
    undefined or infinite OR. Mirrors moe_clinvar_validation.py's
    odds_ratio_ci for the same reason."""
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    or_est = (a * d) / (b * c)
    log_or = np.log(or_est)
    se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return np.exp(log_or - 1.96 * se_log_or), np.exp(log_or + 1.96 * se_log_or)


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
        return df

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
    return df


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
