#!/usr/bin/env python3
"""
Cross-validate the gnomAD_pLI values pulled from dbNSFP
(extract_gnomad_pli.py's output, ultimately a v2.1.1 value) against
gnomAD's own primary gene constraint release -- by default the v4.1
constraint table at release/4.1/constraint/ on gnomAD's AWS bucket.

NOTE ON v4.1 vs v4.1.1: gnomAD's bucket also has a genuinely separate
release/4.1.1/constraint/gnomad.v4.1.1.constraint_metrics.ht/ object
(confirmed present via a bucket listing, dated 2026-03-25 -- NOT the
same object as the v4.1 file used here, and NOT simply a same-path
in-place patch as an earlier version of this docstring incorrectly
claimed). That object is in Hail-native Table format (.ht/), not a
flat TSV, so it cannot be read with plain pandas.read_csv -- reading it
requires the `hail` package. This script does NOT currently fetch the
true v4.1.1 file; it uses v4.1. If you need the actual v4.1.1 values,
install hail and read the .ht table, or ask for that to be added here.

Why this exists: dbNSFP is a secondary aggregator that repackages many
annotation sources, including gnomAD's constraint metrics, and its
gnomAD_pLI column is specifically the older v2.1.1 value. This script
answers the reviewer-relevant question "how does that compare to
gnomAD's own current, primary constraint table?" by fetching it
directly and comparing per-gene pLI values.

IMPORTANT -- v4.1 pLI is NOT expected to closely match v2.1.1 pLI.
gnomAD v4 uses a ~6x larger cohort (807,162 vs 141,456 individuals) and
its own release notes state pLI/missense-Z were recalculated on the
new dataset, with LOEUF now recommended as the *primary* constraint
metric (pLI kept only as a secondary, backward-compatible <0.9 cutoff).
So a low exact-match rate here reflects a genuine update to the
underlying constraint estimate, not an extraction bug -- Spearman rho
between the two versions is reported to show they still rank genes
similarly even where the point estimate moved. Because
ndd_master_table.tsv's gnomAD_pLI column is the v2.1.1 value (via
dbNSFP), using this v4.1 comparison as a drop-in replacement would
change results throughout the pipeline, not just add a robustness
check -- treat this script's output as informational unless you
decide to re-baseline the whole master table on v4.1.

Primary source: gnomAD's own AWS Open Data bucket (ARN
arn:aws:s3:::gnomad-public-us-east-1), object
release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv (95MB,
uncompressed; one row per gene x transcript, both RefSeq- and
Ensembl-named rows for the same transcript). Restrict to
mane_select == True to get one transcript per gene (MANE Select,
introduced in v4.1 specifically to remove the canonical-transcript
ambiguity that v2.1.1's file had); duplicate MANE-select rows for the
same gene (RefSeq NM_ vs Ensembl ENST_ naming of the same transcript)
were checked to always agree on pLI, so are safely deduplicated.
Pass --v2.1.1-url / --v2.1.1-input to instead fetch/use the older
v2.1.1 file for a same-version comparison (see
release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz).
"""

import argparse
import os
import sys
import urllib.request

import pandas as pd

GNOMAD_V41_URL = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
GNOMAD_V211_URL = "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"


def download(url, dest):
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"✓ Downloaded to {dest} ({os.path.getsize(dest) / 1e6:.1f} MB)")


def load_v41(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[df["mane_select"] == True][["gene", "transcript", "lof.pLI"]]
    df = df.rename(columns={"gene": "Gene", "lof.pLI": "pLI_gnomad_direct"})
    return df


def load_v211(path):
    df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    df = df[["gene", "transcript", "pLI"]].rename(
        columns={"gene": "Gene", "pLI": "pLI_gnomad_direct"}
    )
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", choices=["4.1", "2.1.1"], default="4.1",
                     help="Which gnomAD constraint release to compare against (default: 4.1)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", help="Path to an already-downloaded constraint file for --version")
    src.add_argument("--url", nargs="?", const="__default__",
                      help="Download the file for --version (default if flag given with no value: "
                           "gnomAD's own URL for that version)")
    ap.add_argument("--download-to", default=None,
                     help="Where to save the file when --url is used (default: version-specific filename)")
    ap.add_argument("--dbnsfp-pli-table", required=True,
                     help="Path to extract_gnomad_pli.py's output TSV (Gene, gnomAD_pLI, ...)")
    ap.add_argument("--output", default="gnomad_pli_source_comparison.tsv",
                     help="Where to write the per-gene comparison table")
    args = ap.parse_args()

    default_url = GNOMAD_V41_URL if args.version == "4.1" else GNOMAD_V211_URL
    default_name = "gnomad.v4.1.constraint_metrics.tsv" if args.version == "4.1" else "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

    if args.url:
        url = default_url if args.url == "__default__" else args.url
        dest = args.download_to or default_name
        download(url, dest)
        gnomad_path = dest
    else:
        gnomad_path = args.input
        if not os.path.exists(gnomad_path):
            print(f"Error: File not found: {gnomad_path}", file=sys.stderr)
            sys.exit(1)

    print(f"Reading {gnomad_path} (gnomAD v{args.version}) ...")
    gnomad = load_v41(gnomad_path) if args.version == "4.1" else load_v211(gnomad_path)

    print(f"Reading {args.dbnsfp_pli_table} ...")
    dbnsfp = pd.read_csv(args.dbnsfp_pli_table, sep="\t")[["Gene", "gnomAD_pLI"]]

    n_dup_genes = gnomad["Gene"].duplicated().sum()
    print(f"\ngnomAD v{args.version} file: {len(gnomad)} rows, {gnomad['Gene'].nunique()} unique gene symbols "
          f"({n_dup_genes} extra rows per already-seen gene symbol)")

    merged_all = dbnsfp.merge(gnomad.drop_duplicates(subset="Gene", keep="first"), on="Gene", how="inner")
    # "comparable" means both values are actually present -- an inner join on Gene alone
    # does NOT guarantee that, since either source table can carry a NaN value for a gene
    # it still lists a row for (e.g. dbNSFP genes with no dbNSFP-side pLI, or gnomAD genes
    # with an incalculable constraint estimate). Filtering this explicitly, rather than
    # letting NaN silently fall out of the "< 1e-6" exact-match check and get miscounted
    # as a mismatch, is required for n_mismatch below to mean anything.
    n_present_only_one_side = merged_all[["gnomAD_pLI", "pLI_gnomad_direct"]].isna().any(axis=1).sum()
    merged = merged_all.dropna(subset=["gnomAD_pLI", "pLI_gnomad_direct"]).copy()
    merged["abs_diff"] = (merged["gnomAD_pLI"] - merged["pLI_gnomad_direct"]).abs()

    n_total = len(merged)
    n_exact = (merged["abs_diff"] < 1e-6).sum()
    n_mismatch = n_total - n_exact

    print(f"\nGenes matched by symbol: {len(merged_all)} "
          f"({n_present_only_one_side} of those have a missing value on at least one side, excluded below)")
    print(f"Genes comparable (non-missing in both): {n_total}")
    print(f"Exact match (|diff| < 1e-6): {n_exact} ({n_exact / n_total:.2%})")
    print(f"Mismatch: {n_mismatch} ({n_mismatch / n_total:.2%})")

    if args.version == "4.1":
        from scipy.stats import spearmanr
        clean = merged.dropna(subset=["gnomAD_pLI", "pLI_gnomad_direct"])
        rho, p = spearmanr(clean["gnomAD_pLI"], clean["pLI_gnomad_direct"])
        print(f"\nSpearman rho (dbNSFP v2.1.1 vs gnomAD v4.1 MANE-select pLI, n={len(clean)}): {rho:.4f} (p={p:.2e})")
        print("A low exact-match rate against v4.1 is EXPECTED (recalculated on a ~6x larger cohort) "
              "-- see docstring. This does not indicate a problem with the v2.1.1 extraction.")
    else:
        mismatches = merged[merged["abs_diff"] >= 1e-6].sort_values("abs_diff", ascending=False)
        if len(mismatches):
            multi_transcript_genes = set(gnomad[gnomad["Gene"].duplicated(keep=False)]["Gene"])
            mismatches = mismatches.copy()
            mismatches["gene_has_multiple_transcripts_in_gnomad_file"] = mismatches["Gene"].isin(multi_transcript_genes)
            n_explained = mismatches["gene_has_multiple_transcripts_in_gnomad_file"].sum()
            print(f"\nOf the {n_mismatch} mismatches, {n_explained} are genes with more than one "
                  f"transcript row in gnomAD's own file (ambiguous canonical-transcript choice, "
                  f"not a real discrepancy); {n_mismatch - n_explained} are unexplained.")
            print("\nTop mismatches:")
            print(mismatches.head(10).to_string(index=False))

    merged.to_csv(args.output, sep="\t", index=False)
    print(f"\n✓ Full comparison table saved to {args.output}")


if __name__ == "__main__":
    main()
