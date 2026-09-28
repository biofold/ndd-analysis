#!/usr/bin/env python3
"""
Extract Pathogenic/Likely-pathogenic (P/LP) variant records from the raw
NCBI ClinVar variant_summary release and aggregate them to a per-gene
P/LP variant count.

This is the actual source of clinvar_plp_variants.tsv and of the
n_pathogenic_likely_pathogenic / n_plp_ge1star columns in
ndd_master_table.tsv -- previously this step was done interactively and
was not reproducible from the repo. Making it a first-class pipeline
script closes that gap (used for R1.3's independent-validation
analysis, moe_clinvar_validation.py / bar.py / fisher-test.py).

P/LP classification rule (important -- do NOT use ClinSigSimple alone):
  A row counts as P/LP if its ClinicalSignificance TEXT contains
  "pathogenic" (case-insensitive) AND does NOT also contain
  "conflicting" or "uncertain". ClinVar's own ClinSigSimple flag is not
  used as the sole criterion because a single VariationID/AlleleID can
  be linked to multiple conditions with different classifications, and
  ClinSigSimple=1 can still appear on a row whose displayed text reads
  "Uncertain significance" for the condition actually reported on that
  row. Matching manuscript/pipeline convention, this is done at the
  row level (one row per VariationID x Assembly) -- per-condition
  PhenotypeIDs/PhenotypeList breakdowns are not parsed.

Only GRCh38 rows are kept (ClinVar reports both GRCh37 and GRCh38
coordinates for most variants; keeping both would double-count).

A variant is additionally flagged "ge1star" (>=1-star ClinVar review
status) when its ReviewStatus is anything other than ClinVar's
zero-star categories ("no assertion criteria provided", "no assertion
provided", "no classification provided", "no classification for the
individual variant").

Input: NCBI ClinVar's variant_summary.txt.gz
(https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz,
~450MB compressed, updated weekly) -- not re-fetched automatically by
default (large, frequently-updated file); pass --url to download it, or
--input to point at an already-downloaded copy.

Output:
  - <output-prefix>_variants.tsv: filtered per-variant P/LP records
    (AlleleID, Type, Name, GeneSymbol, ClinicalSignificance,
    ClinSigSimple, Assembly, Chromosome, Start, Stop, ReviewStatus,
    VariationID)
  - <output-prefix>_gene_counts.tsv: per-gene summary
    (Gene, n_pathogenic_likely_pathogenic, n_plp_ge1star)
"""

import argparse
import gzip
import os
import sys
import urllib.request

import pandas as pd

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

KEEP_COLUMNS = [
    "AlleleID", "Type", "Name", "GeneSymbol", "ClinicalSignificance",
    "ClinSigSimple", "Assembly", "Chromosome", "Start", "Stop",
    "ReviewStatus", "VariationID",
]

ZERO_STAR_STATUSES = {
    "no assertion criteria provided",
    "no assertion provided",
    "no classification provided",
    "no classification for the individual variant",
}


def is_plp(clinical_significance):
    text = str(clinical_significance).lower()
    return "pathogenic" in text and "conflicting" not in text and "uncertain" not in text


def is_ge1star(review_status):
    return str(review_status).strip().lower() not in ZERO_STAR_STATUSES


def download(url, dest):
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, dest)
    print(f"✓ Downloaded to {dest} ({os.path.getsize(dest) / 1e6:.1f} MB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--input", help="Path to an already-downloaded variant_summary.txt.gz")
    src.add_argument("--url", nargs="?", const=CLINVAR_URL,
                      help=f"Download the release from this URL (default if flag given with no value: {CLINVAR_URL})")
    ap.add_argument("--download-to", default="variant_summary.txt.gz",
                     help="Where to save the file when --url is used (default: variant_summary.txt.gz)")
    ap.add_argument("--assembly", default="GRCh38", help="Assembly to keep (default: GRCh38)")
    ap.add_argument("--output-prefix", default="clinvar_plp",
                     help="Prefix for the two output TSVs (default: clinvar_plp)")
    args = ap.parse_args()

    if args.url:
        download(args.url, args.download_to)
        input_path = args.download_to
    else:
        input_path = args.input
        if not os.path.exists(input_path):
            print(f"Error: File not found: {input_path}", file=sys.stderr)
            sys.exit(1)

    print(f"Reading {input_path} ...")
    opener = gzip.open if input_path.endswith(".gz") else open
    with opener(input_path, "rt", errors="replace") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
    missing = [c for c in KEEP_COLUMNS if c not in header]
    if missing:
        print(f"Error: expected column(s) not found in {input_path}: {missing}", file=sys.stderr)
        sys.exit(1)
    usecols = [header.index(c) for c in KEEP_COLUMNS]

    df = pd.read_csv(input_path, sep="\t", header=0, usecols=usecols, dtype=str, low_memory=False)
    # pandas keeps the leading '#' on the first header cell; normalize
    df.columns = [c.lstrip("#") for c in df.columns]
    df = df[KEEP_COLUMNS]

    n_raw = len(df)
    print(f"Raw rows read: {n_raw}")

    df = df[df["Assembly"] == args.assembly]
    print(f"Rows after Assembly == {args.assembly}: {len(df)}")

    plp_mask = df["ClinicalSignificance"].apply(is_plp)
    df_plp = df[plp_mask].copy()
    print(f"Rows classified P/LP by the text rule: {len(df_plp)}")

    df_plp["ge1star"] = df_plp["ReviewStatus"].apply(is_ge1star)

    variants_out = f"{args.output_prefix}_variants.tsv"
    df_plp.drop(columns=["ge1star"]).to_csv(variants_out, sep="\t", index=False)
    print(f"\n✓ Saved {len(df_plp)} P/LP variant records to {variants_out}")

    # ClinVar's GeneSymbol field is not always HGNC-cased (e.g. "C9orf72"
    # rather than the HGNC-official "C9ORF72"), which silently loses
    # ~10 genes on a case-sensitive join against the master table. Upper-case
    # to match HGNC-style symbols; this is only a cosmetic normalization at
    # the ORF-suffix, not a semantic change.
    df_plp["GeneSymbol"] = df_plp["GeneSymbol"].str.upper()

    gene_counts = (
        df_plp.groupby("GeneSymbol")
        .agg(
            n_pathogenic_likely_pathogenic=("VariationID", "count"),
            n_plp_ge1star=("ge1star", "sum"),
        )
        .reset_index()
        .rename(columns={"GeneSymbol": "Gene"})
    )
    counts_out = f"{args.output_prefix}_gene_counts.tsv"
    gene_counts.to_csv(counts_out, sep="\t", index=False)
    print(f"✓ Saved per-gene P/LP counts for {len(gene_counts)} genes to {counts_out}")


if __name__ == "__main__":
    main()
