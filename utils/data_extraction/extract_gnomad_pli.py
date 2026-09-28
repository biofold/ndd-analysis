#!/usr/bin/env python3
"""
Extract gnomAD pLI (probability of loss-of-function intolerance) and a
handful of related haploinsufficiency/constraint scores from the dbNSFP
gene-level annotation table.

This is the actual source of the gnomAD_pLI, RVIS_percentile_ExAC,
HIPred, GHIS and ClinGen_Haploinsufficiency_Score columns that end up in
ndd_master_table.tsv -- previously that step was done interactively and
was not reproducible from the repo. This script makes it a first-class,
version-controlled pipeline step (used for R1.3's independent-validation
analysis, moe_pli_validation.py).

Input: dbNSFP5.2_gene.gz (gene-level summary file distributed with
dbNSFP, https://dbnsfp.s3.amazonaws.com/ or https://sites.google.com/site/jpopgen/dbNSFP).
It is NOT re-fetched by this script (dbNSFP requires an academic license
click-through and is a multi-GB download for the full annotation set;
only the small gene-level summary file, ~20-25MB gzipped, is needed
here) -- point --dbnsfp-gene-file at a local copy.

Output: a compact per-gene TSV (Gene, gnomAD_pLI, RVIS_percentile_ExAC,
HIPred, GHIS, ClinGen_Haploinsufficiency_Score) suitable for a left-join
onto the gene master table on the Gene_name column.
"""

import argparse
import gzip
import sys

import pandas as pd

COLUMNS_TO_KEEP = [
    "Gene_name",
    "gnomAD_pLI",
    "RVIS_percentile_ExAC",
    "HIPred",
    "GHIS",
    "ClinGen_Haploinsufficiency_Score",
]


def load_dbnsfp_gene_table(path):
    """Stream-read only the needed columns from the (possibly large) dbNSFP
    gene-level gzip file, avoiding loading every one of its ~700 columns."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    missing = [c for c in COLUMNS_TO_KEEP if c not in header]
    if missing:
        print(f"Error: expected column(s) not found in {path}: {missing}", file=sys.stderr)
        sys.exit(1)
    usecols = [header.index(c) for c in COLUMNS_TO_KEEP]
    df = pd.read_csv(path, sep="\t", usecols=usecols, na_values=".", dtype=str)
    df = df[COLUMNS_TO_KEEP]
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dbnsfp-gene-file", required=True,
                     help="Path to dbNSFP5.x_gene.gz (gene-level summary file)")
    ap.add_argument("--output", default="gene_gnomad_pli.tsv",
                     help="Output TSV path (default: gene_gnomad_pli.tsv)")
    args = ap.parse_args()

    print(f"Reading {args.dbnsfp_gene_file} ...")
    df = load_dbnsfp_gene_table(args.dbnsfp_gene_file)
    df = df.rename(columns={"Gene_name": "Gene"})

    for col in ["gnomAD_pLI", "RVIS_percentile_ExAC", "HIPred", "GHIS", "ClinGen_Haploinsufficiency_Score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    n_total = len(df)
    n_with_pli = df["gnomAD_pLI"].notna().sum()
    print(f"Genes in dbNSFP gene table: {n_total}")
    print(f"Genes with a non-missing gnomAD_pLI value: {n_with_pli} ({n_with_pli / n_total:.1%})")

    df.to_csv(args.output, sep="\t", index=False)
    print(f"\n✓ Saved {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()
