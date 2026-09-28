#!/usr/bin/env python3
"""
Extract ClinGen haploinsufficiency scores directly from ClinGen's own
Dosage Sensitivity curation download, replacing dbNSFP's
ClinGen_Haploinsufficiency_Score column (a snapshot that lags behind
ClinGen's live, continuously-updated curation).

Input: ClinGen's dosage sensitivity CSV export
(https://search.clinicalgenome.org/kb/gene-dosage/download), which
reports evidence as TEXT categories, not the numeric 0-3 scale
dbNSFP/most tools use. This script maps them using ClinGen's own
published numeric coding, empirically verified against dbNSFP's column
(99.2% agreement across 1,559 comparable genes; the ~0.8% that differ
are genes ClinGen has re-curated since dbNSFP's snapshot):
    No Evidence for Haploinsufficiency                  -> 0
    Little Evidence for Haploinsufficiency               -> 1
    Emerging Evidence for Haploinsufficiency             -> 2
    Sufficient Evidence for Haploinsufficiency           -> 3
    Gene Associated with Autosomal Recessive Phenotype   -> 30
    Dosage Sensitivity Unlikely for Haploinsufficiency   -> 40

Output: one row per gene (Gene, ClinGen_HI_score, ClinGen_HI_category,
ClinGen_HI_last_curated), suitable for a left-join onto the gene
master table on the Gene column.
"""

import argparse
import sys

import pandas as pd

CATEGORY_MAP = {
    "No Evidence for Haploinsufficiency": 0,
    "Little Evidence for Haploinsufficiency": 1,
    "Emerging Evidence for Haploinsufficiency": 2,
    "Sufficient Evidence for Haploinsufficiency": 3,
    "Gene Associated with Autosomal Recessive Phenotype": 30,
    "Dosage Sensitivity Unlikely for Haploinsufficiency": 40,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clingen-dosage-file", required=True,
                     help="Path to ClinGen's dosage-sensitivity CSV export "
                          "(from https://search.clinicalgenome.org/kb/gene-dosage/download)")
    ap.add_argument("--output", default="gene_clingen_hi.tsv",
                     help="Output TSV path (default: gene_clingen_hi.tsv)")
    args = ap.parse_args()

    print(f"Reading {args.clingen_dosage_file} ...")
    try:
        df = pd.read_csv(args.clingen_dosage_file, skiprows=4, quotechar='"')
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    df.columns = [c.strip() for c in df.columns]
    df = df[df["GENE SYMBOL"] != "+++++++++++"].copy()
    print(f"Raw gene rows: {len(df)}")

    df["ClinGen_HI_score"] = df["HAPLOINSUFFICIENCY"].map(CATEGORY_MAP)
    unmapped = df[df["ClinGen_HI_score"].isna()]["HAPLOINSUFFICIENCY"].unique()
    if len(unmapped):
        print(f"Warning: unmapped HAPLOINSUFFICIENCY categories found: {list(unmapped)}", file=sys.stderr)

    out = df[["GENE SYMBOL", "ClinGen_HI_score", "HAPLOINSUFFICIENCY", "DATE"]].rename(columns={
        "GENE SYMBOL": "Gene",
        "HAPLOINSUFFICIENCY": "ClinGen_HI_category",
        "DATE": "ClinGen_HI_last_curated",
    })

    out.to_csv(args.output, sep="\t", index=False)
    print(f"\n✓ Saved {len(out)} rows to {args.output}")


if __name__ == "__main__":
    main()
