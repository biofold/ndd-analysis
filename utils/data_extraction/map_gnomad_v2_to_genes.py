#!/usr/bin/env python3
"""
Attach gnomAD v2.1.1 constraint (pLI, LOEUF) to the 19,354-gene universe.

WHY THIS SCRIPT EXISTS
----------------------
The MOE validation analyses (R1.3, R2.10) used pLI from dbNSFP 5.2, which joined
gnomAD to genes by SYMBOL. Against gnomAD's own v2.1.1 file that left 1,385 scored
genes empty (almost all HGNC renames after gnomAD's 2018 naming, e.g. AARS1/AARS)
and gave 18 genes another gene's value (TUBB3, PI4K2A). dbNSFP also carries no
LOEUF. This script joins gnomAD's own table (data/gene_gnomad_v2_constraint.tsv,
from extract_gnomad_constraint.py) to the gene universe with the same rule the
iNDDx loader uses:

  1. Ensembl gene ID: our HGNC symbol -> HGNC complete set -> ensembl_gene_id ->
     gnomAD gene_id.
  2. Otherwise the gene symbol, but only when that symbol names exactly ONE gnomAD
     gene (98 v2.1.1 symbols map to more than one gene_id and are never used).

Each gnomAD gene is assigned to at most one of our genes. The match method is
recorded per gene so the join can be audited.

Output columns: Gene, Ensembl_ID, gnomad_v2_match, gnomAD_pLI, gnomAD_LOEUF,
gnomAD_LOEUF_decile, gnomAD_oe_lof, gnomAD_constraint_flag.

Usage:
  python3 utils/data_extraction/map_gnomad_v2_to_genes.py \\
      --genes data/gene_all_score.txt \\
      --hgnc data/raw/hgnc_complete_set_20260930.txt.gz \\
      --gnomad data/gene_gnomad_v2_constraint.tsv \\
      --output data/gene_gnomad_v2_mapped.tsv
"""

import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="gene universe (first column = HGNC symbol)")
    ap.add_argument("--hgnc", required=True, help="HGNC hgnc_complete_set.txt[.gz]")
    ap.add_argument("--gnomad", required=True, help="data/gene_gnomad_v2_constraint.tsv")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t", comment=None, dtype=str)
    genes = genes.rename(columns={genes.columns[0]: "Gene"})[["Gene"]]
    genes["key"] = genes["Gene"].str.upper()

    hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc[hgnc["status"] == "Approved"][["symbol", "ensembl_gene_id"]]
    hgnc["key"] = hgnc["symbol"].str.upper()
    hgnc = hgnc.drop_duplicates("key")

    gn = pd.read_csv(args.gnomad, sep="\t", dtype=str)
    gn["key"] = gn["Gene"].str.upper()
    value_cols = ["gnomad_v2_pli", "gnomad_v2_loeuf", "gnomad_v2_loeuf_decile",
                  "gnomad_v2_oe_lof", "gnomad_v2_constraint_flag"]

    # 1. Ensembl ID via HGNC
    g = genes.merge(hgnc[["key", "ensembl_gene_id"]], on="key", how="left")
    by_ens = gn.drop_duplicates("Ensembl_ID").set_index("Ensembl_ID")
    g["Ensembl_ID"] = g["ensembl_gene_id"]
    hit = g["Ensembl_ID"].isin(by_ens.index)
    g.loc[hit, "gnomad_gene_id"] = g.loc[hit, "Ensembl_ID"]
    g.loc[hit, "gnomad_v2_match"] = "ensembl"

    # 2. unambiguous symbol, for genes still unmatched, never re-using a gnomAD gene
    sym_counts = gn.groupby("key")["Ensembl_ID"].nunique()
    unambiguous = gn[gn["key"].map(sym_counts) == 1].drop_duplicates("key").set_index("key")
    claimed = set(g["gnomad_gene_id"].dropna())
    todo = g["gnomad_gene_id"].isna() & g["key"].isin(unambiguous.index)
    for i in g.index[todo]:
        gid = unambiguous.loc[g.at[i, "key"], "Ensembl_ID"]
        if gid not in claimed:
            g.at[i, "gnomad_gene_id"] = gid
            g.at[i, "gnomad_v2_match"] = "symbol"
            claimed.add(gid)

    vals = by_ens[value_cols]
    g = g.join(vals, on="gnomad_gene_id")
    dup = g["gnomad_gene_id"].dropna().duplicated().sum()
    if dup:
        sys.exit(f"Error: {dup} gnomAD genes assigned to more than one gene")

    out = g.rename(columns={
        "gnomad_v2_pli": "gnomAD_pLI", "gnomad_v2_loeuf": "gnomAD_LOEUF",
        "gnomad_v2_loeuf_decile": "gnomAD_LOEUF_decile", "gnomad_v2_oe_lof": "gnomAD_oe_lof",
        "gnomad_v2_constraint_flag": "gnomAD_constraint_flag",
    })[["Gene", "Ensembl_ID", "gnomad_v2_match", "gnomAD_pLI", "gnomAD_LOEUF",
        "gnomAD_LOEUF_decile", "gnomAD_oe_lof", "gnomAD_constraint_flag"]]
    out.to_csv(args.output, sep="\t", index=False)

    scored = out["gnomAD_pLI"].notna()
    print(f"genes: {len(out):,}")
    print(f"matched to a gnomAD v2.1.1 gene: {out['gnomad_v2_match'].notna().sum():,} "
          f"(ensembl {int((out['gnomad_v2_match'] == 'ensembl').sum()):,}, "
          f"symbol {int((out['gnomad_v2_match'] == 'symbol').sum()):,})")
    print(f"with pLI and LOEUF: {int(scored.sum()):,}")
    print(f"written: {args.output}")


if __name__ == "__main__":
    main()
