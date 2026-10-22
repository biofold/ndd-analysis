#!/usr/bin/env python3
"""
Build the per-gene TADA-NDD outcome for the candidate set (R1.3).

Source: Fu et al. 2022, Nat Genet 54:1320-1331 (doi:10.1038/s41588-022-01104-0),
Supplementary Table 11 ("Gene association results"), stored verbatim as
data/raw/fu2022_supplementary_table11.tsv. The table gives, for 18,128
autosomal genes, the false discovery rate of the extended TADA model for
ASD, DD and the combined NDD analysis (FDR_TADA_NDD).

Gene matching, in order (no value is assigned to more than one gene):
  1. Ensembl gene ID (gene_id, version stripped) -> our gene symbol, using the
     Ensembl IDs of data/gene_gnomad_v2_mapped.tsv;
  2. otherwise the GENCODE v33 symbol, then the Fu symbol, resolved through
     data/gene_id_resolution.tsv.gz (status == "resolved" only).
If several Fu rows map to the same gene, the row with the smallest
p_TADA_NDD is kept.

Output (candidate genes with a TADA-NDD result):
  Gene, MOE_score, FDR_TADA_NDD, TADA_NDD (1 if FDR <= --fdr, default 0.1),
  TADA_NDD_fdr005 (1 if FDR <= 0.05).
"""

import argparse
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv (Gene, Class, MOE_score)")
    ap.add_argument("--fu-table", default="data/raw/fu2022_supplementary_table11.tsv")
    ap.add_argument("--ensembl-map", default="data/gene_gnomad_v2_mapped.tsv",
                    help="Table with Gene and Ensembl_ID columns")
    ap.add_argument("--id-resolution", default="data/gene_id_resolution.tsv.gz")
    ap.add_argument("--fdr", type=float, default=0.1, help="FDR threshold defining NDD-associated genes")
    ap.add_argument("--class-filter", default="candidate")
    ap.add_argument("--output", default="moe_tada_input.tsv")
    args = ap.parse_args()

    for path in (args.master_table, args.fu_table, args.ensembl_map, args.id_resolution):
        if not os.path.exists(path):
            print(f"Error: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    fu = pd.read_csv(args.fu_table, sep="\t", float_precision="round_trip")
    ens = pd.read_csv(args.ensembl_map, sep="\t", usecols=["Gene", "Ensembl_ID"])
    res = pd.read_csv(args.id_resolution, sep="\t", dtype=str)

    e2g = ens.dropna().drop_duplicates("Ensembl_ID").set_index("Ensembl_ID")["Gene"]
    fu["Gene"] = fu["gene_id"].astype(str).str.split(".").str[0].map(e2g)
    resolved = (res[(res["status"] == "resolved") & res["gene"].notna()]
                .drop_duplicates("query_key").set_index("query_key")["gene"])
    for col in ("gene_gencodeV33", "gene"):
        missing = fu["Gene"].isna()
        fu.loc[missing, "Gene"] = fu.loc[missing, col].map(resolved)
    n_mapped = int(fu["Gene"].notna().sum())
    fu = fu.dropna(subset=["Gene"]).sort_values("p_TADA_NDD").drop_duplicates("Gene")

    master = pd.read_csv(args.master_table, sep="\t", dtype=str)
    master["MOE_score"] = pd.to_numeric(master["MOE_score"], errors="coerce")
    if args.class_filter:
        master = master[master["Class"] == args.class_filter]
    out = master[["Gene", "MOE_score"]].dropna().merge(
        fu[["Gene", "FDR_TADA_NDD"]], on="Gene", how="inner")
    out["TADA_NDD"] = (out["FDR_TADA_NDD"] <= args.fdr).astype(int)
    out["TADA_NDD_fdr005"] = (out["FDR_TADA_NDD"] <= 0.05).astype(int)
    out = out.sort_values("Gene")
    out.to_csv(args.output, sep="\t", index=False, float_format="%.6g")

    print(f"Fu et al. rows: {len(pd.read_csv(args.fu_table, sep=chr(9), usecols=['gene'])):,}; "
          f"mapped to our genes: {n_mapped:,}")
    print(f"{args.class_filter or 'all'} genes with a TADA-NDD result: {len(out):,}; "
          f"FDR <= {args.fdr}: {int(out['TADA_NDD'].sum()):,}; "
          f"FDR <= 0.05: {int(out['TADA_NDD_fdr005'].sum()):,}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
