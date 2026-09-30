#!/usr/bin/env python
"""
The project's single gnomAD constraint extractor: per-gene pLI, LOEUF, o/e LoF and
LOEUF decile from gnomAD's OWN v2.1.1 gene constraint release.

It replaces extract_gnomad_v4_constraint.py (gnomAD v4.x is not used anywhere in
iNDDx) and supersedes extract_gnomad_pli.py as the source of pLI: that script takes
pLI from the dbNSFP 5.2 gene table, which matched gnomAD by SYMBOL. Against gnomAD's
own file (joined on Ensembl ID) the values agree for 99.89% of the 17,056 genes both
score, but dbNSFP
  * holds ANOTHER gene's value for 18 genes: symbols with two gnomAD rows, where
    dbNSFP took the other row (e.g. TUBB3 0.000 vs 0.968, PI4K2A 0.000 vs 0.804);
  * is EMPTY for 1,385 genes gnomAD scores, 1,382 of them because the HGNC symbol
    changed after gnomAD's 2018 naming (AARS1/AARS, ABRAXAS1/FAM175A, ACKR1/DARC).
dbNSFP also carries no LOEUF.

CONSUMERS
---------
The output feeds the iNDDx API `gnomad_v2` block (api/db/load_mongo.py) and the
website's pLI/LOEUF (tools/make_ndd_data_js.py), and is the input for recalculating
the pLI analyses on gnomAD v2.1.1. Every pipeline uses gnomAD v2.1.1 only; the
dbNSFP-derived data/gene_gnomad_pli.tsv is kept as a record and is not read.

INPUT
-----
data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz -- one row per gene (canonical
transcript), 19,704 rows, GRCh37 / GENCODE v19. It covers chrX and chrY (723 and
23 genes with a value).

COLUMNS EMITTED
---------------
Ensembl_ID (gnomAD gene_id), Gene (symbol, upper-cased), transcript,
gnomad_v2_pli, gnomad_v2_loeuf (oe_lof_upper), gnomad_v2_oe_lof,
gnomad_v2_loeuf_decile (oe_lof_upper_bin, 0 = most constrained),
gnomad_v2_constraint_flag (gnomAD's own reliability flag; emitted, not used to
drop rows).

Joining to iNDDx is done by the consumer (api/db/load_mongo.py), on Ensembl ID
first and an UNAMBIGUOUS symbol second, because 98 v2.1.1 symbols (88 among
genes with a LOEUF) map to more than one gene_id.

Usage
-----
python utils/data_extraction/extract_gnomad_constraint.py \\
    --input data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \\
    --output data/gene_gnomad_v2_constraint.tsv
"""
import argparse
import sys

import pandas as pd

MAP = {
    "gene_id": "Ensembl_ID",
    "gene": "Gene",
    "transcript": "transcript",
    "pLI": "gnomad_v2_pli",
    "oe_lof_upper": "gnomad_v2_loeuf",
    "oe_lof": "gnomad_v2_oe_lof",
    "oe_lof_upper_bin": "gnomad_v2_loeuf_decile",
    "constraint_flag": "gnomad_v2_constraint_flag",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", compression="gzip", low_memory=False)
    missing = [c for c in MAP if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: expected columns absent from input: {missing}")
    print(f"{args.input}: {len(df):,} rows")

    out = df[list(MAP)].rename(columns=MAP)
    if not out["Ensembl_ID"].is_unique:
        sys.exit("ERROR: gene_id is not unique -- the file is expected to hold one row per gene")
    out["Gene"] = out["Gene"].astype(str).str.upper()
    out["gnomad_v2_loeuf_decile"] = out["gnomad_v2_loeuf_decile"].astype("Int64")
    # pLI and LOEUF are computed together; one present without the other means a
    # format change upstream, not a real gene state.
    both = out["gnomad_v2_pli"].notna() == out["gnomad_v2_loeuf"].notna()
    if not both.all():
        sys.exit(f"ERROR: {int((~both).sum())} genes have pLI xor LOEUF")
    out = out.sort_values("Ensembl_ID")

    print(f"  genes: {len(out):,} | with LOEUF: {out['gnomad_v2_loeuf'].notna().sum():,}")
    print(f"  symbols mapping to >1 gene_id: {int(out['Gene'].duplicated(keep=False).sum())}")
    print(f"  with a constraint_flag: {int(out['gnomad_v2_constraint_flag'].notna().sum()):,} (emitted, not dropped)")
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
