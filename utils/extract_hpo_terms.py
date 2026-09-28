#!/usr/bin/env python3
"""
Extract per-gene HPO (Human Phenotype Ontology) term lists directly from
the HPO project's own gene-to-phenotype annotation file.

This replaces the prior source of ndd_master_table.tsv's HPO_id/HPO_name
columns, which came from dbNSFP 5.2's flattened re-export of HPO data
rather than from HPO itself -- dbNSFP repackages the file on its own
snapshot schedule with its own flattening logic, and was found not to
exactly match a fresh pull from HPO (e.g. SCN2A: 208 unique HPO terms
extracted here vs 162 from the dbNSFP-derived column). Going to the
primary source removes that dependency, mirroring what
extract_gnomad_pli.py / extract_clinvar_plp.py already do for pLI and
ClinVar.

Input: HPO's genes_to_phenotype.txt
(https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt
or https://hpo.jax.org/data/annotations), columns ncbi_gene_id,
gene_symbol, hpo_id, hpo_name, frequency, disease_id -- one row per
(gene, HPO term, source disease) triple, so the SAME (gene, hpo_id)
pair legitimately repeats once per disease that gene causes with that
phenotype (18.0% of all rows are such duplicates; deduplicating by
(gene_symbol, hpo_id) is mandatory, not optional).

HPO's gene_symbol field is not always HGNC-cased (e.g. "C9orf72" rather
than the HGNC-official "C9ORF72"), the same issue found and fixed for
ClinVar's GeneSymbol field -- upper-cased here for the same reason.
A further ~67 genes present in the HPO file cannot be matched even
after that fix; inspection shows these are mitochondrial tRNA genes,
microRNAs, snoRNAs and lncRNAs, i.e. legitimately outside the
19,354-gene HGNC protein-coding universe used throughout this
pipeline, not a bug.

Output: one row per gene (Gene, HPO_id, HPO_name, n_hpo_terms), with
HPO_id/HPO_name semicolon-joined in first-occurrence order, suitable
for a left-join onto the gene master table on the Gene column.
"""

import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-to-phenotype-file", required=True,
                     help="Path to HPO's genes_to_phenotype.txt")
    ap.add_argument("--output", default="gene_hpo_terms.tsv",
                     help="Output TSV path (default: gene_hpo_terms.tsv)")
    args = ap.parse_args()

    print(f"Reading {args.genes_to_phenotype_file} ...")
    try:
        df = pd.read_csv(args.genes_to_phenotype_file, sep="\t")
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    n_raw = len(df)
    print(f"Raw rows: {n_raw}")

    df["gene_symbol"] = df["gene_symbol"].str.upper()

    n_dup = df.duplicated(subset=["gene_symbol", "hpo_id"]).sum()
    print(f"Duplicate (gene_symbol, hpo_id) rows (same term via >1 source disease): "
          f"{n_dup} ({n_dup / n_raw:.1%})")

    df = df.drop_duplicates(subset=["gene_symbol", "hpo_id"], keep="first")

    def agg_gene(g):
        return pd.Series({
            "HPO_id": ";".join(g["hpo_id"]),
            "HPO_name": ";".join(g["hpo_name"]),
            "n_hpo_terms": len(g),
        })

    gene_hpo = df.groupby("gene_symbol", sort=True).apply(agg_gene, include_groups=False)
    gene_hpo = gene_hpo.reset_index().rename(columns={"gene_symbol": "Gene"})

    print(f"\nGenes with >=1 HPO term: {len(gene_hpo)}")
    print(f"Median terms/gene: {gene_hpo['n_hpo_terms'].median():.0f}")

    gene_hpo.to_csv(args.output, sep="\t", index=False)
    print(f"\n✓ Saved {len(gene_hpo)} rows to {args.output}")


if __name__ == "__main__":
    main()
