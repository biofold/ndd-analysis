#!/usr/bin/env python3
"""
extract_mondo.py

Full pipeline for gene-disease mapping:
  Step 1 (add_disease_code)       -> adds disease codes from ClinVar gene_condition_source_id
  Step 2 (retrieve_missing_mondo) -> recovers missing MONDO codes from MONDO_CUI_history
  Step 3 (retrieve_parent_mondo)  -> maps each MONDO code to its top-level parent category,
                                      using the MONDO ontology graph (mondo.json), and cleans
                                      up rows with deprecated/obsolete MONDO codes
  Step 4 (generate_gmt)           -> builds two .gmt files: one per specific MONDO disease
                                      term and one per top-level MONDO parent category
 
Usage:
    python utils/data_extraction/extract_mondo.py \
        --genes list_NDD_score.tsv \
        --clinvar data/raw/gene_condition_source_id_2026-07-14.tsv \
        --mondo-history data/raw/MONDO_CUI_history_2026-07-16.txt \
        --mondo-json data/raw/mondo_2026-07-06.json.gz \
        --gmt-prefix MONDO_2026 \
        --outdir results/
 
Input files (data/raw/):
  gene_condition_source_id_2026-07-14.tsv  ClinVar gene_condition_source_id, accessed 14 July 2026
  MONDO_CUI_history_2026-07-16.txt         MedGen MONDO_CUI_history, accessed 16 July 2026
  mondo_2026-07-06.json.gz                 MONDO ontology, release 2026-07-06 (gzip or plain JSON)

Requires add_disease_code.py, retrieve_missing_mondo.py,
retrieve_parent_mondo.py and generate_gmt.py in the same directory.
"""
import argparse
import gzip
import json
from pathlib import Path
import pandas as pd
 
from add_disease_code import add_disease_codes
from retrieve_missing_mondo import MONDO_CUI_dict, extract_missing_CUI
from retrieve_parent_mondo import add_mondo_parents, clean_deprecated_mondo
from generate_gmt import build_disease_gmt_lines, build_group_gmt_lines, save_gmt

def count_genes_with_mondo(df: pd.DataFrame) -> int:
    """Count unique genes that have at least one MONDO-linked annotation."""
    return df.loc[df['Source_Disease'] == 'MONDO', 'Gene_Symbol'].nunique()
 
 
def count_mondo_rows(df: pd.DataFrame) -> int:
    """Count total rows (disease annotations) linked to MONDO."""
    return len(df.loc[df['Source_Disease'] == 'MONDO'])
 
def run_pipeline(genes_file, clinvar_file, mondo_history_file, mondo_json_file,
                  outdir="results", keep_intermediate=True, gmt_prefix="MONDO_2026"):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
 
    # ---- Step 1: ClinVar -> MONDO/OMIM/... ----
    print("=== STEP 1: adding disease codes from ClinVar ===")
    genes_df = pd.read_csv(genes_file, sep='\t')
    clinvar_df = pd.read_csv(clinvar_file, sep='\t')
    step1_df = add_disease_codes(genes_df, clinvar_df)
 
    n_genes_step1 = count_genes_with_mondo(step1_df)
    n_annotations_step1 = count_mondo_rows(step1_df)
    print(f"Genes with at least one MONDO annotation after step 1: {n_genes_step1}, "
          f"with a total of {n_annotations_step1} MONDO annotations")
 
    step1_path = outdir / 'list_NDD_diseases.tsv'
    if keep_intermediate:
        step1_df.to_csv(step1_path, sep='\t', index=False)
        print(f"Intermediate file saved: {step1_path}")
 
    # ---- Step 2: recover missing MONDO codes ----
    print("=== STEP 2: recovering missing MONDO codes ===")
    mondo_dict = MONDO_CUI_dict(mondo_history_file)
    print(f"Loaded {len(mondo_dict)} current CUI->MONDO mappings (IsCurrent==1)")
    final_df = extract_missing_CUI(step1_df, mondo_dict)
 
    n_genes_step2 = count_genes_with_mondo(final_df)
    n_annotations_step2 = count_mondo_rows(final_df)
    print(f"Genes with at least one MONDO annotation after step 2: {n_genes_step2} "
          f"(+{n_genes_step2 - n_genes_step1} genes), with a total of {n_annotations_step2} "
          f"MONDO annotations (+{n_annotations_step2 - n_annotations_step1})")
 
    final_path = outdir / 'updated_list_NDD_diseases.tsv'
    final_df.to_csv(final_path, sep='\t', index=False)
    print(f"Step 2 output saved: {final_path}")
 
    # ---- Step 3: map MONDO codes to their top-level parent category ----
    print("=== STEP 3: mapping MONDO codes to parent categories ===")
    opener = gzip.open if str(mondo_json_file).endswith('.gz') else open
    with opener(mondo_json_file, 'rt') as f:
        mondo_data = json.load(f)
    with_parents_df = add_mondo_parents(final_df, mondo_data)
 
    with_parents_path = outdir / 'updated_list_NDD_diseases_with_parents.tsv'
    if keep_intermediate:
        with_parents_df.to_csv(with_parents_path, sep='\t', index=False)
        print(f"Intermediate file saved: {with_parents_path}")
 
    cleaned_df = clean_deprecated_mondo(with_parents_df)
 
    final_path = outdir / 'updated_list_NDD_diseases_with_parents_cleaned.tsv'
    cleaned_df.to_csv(final_path, sep='\t', index=False)
    print(f"Step 3 output saved: {final_path}")
 
    # ---- Step 4: build .gmt files (per-disease and per-parent-group) ----
    print("=== STEP 4: building .gmt files ===")
    disease_lines = build_disease_gmt_lines(cleaned_df)
    disease_gmt_path = outdir / f"{gmt_prefix}.gmt"
    save_gmt(disease_lines, disease_gmt_path)
    print(f"{len(disease_lines)} disease terms written to {disease_gmt_path}")
 
    group_lines = build_group_gmt_lines(cleaned_df)
    group_gmt_path = outdir / f"{gmt_prefix}_GROUPS.gmt"
    save_gmt(group_lines, group_gmt_path)
    print(f"{len(group_lines)} parent groups written to {group_gmt_path}")
 
    print(f"Pipeline completed. Final outputs: {disease_gmt_path}, {group_gmt_path}")
 
    return cleaned_df
 
 
def parse_args():
    p = argparse.ArgumentParser(description="Gene-disease mapping pipeline (ClinVar + MONDO)")
    p.add_argument('--genes', required=True, help="TSV with the NDD gene list and score (Gene_Name column)")
    p.add_argument('--clinvar', required=True, help="ClinVar gene_condition_source_id TSV")
    p.add_argument('--mondo-history', required=True, help="MONDO_CUI_history.txt file")
    p.add_argument('--mondo-json', required=True, help="MONDO ontology JSON file (mondo.json or mondo.json.gz)")
    p.add_argument('--gmt-prefix', default='MONDO_2026',
                    help="Prefix for the output .gmt files (default: MONDO_2026 -> MONDO_2026.gmt / MONDO_2026_GROUPS.gmt)")
    p.add_argument('--outdir', default='results', help="Output folder (default: ./results)")
    p.add_argument('--no-keep-intermediate', action='store_true',
                    help="Do not save the step 1 intermediate file to disk")
    return p.parse_args()
 
 
def main():
    args = parse_args()
    run_pipeline(
        genes_file=args.genes,
        clinvar_file=args.clinvar,
        mondo_history_file=args.mondo_history,
        mondo_json_file=args.mondo_json,
        outdir=args.outdir,
        keep_intermediate=not args.no_keep_intermediate,
        gmt_prefix=args.gmt_prefix,
    )
 
 
if __name__ == '__main__':
    main()