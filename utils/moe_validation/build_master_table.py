#!/usr/bin/env python3
"""
Assemble the per-gene master table (ndd_master_table.tsv) that the
reviewer-requested validation analyses run on, and emit the two small
derived input tables those analyses consume.

WHY THIS SCRIPT EXISTS
----------------------
The validation utilities written in response to reviewer points R1.3
(independent validation of the MOE >= 4 threshold), R2.9 (is MOE just
tracking annotation density?) and R2.10 (redundancy among MOE
components) -- moe_clinvar_validation.py, moe_pli_validation.py,
moe_annotation_bias_baseline.py, plus the generic violin/ks-test/
bar/fisher-test/roc-pr plotters -- all read a single merged per-gene
table. Until now that table existed only as a hand-assembled file
outside the repository, which meant the figures in the response to the
reviewers could not be regenerated from `run_ndd_analysis.py` alone.
Reviewers of a revision are entitled to expect that every number and
figure in the response is reproducible from the deposited pipeline, so
the merge is done here, from files that are all under version control.

WHAT IS MERGED
--------------
Gene classification and source membership come from the combined score
file (data/gene_all_score.txt), whose `Annotation` field is a
pipe-separated list of the resources that nominated each gene
(hc, candidate, neuro, develop, sanchis, sfari, sfari+). These are
expanded into one boolean column per source so that the per-source
contribution asked for in R2.3 can be recomputed from the table.

The MOE score and its five per-library indicator components come from
the pipeline's own supercandidate.tsv (scripts/2_supercandidate.py).
Note that the MOE score is defined only for the candidate class: it
scores candidate genes against the term sets that are significantly
enriched in the high-confidence (curated) class, so genes outside the
candidate set get an empty MOE_score and are dropped by the validation
scripts. This is intentional -- scoring the curated genes against
enrichment derived from the curated genes themselves would be circular.

The independent validation labels are joined in from data/:
  - gnomAD_pLI (gene_gnomad_pli.tsv) -- population constraint, derived
    from allele frequencies alone, so independent of literature
    curation and of the GO/pathway annotation MOE is built from.
  - ClinVar P/LP counts (clinvar_plp_gene_counts.tsv) -- clinical
    variant curation, independent of pathway enrichment.
  - GO annotation counts and PubMed counts (gene_go_annotation_counts.tsv,
    gene_pubmed_counts.tsv) -- the two "how well studied is this gene"
    proxies used as competing baselines for R2.9.
  - SysNDD score and MONDO terms -- the disease-level annotations
    already used in the manuscript.

DERIVED INPUT TABLES
--------------------
Two narrow tables are written alongside the master table, restricted to
genes with a defined MOE score, in the 1-based column layout the generic
plotters expect (category/predictor in column 2, value in column 3):
  moe_pli_input.tsv      Gene, MOE_score, gnomAD_pLI      (violin, ks-test)
  moe_clinvar_input.tsv  Gene, MOE_score, ClinVar_PLP     (bar, fisher-test, roc-pr)

Usage:
  python3 utils/build_master_table.py \
      --data-dir data \
      --supercandidate results/main/supercandidate.tsv \
      --output-dir results/main
"""

import argparse
import os
import sys

import pandas as pd

# Annotation token -> master-table column name. The tokens are the ones used
# in data/gene_all_score.txt and in the subgroup lists passed to venn-plot.py
# and upset-plot.py, so the two stay consistent.
SOURCE_COLUMNS = [
    ("hc", "GeneTrek_HC"),
    ("candidate", "GeneTrek_LC"),
    ("neuro", "Orphanet_neuro"),
    ("develop", "Orphanet_develop"),
    ("sanchis", "SanchisJuan"),
    ("sfari+", "SFARI_plus"),
    ("sfari", "SFARI_other"),
]

# MOE library -> component column. Mirrors config.yml's libraries_supercandidate.
MOE_COMPONENTS = [
    ("GO_Biological_Process_2026", "MOE_GOBP"),
    ("GO_Cellular_Component_2026", "MOE_GOCC"),
    ("GO_Molecular_Function_2026", "MOE_GOMF"),
    ("KEGG_2021_Human", "MOE_KEGG"),
    ("Reactome_Pathways_2024", "MOE_Reactome"),
]


def read_gene_list(path):
    """Read a one-symbol-per-line gene list into a set."""
    if not os.path.exists(path):
        sys.stderr.write(f"Warning: gene list not found, skipping: {path}\n")
        return set()
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def read_table(path, **kwargs):
    """Read a TSV if it exists, otherwise warn and return None."""
    if not os.path.exists(path):
        sys.stderr.write(f"Warning: table not found, column(s) will be empty: {path}\n")
        return None
    return pd.read_csv(path, sep="\t", dtype=str, **kwargs)


def merge_on_symbol(master, right, columns, gene_col="Gene"):
    """Left-join `columns` from `right` onto `master`, matching gene symbols
    case-insensitively.

    Gene symbols are not consistently cased across the source databases: HGNC
    (and therefore data/gene_all_score.txt, which is the spine of this table)
    writes the open-reading-frame genes as C10orf105, while ClinVar, the GO
    annotation counts, the HPO extract and every .gmt library in libs/ upper-case
    them to C10ORF105. A plain merge on the symbol silently drops those genes'
    values -- they come back as missing rather than as an error, which is the
    worst failure mode for a validation table. Matching on the upper-cased
    symbol fixes it while keeping the HGNC spelling in the output.

    scripts/2_supercandidate.py already upper-cases symbols when it scores
    genes against the libraries, so this makes the master table consistent with
    how the MOE score itself is computed.
    """
    right = right.copy()
    right["_key"] = right[gene_col].astype(str).str.strip().str.upper()
    right = right.drop_duplicates("_key")[["_key"] + columns]
    master = master.copy()
    master["_key"] = master["Gene"].astype(str).str.upper()
    merged = master.merge(right, on="_key", how="left").drop(columns="_key")
    return merged


def build(data_dir, supercandidate_file, sfari_file, output_dir):
    combined_file = os.path.join(data_dir, "gene_all_score.txt")
    if not os.path.exists(combined_file):
        sys.exit(f"Error: required input not found: {combined_file}")

    combined = pd.read_csv(combined_file, sep="\t", dtype=str)
    combined = combined.rename(columns={combined.columns[0]: "Gene"})
    combined["Gene"] = combined["Gene"].str.lstrip("#").str.strip()
    master = combined[["Gene", "Class"]].copy()

    # --- source membership, expanded from the pipe-separated Annotation field ---
    annotation = combined["Annotation"].fillna("none")
    tokens = annotation.apply(lambda s: set(s.split("|")) if s != "none" else set())
    for token, column in SOURCE_COLUMNS:
        master[column] = tokens.apply(lambda t, tok=token: tok in t)
    master["n_sources"] = tokens.apply(len)

    # --- MOE score and per-library components ---
    if not os.path.exists(supercandidate_file):
        sys.exit(f"Error: required input not found: {supercandidate_file}")
    sc = pd.read_csv(supercandidate_file, sep="\t", dtype=str)
    sc = sc.rename(columns={"Score": "MOE_score", "Terms": "MOE_libraries"})
    sc["MOE_score"] = pd.to_numeric(sc["MOE_score"], errors="coerce")
    sc["_key"] = sc["Gene"].astype(str).str.upper()
    sc_terms = sc.set_index("_key")["MOE_libraries"].fillna("")
    master = merge_on_symbol(master, sc, ["MOE_score"])
    key = master["Gene"].astype(str).str.upper()
    for library, column in MOE_COMPONENTS:
        hit = sc_terms.apply(lambda s, lib=library: int(lib in s.split("|")))
        master[column] = key.map(hit)

    # --- SFARI gene score (categorical evidence tier, 1 = highest) ---
    if os.path.exists(sfari_file):
        sfari = pd.read_csv(sfari_file, dtype=str)
        sfari = sfari.rename(columns={"gene-symbol": "Gene", "gene-score": "SFARI_score"})
        master = merge_on_symbol(master, sfari, ["SFARI_score"])
    else:
        sys.stderr.write(f"Warning: SFARI file not found: {sfari_file}\n")
        master["SFARI_score"] = pd.NA

    # --- cancer gene membership (KEGG pathways in cancer) ---
    cancer = read_gene_list(os.path.join(data_dir, "cancer.txt"))
    master["Cancer_KEGG"] = master["Gene"].isin(cancer)

    # --- SysNDD score ---
    sysndd = read_table(os.path.join(data_dir, "SysNDD_all.txt"))
    if sysndd is not None:
        sysndd = sysndd.rename(columns={"Score": "SysNDD_score"})
        master = merge_on_symbol(master, sysndd, ["SysNDD_score"])
    else:
        master["SysNDD_score"] = pd.NA

    # --- MONDO disease terms (headerless: gene <tab> pipe-separated terms) ---
    mondo_path = os.path.join(data_dir, "mondo.txt")
    if os.path.exists(mondo_path):
        mondo = pd.read_csv(mondo_path, sep="\t", dtype=str, header=None,
                            names=["Gene", "MONDO_terms"])
        master = merge_on_symbol(master, mondo, ["MONDO_terms"])
    else:
        sys.stderr.write(f"Warning: MONDO file not found: {mondo_path}\n")
        master["MONDO_terms"] = pd.NA

    # --- independent validation labels ---
    pli = read_table(os.path.join(data_dir, "gene_gnomad_pli.tsv"))
    if pli is not None:
        master = merge_on_symbol(master, pli, ["gnomAD_pLI"])
    else:
        master["gnomAD_pLI"] = pd.NA

    clinvar = read_table(os.path.join(data_dir, "clinvar_plp_gene_counts.tsv"))
    if clinvar is not None:
        clinvar = clinvar[clinvar["Gene"] != "-"]
        master = merge_on_symbol(
            master, clinvar,
            ["n_pathogenic_likely_pathogenic", "n_plp_ge1star"])
        for col in ("n_pathogenic_likely_pathogenic", "n_plp_ge1star"):
            master[col] = pd.to_numeric(master[col], errors="coerce").fillna(0).astype(int)
    else:
        master["n_pathogenic_likely_pathogenic"] = 0
        master["n_plp_ge1star"] = 0

    go_counts = read_table(os.path.join(data_dir, "gene_go_annotation_counts.tsv"))
    if go_counts is not None:
        master = merge_on_symbol(master, go_counts, ["n_go_annotations"])
        master["n_go_annotations"] = (
            pd.to_numeric(master["n_go_annotations"], errors="coerce").fillna(0).astype(int)
        )
    else:
        master["n_go_annotations"] = 0

    pubmed = read_table(os.path.join(data_dir, "gene_pubmed_counts.tsv"))
    if pubmed is not None:
        pubmed = pubmed.rename(columns={"Symbol": "Gene"})
        pubmed["n_pubmed"] = pd.to_numeric(pubmed["n_pubmed"], errors="coerce")
        # gene2pubmed keys on GeneID, so a few symbols occur twice; collapse to
        # one row per symbol keeping the largest count, matching the policy in
        # utils/moe_annotation_bias_baseline.py so the two agree on n.
        pubmed = pubmed.groupby("Gene", as_index=False)["n_pubmed"].max()
        master = merge_on_symbol(master, pubmed, ["n_pubmed"])
        master["n_pubmed"] = (
            pd.to_numeric(master["n_pubmed"], errors="coerce").fillna(0).astype(int)
        )
    else:
        master["n_pubmed"] = 0

    # --- HPO phenotype breadth and ClinGen haploinsufficiency ---
    # Only the counts/score are carried here; the full HPO term lists are large
    # and live in data/gene_hpo_terms.tsv for anyone who needs them.
    hpo = read_table(os.path.join(data_dir, "gene_hpo_terms.tsv"))
    if hpo is not None and "n_hpo_terms" in hpo.columns:
        hpo = hpo[hpo["Gene"] != "-"]
        master = merge_on_symbol(master, hpo, ["n_hpo_terms"])
        master["n_hpo_terms"] = (
            pd.to_numeric(master["n_hpo_terms"], errors="coerce").fillna(0).astype(int)
        )
    else:
        master["n_hpo_terms"] = 0

    clingen = read_table(os.path.join(data_dir, "gene_clingen_hi.tsv"))
    if clingen is not None and "ClinGen_HI_score" in clingen.columns:
        master = merge_on_symbol(master, clingen, ["ClinGen_HI_score"])
    else:
        master["ClinGen_HI_score"] = pd.NA

    # --- write the master table ---
    os.makedirs(output_dir, exist_ok=True)
    master_path = os.path.join(output_dir, "ndd_master_table.tsv")
    master.to_csv(master_path, sep="\t", index=False)
    print(f"Master table written: {master_path}  ({len(master):,} genes, "
          f"{len(master.columns)} columns)")

    # --- derived validation inputs, restricted to genes with a defined MOE score ---
    scored = master[master["MOE_score"].notna()].copy()
    scored["MOE_score"] = scored["MOE_score"].astype(int)

    pli_tab = scored[["Gene", "MOE_score", "gnomAD_pLI"]].copy()
    pli_tab["gnomAD_pLI"] = pd.to_numeric(pli_tab["gnomAD_pLI"], errors="coerce")
    pli_tab = pli_tab.dropna(subset=["gnomAD_pLI"])
    pli_path = os.path.join(output_dir, "moe_pli_input.tsv")
    pli_tab.to_csv(pli_path, sep="\t", index=False)
    print(f"pLI validation input written: {pli_path}  ({len(pli_tab):,} genes)")

    clinvar_tab = scored[["Gene", "MOE_score"]].copy()
    clinvar_tab["ClinVar_PLP"] = (
        scored["n_pathogenic_likely_pathogenic"].astype(int) > 0
    ).astype(int)
    clinvar_path = os.path.join(output_dir, "moe_clinvar_input.tsv")
    clinvar_tab.to_csv(clinvar_path, sep="\t", index=False)
    print(f"ClinVar validation input written: {clinvar_path}  "
          f"({len(clinvar_tab):,} genes, {int(clinvar_tab['ClinVar_PLP'].sum()):,} P/LP positive)")

    return master_path, pli_path, clinvar_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--data-dir", default="data",
                    help="Directory holding the source data files (default: data)")
    ap.add_argument("--supercandidate", default="results/main/supercandidate.tsv",
                    help="MOE score file from scripts/2_supercandidate.py")
    ap.add_argument("--sfari-file",
                    default=os.path.join(
                        "data", "SFARI-Gene_genes_03-28-2024release_05-17-2024.csv"),
                    help="SFARI gene release CSV")
    ap.add_argument("--output-dir", default="results/main",
                    help="Where to write ndd_master_table.tsv and the derived inputs")
    args = ap.parse_args()

    build(args.data_dir, args.supercandidate, args.sfari_file, args.output_dir)
