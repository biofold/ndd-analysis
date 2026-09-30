## Collection of tables of NDD associated Genes

1. SFARI Database Authism related genes			SFARI-Gene_genes_03-28-2024release_05-17-2024.csv
2. Orphanet developmental disorders			develop_genename.txt
3. Orphanet neurological disorders			neuro_genename.txt
4. Data from publication PMID:37541188			list-ndd-37541188.txt
5. Combination of datasets 1-4				list-ndd-4sets.txt
6. HC genes from PMID:33932580				genelist_hc_pmid-33932580-ndd.txt	
7. Candidate genese from PMID:33932580			candidate_genes_pmid-33932580-ndd.txt
8. Full list genes from PMID:33932580			pmid-33932580-ndd.txt

The Orphanet genes are retrieved from UP000005640_9606.idmapping.tar.gz release Jan 2024.


## gnomAD pLI and ClinVar P/LP status (R1.3 independent validation)

Raw source files (tracked via Git LFS, see .gitattributes) live in `data/raw/`:

1. `data/raw/dbNSFP5.2_gene.gz` -- dbNSFP5.2 gene-level annotation table
   (downloaded 2026-09-27 from https://dbnsfp.s3.amazonaws.com/, ~22MB).
   Source of `gene_gnomad_pli.tsv` via `utils/data_extraction/extract_gnomad_pli.py`.
2. `data/raw/variant_summary_20260924.txt.gz` -- NCBI ClinVar's weekly
   variant_summary release (downloaded 2026-09-28; NCBI's
   Last-Modified header on the file itself reads 2026-09-24, ~430MB).
   Source of `clinvar_plp_variants.tsv` / `clinvar_plp_gene_counts.tsv`
   via `utils/data_extraction/extract_clinvar_plp.py`.

Derived per-gene/per-variant tables (regenerate with the commands below):

3. `data/gene_gnomad_pli.tsv` -- one row per dbNSFP gene (40,225 rows):
   Gene, gnomAD_pLI, RVIS_percentile_ExAC, HIPred, GHIS,
   ClinGen_Haploinsufficiency_Score.
   `python utils/data_extraction/extract_gnomad_pli.py --dbnsfp-gene-file data/raw/dbNSFP5.2_gene.gz --output data/gene_gnomad_pli.tsv`
4. `data/clinvar_plp_variants.tsv` -- 366,744 ClinVar variant records
   classified Pathogenic/Likely-pathogenic by the text-based rule
   (see utils/data_extraction/extract_clinvar_plp.py docstring for why ClinSigSimple
   alone is not used), GRCh38 only.
5. `data/clinvar_plp_gene_counts.tsv` -- per-gene rollup of (4):
   Gene, n_pathogenic_likely_pathogenic, n_plp_ge1star.
   `python utils/data_extraction/extract_clinvar_plp.py --input data/raw/variant_summary_20260924.txt.gz --output-prefix data/clinvar_plp`
   (or pass `--url` with no value to fetch the current release fresh
   from NCBI instead of using the checked-in copy).

Both derived tables were cross-checked against ndd_master_table.tsv and
match exactly (0/19,354 gene-level mismatches on the ClinVar counts;
exact value match on gnomAD_pLI for spot-checked genes).

### Provenance check: dbNSFP's gnomAD_pLI vs gnomAD's own release

`gnomAD_pLI` (as pulled from dbNSFP above) is specifically gnomAD's
**v2.1.1** pLI value. `utils/data_extraction/verify_gnomad_pli_source.py` fetches
gnomAD's own gene constraint table directly (from gnomAD's AWS Open
Data bucket, arn:aws:s3:::gnomad-public-us-east-1) and compares:

- `data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` (same-version
  check): of 17,268 genes comparable (non-missing on both sides),
  99.86% (17,243) match dbNSFP's value exactly; ALL 25 non-matches are
  genes with more than one transcript row in gnomAD's own v2.1.1 file
  (no unambiguous canonical-transcript flag in that release) -- i.e.
  dbNSFP is confirmed to be a faithful pass-through of gnomAD v2.1.1,
  not an independent recomputation, with zero unexplained discrepancy.
- `data/raw/gnomad.v4.1.constraint_metrics.tsv.gz` (cross-version
  check, MANE-Select-filtered to one row/gene): only 28.9% exact match
  against v2.1.1, but Spearman rho=0.834 (n=15,875) -- gnomAD v4.1 used
  a ~6x larger cohort and explicitly recalculated pLI (recommending
  LOEUF as the primary constraint metric going forward), so the two
  versions ranking genes similarly while disagreeing on point
  estimates is expected, not a bug. See
  `data/gnomad_pli_v41_vs_v211_comparison.tsv` for the full per-gene
  comparison. ndd_master_table.tsv's gnomAD_pLI column remains the
  v2.1.1 value throughout this pipeline; switching to v4.1 would
  require re-baselining every downstream analysis, not just adding a
  check.
  `python utils/data_extraction/verify_gnomad_pli_source.py --version 4.1 --input data/raw/gnomad.v4.1.constraint_metrics.tsv.gz --dbnsfp-pli-table data/gene_gnomad_pli.tsv`

## gnomAD v2.1.1 constraint from gnomAD's own release (pLI, LOEUF)

The single gnomAD extractor. gnomAD v4.x is not used anywhere in iNDDx.

- Input: `data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` (19,704 genes, one row per
  gene, canonical transcript, GRCh37 / GENCODE v19; covers chrX and chrY).
- `python utils/data_extraction/extract_gnomad_constraint.py --input data/raw/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz --output data/gene_gnomad_v2_constraint.tsv`
- Output `data/gene_gnomad_v2_constraint.tsv`: `Ensembl_ID, Gene, transcript, gnomad_v2_pli,
  gnomad_v2_loeuf, gnomad_v2_oe_lof, gnomad_v2_loeuf_decile, gnomad_v2_constraint_flag`.
  19,197 genes have pLI and LOEUF (always both or neither; the script refuses otherwise).
  LOEUF is `oe_lof_upper`; lower = more constrained; decile 0 = most constrained.
  `constraint_flag` (e.g. `no_exp_lof`, `syn_outlier`) is gnomAD's reliability warning,
  emitted, not used to drop genes.
- Join to iNDDx genes on **Ensembl gene ID**, then on a symbol only when it names exactly
  one gene_id (98 v2.1.1 symbols map to more than one). This covers 18,311 iNDDx genes.

**Why not dbNSFP's pLI (`gene_gnomad_pli.tsv`, `extract_gnomad_pli.py`).** dbNSFP 5.2 passes
gnomAD v2.1.1 pLI through unchanged -- 99.89% exact agreement on the 17,056 genes both score
-- but it matched gnomAD by gene **symbol**, which fails two ways:

| failure | genes | cause | example |
|---|---|---|---|
| another gene's value | 18 | symbol has two gnomAD rows; dbNSFP took the other | *TUBB3* 0.000 vs 0.968, *PI4K2A* 0.000 vs 0.804, *MDGA2* 0.665 vs 0.999 |
| empty although gnomAD scores it | 1,385 | 1,382 are HGNC renames since gnomAD's 2018 naming | *AARS1*/AARS, *ABRAXAS1*/FAM175A, *ACKR1*/DARC |

Among iNDDx genes the empty cases hit 117 high-confidence, 287 candidate and 742
no-evidence genes. dbNSFP also carries no LOEUF. The same symbol-matching gap affects
dbNSFP's RVIS and GDI (99.6% / 99.8% missing among renamed genes vs 12.1% / 5.9% otherwise),
which have no reachable primary source to repair from.

**Status of the published analyses.** The pLI analyses (Figure 5, the KS matrix, the
MOE-tier comparisons, `utils/moe_validation/`) were computed on `gene_gnomad_pli.tsv`
(dbNSFP) and still read it. Switch them to `gene_gnomad_v2_constraint.tsv` deliberately
when recalculating; nothing has been switched yet.

The v4.1 raw file, `verify_gnomad_pli_source.py` and
`gnomad_pli_v41_vs_v211_comparison.tsv` are kept only as the record of the R1.3 provenance
check (v2.1.1 same-version 99.86%; v4.1 cross-version Spearman rho 0.834).

## HPO and ClinGen haploinsufficiency (moved off dbNSFP to primary sources)

ndd_master_table.tsv's HPO_id/HPO_name and ClinGen_Haploinsufficiency_Score
columns were originally dbNSFP5.2 re-exports; both are now sourced
directly from their own primary databases instead:

- `data/raw/genes_to_phenotype_20260928.txt.gz` -- HPO's own
  gene-to-phenotype annotation file (downloaded 2026-09-28 from
  https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt).
  `utils/data_extraction/extract_hpo_terms.py` dedups by (gene_symbol, hpo_id) -- 18.0%
  of raw rows are duplicates via >1 source disease -- and upper-cases
  gene symbols to HGNC convention. Result: `data/gene_hpo_terms.tsv`.
  SCN2A: 208 terms direct vs 162 from dbNSFP's stale snapshot.
  `python utils/data_extraction/extract_hpo_terms.py --genes-to-phenotype-file data/raw/genes_to_phenotype_20260928.txt.gz --output data/gene_hpo_terms.tsv`
- `data/raw/clingen_dosage_sensitivity_20260928.tsv` -- ClinGen's own
  dosage-sensitivity curation export (downloaded 2026-09-28 from
  https://search.clinicalgenome.org/kb/gene-dosage/download).
  `utils/data_extraction/extract_clingen_hi.py` maps ClinGen's text evidence categories
  to the standard 0-3/30/40 numeric scale (verified empirically against
  dbNSFP: 99.2% agreement across 1,559 comparable genes; the 0.8% that
  differ are genes ClinGen has re-curated since dbNSFP's snapshot).
  Result: `data/gene_clingen_hi.tsv`.
  `python utils/data_extraction/extract_clingen_hi.py --clingen-dosage-file data/raw/clingen_dosage_sensitivity_20260928.tsv --output data/gene_clingen_hi.tsv`

Both were deployed to the live master table and MongoDB on 2026-09-28
(verified live: SCN2A now shows 208 HPO terms and ClinGen HI score 3.0
via https://data.biofold.org/inddx/api/gene/SCN2A).

Still sourced from dbNSFP5.2 (`data/raw/dbNSFP5.2_gene.gz`, see
extract_gnomad_pli.py above): gnomAD_pLI, RVIS_percentile_ExAC, HIPred,
GHIS. RVIS's primary source (genic-intolerance.org) now redirects to
an unrelated commercial site (domain squatted/abandoned) and GDI's
primary source (hgidsoft.rockefeller.edu) times out on TLS handshake
(server unreachable) -- neither is currently extractable directly.
gnomAD_pLI has a working direct source (see the provenance-check
section above) but has not yet been swapped in as the master table's
live source the way HPO/ClinGen were.

## GO Slim libraries (R1.2 redundancy reduction)

Enrichr provides no GO Slim library, so `libs/GOslim_{Biological_Process,
Cellular_Component,Molecular_Function}_2026.gmt` are built from GO's own files
with the GO Consortium's map2slim (OWLTools):

1. `data/raw/go-basic_2026-07-26.obo.gz` -- GO release 2026-07-26
   (https://current.geneontology.org/ontology/go-basic.obo, downloaded 2026-09-29).
2. `data/raw/goslim_generic_2026-07-26.obo.gz` -- generic GO Slim, 140 terms, same
   release (https://current.geneontology.org/ontology/subsets/goslim_generic.obo),
   kept for reference; the subset membership used is the `subset: goslim_generic`
   tag in (1).
3. `data/raw/goa_human_2026-05-21.gaf.gz` -- human GO annotations, date-generated
   2026-05-21 (https://current.geneontology.org/annotations/goa_human.gaf.gz,
   downloaded 2026-09-29).

map2slim assigns each annotation to its most specific slim term(s); the script
then adds every gene to all slim terms above those (is_a + part_of), so each slim
gene set contains every gene annotated to it at any depth. NOT-qualified and ND
annotations are removed. Needs Java and OWLTools 2024-06-12
(https://github.com/owlcollab/owltools/releases/download/2024-06-12/owltools);
neither is needed to run the pipeline itself.

  `python utils/enrichment/build_goslim_gmt.py --owltools /path/to/owltools --go-obo data/raw/go-basic_2026-07-26.obo.gz --gaf data/raw/goa_human_2026-05-21.gaf.gz --outdir libs --report libs/GOslim_2026_build_report.tsv`

Result: BP 69, CC 25, MF 38 slim terms (counts per library in
`libs/GOslim_2026_build_report.tsv`).

## Useful links
1. DDG2P   https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz \
           https://ftp.ebi.ac.uk/pub/databases/gene2phenotype
2. SFARI   https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes \
           https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=human-gene-scores
3. GeneTrek https://genetrek.pasteur.fr/downloadAllData?filetype=tsv
4. SynNDD  https://sysndd.dbmr.unibe.ch/Genes
5. DBD     https://dbd.geisingeradmi.org/downloads/DBD-Genes-Full-Data.csv
