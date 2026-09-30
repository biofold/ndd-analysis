# utils/

Every script here is invoked by `run_ndd_analysis.py` (or, for
`data_extraction/`, run by hand when a source database is refreshed — see
`data/README.md`). Grouped by task rather than left flat:

- **`data_extraction/`** — pull per-gene tables out of raw third-party
  downloads (`data/raw/`) and into `data/*.tsv`: gnomAD pLI, ClinVar P/LP
  variants, HPO phenotype terms, ClinGen haploinsufficiency, plus
  `verify_gnomad_pli_source.py`, which cross-checks the dbNSFP-derived pLI
  against gnomAD's own release.

- **`enrichment/`** — build or post-process Enrichr/gseapy enrichment
  results: `build_goslim_gmt.py` builds the GO Slim libraries from GO's human annotations with the GO Consortium's map2slim (see data/README.md); `normalize_gmt_case.py` upper-cases the gene symbols of a .gmt library (Enrichr convention, used by every library in libs/; `--check` exits 1 if any library is not upper case). Per-gene data files keyed to the HGNC spine (e.g. data/mondo.txt) keep HGNC case;
  `enrichr_odds_ratio_ci.py` appends a 95% CI to the Odds Ratio column
  (wired into `scripts/1_enrichr_all.py` — every enrichment table the
  pipeline produces gets `OR_CI_low`/`OR_CI_high` automatically).

- **`moe_validation/`** — independent validation of the MOE prioritization
  score against reviewer points R1.3/R2.9/R2.10: `build_master_table.py`
  assembles the merged per-gene table these run on;
  `moe_clinvar_validation.py` and `moe_pli_validation.py` test it against
  ClinVar and gnomAD constraint; `moe_annotation_bias_baseline.py` tests it
  against annotation-volume baselines; `mondo_hpo_by_moe.py` compares MONDO
  and HPO term coverage across MOE tiers; `moe_component_analysis.py`
  quantifies redundancy between the five MOE components (phi, Jaccard, PCA),
  compares leave-one-component-out and HPO-augmented scores by DeLong's test,
  reports Youden J / F1 per MOE cutoff, ClinGen haploinsufficiency by tier and
  MOE discrimination within PubMed-count quartiles (Tables S15-S20).

- **`plotting/`** — generic statistical plotting primitives shared across
  pipeline steps (violin, bar, pairwise KS/Fisher-test matrices, ROC/PR,
  Venn, UpSet, scatter). Take a delimited file plus 1-based column indices;
  none of them are specific to any one dataset.

- **`table_generation/`** — shell scripts that compute the manuscript's
  Table 2 and Tables S1-S3/compara from pipeline outputs, plus
  `tsv2excel.py`, which packages TSVs into the supplementary `.xlsx` files,
  and `source_contribution.py`, which checks the Table 1 classification rule
  against the pipeline classes, counts genes per rule and flags per gene,
  reclassifies with each resource removed, lists SFARI-only candidates and
  draws the classification decision tree (Tables S21-S25; step 6, Part E.4).

- **`cancer_comparison/`** — the cancer-vs-non-cancer GO-term comparison
  pipeline. `run_cancer_comparison.sh` calls `aggregate_pvals.py` and
  `fisher_bh.py` by path relative to its own location, so the three stay
  together.

- **`file_ops/`** — generic gene-list/set manipulation shared by several
  pipeline steps (`set_operations.sh`: intersect/diff/union on keyed files;
  `lib2list.sh`: flatten a `.gmt` library to a gene list).

No script in one subfolder imports from another subfolder or from a sibling
script outside `cancer_comparison/` — they are called only from
`run_ndd_analysis.py` (or, for `data_extraction/`, from the command line),
each with an absolute path built via `get_absolute_path("utils/<subfolder>/<script>")`.
