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
   Source of `gene_gnomad_pli.tsv` via `utils/extract_gnomad_pli.py`.
2. `data/raw/variant_summary_20260924.txt.gz` -- NCBI ClinVar's weekly
   variant_summary release (downloaded 2026-09-28; NCBI's
   Last-Modified header on the file itself reads 2026-09-24, ~430MB).
   Source of `clinvar_plp_variants.tsv` / `clinvar_plp_gene_counts.tsv`
   via `utils/extract_clinvar_plp.py`.

Derived per-gene/per-variant tables (regenerate with the commands below):

3. `data/gene_gnomad_pli.tsv` -- one row per dbNSFP gene (40,225 rows):
   Gene, gnomAD_pLI, RVIS_percentile_ExAC, HIPred, GHIS,
   ClinGen_Haploinsufficiency_Score.
   `python utils/extract_gnomad_pli.py --dbnsfp-gene-file data/raw/dbNSFP5.2_gene.gz --output data/gene_gnomad_pli.tsv`
4. `data/clinvar_plp_variants.tsv` -- 366,744 ClinVar variant records
   classified Pathogenic/Likely-pathogenic by the text-based rule
   (see utils/extract_clinvar_plp.py docstring for why ClinSigSimple
   alone is not used), GRCh38 only.
5. `data/clinvar_plp_gene_counts.tsv` -- per-gene rollup of (4):
   Gene, n_pathogenic_likely_pathogenic, n_plp_ge1star.
   `python utils/extract_clinvar_plp.py --input data/raw/variant_summary_20260924.txt.gz --output-prefix data/clinvar_plp`
   (or pass `--url` with no value to fetch the current release fresh
   from NCBI instead of using the checked-in copy).

Both derived tables were cross-checked against ndd_master_table.tsv and
match exactly (0/19,354 gene-level mismatches on the ClinVar counts;
exact value match on gnomAD_pLI for spot-checked genes).

## Useful links
1. DDG2P   https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz \
           https://ftp.ebi.ac.uk/pub/databases/gene2phenotype
2. SFARI   https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes \
           https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=human-gene-scores
3. GeneTrek https://genetrek.pasteur.fr/downloadAllData?filetype=tsv
4. SynNDD  https://sysndd.dbmr.unibe.ch/Genes
5. DBD     https://dbd.geisingeradmi.org/downloads/DBD-Genes-Full-Data.csv
