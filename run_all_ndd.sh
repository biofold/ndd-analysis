#!/bin/bash
prog_dir=/Users/emidio/projects/ndd-analysis
cancer_dir=$prog_dir/cancer
cancer_out=$prog_dir/out_cancer

cd $prog_dir
mkdir -p $prog_dir/results
mkdir -p $prog_dir/out_cancer

source /Users/emidio/miniconda/etc/profile.d/conda.sh
conda activate ndd_analysis

echo "- Enrichment Analysis"
python3 run_ndd_analysis.py -c config.yaml

#Cancer analysis
echo "- Cancer Analysis"
python3 $prog_dir/utils/eqset.py $prog_dir/data/gene_set2.txt $prog_dir/data/cancer.txt 1 1 False >$cancer_dir/cancer_gs2.txt 
python3 $prog_dir/utils/eqset.py $prog_dir/data/gene_set1.txt $prog_dir/data/cancer.txt 1 1 False >$cancer_dir/cancer_gs1.txt 
python3 $prog_dir/utils/eqset.py $prog_dir/data/gene_set0.txt $prog_dir/data/cancer.txt 1 1 False >$cancer_dir/cancer_gs0.txt 
python3 $prog_dir/utils/difset.py $prog_dir/data/gene_set1.txt $cancer_dir/cancer_gs1.txt 1 1 >$cancer_dir/noncancer_gs1.txt 
python3 $prog_dir/utils/difset.py $prog_dir/data/gene_set2.txt $cancer_dir/cancer_gs2.txt 1 1 >$cancer_dir/noncancer_gs2.txt 
python3 $prog_dir/utils/difset.py $prog_dir/data/gene_set0.txt $cancer_dir/cancer_gs0.txt 1 1 >$cancer_dir/noncancer_gs0.txt 
echo "  1. Cancer Gene Enrichment"
python3 $prog_dir/scripts/1_enrichr_all.py $cancer_dir/cancer_gs2.txt $cancer_dir/noncancer_gs2.txt $prog_dir/data/gene_set2.txt $prog_dir/data/gene_all.txt --output_dir $cancer_out
python3 $prog_dir/scripts/1_enrichr_all.py $cancer_dir/cancer_gs1.txt $cancer_dir/noncancer_gs1.txt $prog_dir/data/gene_set1.txt $prog_dir/data/gene_all.txt --output_dir $cancer_out
python3 $prog_dir/scripts/1_enrichr_all.py $cancer_dir/cancer_gs0.txt $cancer_dir/noncancer_gs0.txt $prog_dir/data/gene_set0.txt $prog_dir/data/gene_all.txt --output_dir $cancer_out 

echo "  2. Aggregate cancer results"
cs1=`wc -l  $cancer_dir/cancer_gs1.txt |awk '{print $1}'`
ncs1=`wc -l  $cancer_dir/noncancer_gs1.txt |awk '{print $1}'`
s1=`wc -l $prog_dir/data/gene_set1.txt |awk '{print $1}'`

python3 $prog_dir/scripts/aggregate_pvals.py $cancer_out/cancer_gs1_GO_Biological_Process_2026.tsv $cancer_out/noncancer_gs1_GO_Biological_Process_2026.tsv $prog_dir/results/gene_set1_GO_Biological_Process_2026.tsv $cs1 $ncs1 $s1 > $cancer_out/aggregate_set1.txt

cs2=`wc -l  $cancer_dir/cancer_gs2.txt |awk '{print $1}'`
ncs2=`wc -l  $cancer_dir/noncancer_gs2.txt |awk '{print $1}'`
s2=`wc -l $prog_dir/data/gene_set2.txt |awk '{print $1}'`

python3 $prog_dir/scripts/aggregate_pvals.py $cancer_out/cancer_gs2_GO_Biological_Process_2026.tsv $cancer_out/noncancer_gs2_GO_Biological_Process_2026.tsv $prog_dir/results/gene_set2_GO_Biological_Process_2026.tsv $cs2 $ncs2 $s2 >$cancer_out/aggregate_set2.txt

echo "  3. Calculate fisher test"
python $prog_dir/utils/fisher_cols.py $cancer_out/aggregate_set1.txt 3,4 5,6 >$cancer_out/cancer_compara_set1.txt
python $prog_dir/utils/fisher_cols.py $cancer_out/aggregate_set2.txt 3,4 5,6 >$cancer_out/cancer_compara_set2.txt

echo "  4. BH correction"
python $prog_dir/utils/bh.py  $cancer_out/cancer_compara_set1.txt 17 |sort -gk 24   >$cancer_out/cancer_compara_set1.adj 
python $prog_dir/utils/bh.py  $cancer_out/cancer_compara_set2.txt 17 |sort -gk 24   >$cancer_out/cancer_compara_set2.adj

echo "  5. Make scatter plot"
python3 $prog_dir/scripts/scatter_cancer.py  <(awk '{if ($11<=0.01) print -log($9)/log(10),-log($10)/log(10),-log($24)/log(10)}' $cancer_out/cancer_compara_set1.adj) 2 1 3  --output $cancer_out/plot-compara-set1.png
python3 $prog_dir/scripts/scatter_cancer.py  <(awk '{if ($11<=0.01) print -log($9)/log(10),-log($10)/log(10),-log($24)/log(10)}' $cancer_out/cancer_compara_set2.adj) 2 1 3  --output $cancer_out/plot-compara-set2.png
