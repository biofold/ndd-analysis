#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

def get_absolute_path(relative_path):
    """Convert relative path to absolute path based on script location"""
    return (Path(__file__).parent / relative_path).resolve()

def run_pipeline():
    # Set up absolute paths
    base_dir = Path(__file__).parent
    data_dir = get_absolute_path("data")
    lib_dir = get_absolute_path("libs")
    results_dir = get_absolute_path("results")
    
    # Create results directory if it doesn't exist
    results_dir.mkdir(exist_ok=True)
    
    # Define file paths
    gene_files = {
        "set0": data_dir / "gene_set0.txt",
        "set1": data_dir / "gene_set1.txt",
        "set2": data_dir / "gene_set2.txt",
        "background": data_dir / "gene_all.txt"
    }
    gmt_file = lib_dir / "MONDO_GROUPS_2025.gmt"
    
    # Check if all input files exist
    for name, path in gene_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing input file: {path}")
    
    print("=== Starting NDD Analysis Pipeline ===")
    
    # Get absolute paths to script files (now in scripts directory)
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    script2 = get_absolute_path("scripts/2_supercandidate.py")
    script3 = get_absolute_path("scripts/3_count_supercandidate.py")
    script4 = get_absolute_path("scripts/4_compare_subsets.py")

    # Step 1: Enrichment analysis
    print("\nRunning enrichment analysis...")
    subprocess.run([
        "python", str(script1),
        str(gene_files["set2"]),
        str(gene_files["set1"]),
        str(gene_files["set0"]),
        str(gene_files["background"]),
        "--output_dir", str(results_dir)
    ], check=True)
    
    # Step 2: Supercandidate identification
    print("\nIdentifying supercandidate genes...")
    subprocess.run([
        "python", str(script2),
        str(gene_files["set1"]),
        str(gene_files["set2"]),
        "--output_dir", str(results_dir)
    ], check=True)
    
    # Step 3: Score distribution analysis
    print("\nAnalyzing score distributions...")
    mondo_enrichment = results_dir / "gene_set2_MONDO_GROUPS_2025.tsv"
    subprocess.run([
        "python", str(script3),
        str(results_dir / "supercandidate.tsv"),
        str(gene_files["set1"]),
        str(mondo_enrichment),
        str(gmt_file),
        "--plot", str(results_dir / "dist_mondo_supercandidate.png"),
        "--output", str(results_dir / "dist_mondo_supercandidate.txt")
    ], check=True)
    
    # Step 4: Subset comparison
    print("\nComparing subsets...")
    subprocess.run([
        "python", str(script4),
        str(results_dir / "dist_mondo_supercandidate.txt"),
        "2", "3",
        "--label_col", "1",
        "--matrix", str(results_dir / "mondo_supercandidate_matrix.txt"),
        "--dendrogram", str(results_dir / "mondo_supercandidate_fisher.png"),
        "--output_file", str(results_dir / "mondo_supercandidate_fisher.txt")
    ], check=True)
    
    print("\n=== Pipeline completed successfully ===")
    print(f"Results saved to: {results_dir}")

if __name__ == "__main__":
    run_pipeline()
