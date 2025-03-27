#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

def get_absolute_path(relative_path):
    """Convert relative path to absolute path based on script location"""
    return (Path(__file__).parent / relative_path).resolve()

def run_command(command, step_name=None):
    """
    Run a command with proper error handling and output capture
    
    Args:
        command (list): Command to run as list of strings
        step_name (str): Optional description of the step
    
    Returns:
        subprocess.CompletedProcess: The completed process object
    """
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    # Print the command to stderr
    print(f"Running: {' '.join(command)}", file=sys.stderr)
    
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        return result
    except subprocess.CalledProcessError as e:
        print("\nCommand failed with error:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        raise
    except Exception as e:
        print(f"\nUnexpected error running command: {e}", file=sys.stderr)
        raise

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
    
    print("=== Starting NDD Analysis Pipeline ===", file=sys.stderr)
    
    # Get absolute paths to script files (now in scripts directory)
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    script2 = get_absolute_path("scripts/2_supercandidate.py")
    script3 = get_absolute_path("scripts/3_count_supercandidate.py")
    script4 = get_absolute_path("scripts/4_compare_subsets.py")

    try:
        # Step 1: Enrichment analysis
        run_command([
            "python", str(script1),
            str(gene_files["set2"]),
            str(gene_files["set1"]),
            str(gene_files["set0"]),
            str(gene_files["background"]),
            "--output_dir", str(results_dir)
        ], step_name="Running enrichment analysis")
        
        # Step 2: Supercandidate identification
        run_command([
            "python", str(script2),
            str(gene_files["set1"]),
            str(gene_files["set2"]),
            "--output_dir", str(results_dir)
        ], step_name="Identifying supercandidate genes")
        
        # Step 3: Score distribution analysis
        mondo_enrichment = results_dir / "gene_set2_MONDO_GROUPS_2025.tsv"
        run_command([
            "python", str(script3),
            str(results_dir / "supercandidate.tsv"),
            str(gene_files["set1"]),
            str(mondo_enrichment),
            str(gmt_file),
            "--plot", str(results_dir / "dist_mondo_supercandidate.png"),
            "--output", str(results_dir / "dist_mondo_supercandidate.txt")
        ], step_name="Analyzing score distributions")
        
        # Step 4: Subset comparison
        run_command([
            "python", str(script4),
            str(results_dir / "dist_mondo_supercandidate.txt"),
            "2", "3",
            "--label_col", "1",
            "--matrix", str(results_dir / "mondo_supercandidate_matrix.txt"),
            "--dendrogram", str(results_dir / "mondo_supercandidate_fisher.png"),
            "--output_file", str(results_dir / "mondo_supercandidate_fisher.txt")
        ], step_name="Comparing subsets")
        
        print("\n=== Pipeline completed successfully ===", file=sys.stderr)
        print(f"Results saved to: {results_dir}", file=sys.stderr)
    
    except subprocess.CalledProcessError:
        print("\n=== Pipeline failed ===", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n=== Pipeline failed with unexpected error: {e} ===", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    run_pipeline()
