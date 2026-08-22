#!/usr/bin/env python3
import os
import argparse
import subprocess
import sys
import yaml
from pathlib import Path

def get_absolute_path(relative_path):
    """Convert relative path to absolute path based on script location"""
    return (Path(__file__).parent / relative_path).resolve()

def get_conda_python(env_name="ndd_analysis"):
    """
    Get the Python executable path for a conda environment.
    
    Args:
        env_name (str): Name of the conda environment
    
    Returns:
        str: Path to the Python executable in the conda environment
    """
    try:
        # Try to get the Python path from conda
        result = subprocess.run(
            ["conda", "run", "-n", env_name, "which", "python"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        python_path = result.stdout.strip()
        return python_path
    except subprocess.CalledProcessError:
        # Fallback: try to construct the path manually
        conda_prefix = os.environ.get("CONDA_PREFIX", None)
        if conda_prefix:
            # If we're already in a conda environment, use its parent
            conda_base = os.path.dirname(conda_prefix)
        else:
            # Default conda installation paths
            possible_bases = [
                os.path.expanduser("~/miniconda3"),
                os.path.expanduser("~/anaconda3"),
                "/opt/miniconda3",
                "/opt/anaconda3"
            ]
            conda_base = None
            for base in possible_bases:
                if os.path.exists(base):
                    conda_base = base
                    break
        
        if conda_base:
            python_path = os.path.join(conda_base, "envs", env_name, "bin", "python")
            if os.path.exists(python_path):
                return python_path
        
        # If all else fails, use system python
        sys.stderr.write(f"Warning: Could not find conda environment '{env_name}', using system python\n")
        return "python"

def run_command(command, step_name=None, env_name="ndd_analysis"):
    """
    Run a command with proper error handling and output capture.
    Uses the specified conda environment.
    
    Args:
        command (list): Command to run as list of strings
        step_name (str): Optional description of the step
        env_name (str): Name of the conda environment to use
    
    Returns:
        subprocess.CompletedProcess: The completed process object
    """
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    # If the first element is 'python', replace it with the conda environment's python
    if command[0] in ["python", "python3"]:
        command[0] = get_conda_python(env_name)
    
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

def run_command_conda(command, step_name=None, env_name="ndd_analysis"):
    """
    Run a command using conda run to ensure it executes in the correct environment.
    
    Args:
        command (list): Command to run as list of strings (without 'python')
        step_name (str): Optional description of the step
        env_name (str): Name of the conda environment to use
    
    Returns:
        subprocess.CompletedProcess: The completed process object
    """
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    # Build conda run command
    conda_command = ["conda", "run", "-n", env_name] + command
    
    # Print the command to stderr
    print(f"Running: {' '.join(conda_command)}", file=sys.stderr)
    print(f"(Using conda environment: {env_name})", file=sys.stderr)
    
    try:
        result = subprocess.run(
            conda_command,
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

def load_yaml_config(config_path):
    """Load configuration from YAML file"""
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def check_conda_env(env_name="ndd_analysis"):
    """
    Check if the conda environment exists.
    
    Args:
        env_name (str): Name of the conda environment
    
    Returns:
        bool: True if environment exists, False otherwise
    """
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        # Check if environment is in the list
        for line in result.stdout.split('\n'):
            if line.strip() and not line.startswith('#'):
                if line.split()[0] == env_name:
                    return True
        return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.stderr.write("Warning: conda not found or not accessible\n")
        return False

def run_pipeline(data_dir=None, lib_dir=None, output_dir=None, 
                 set0=None, set1=None, set2=None, background=None,
                 config_file=None, libraries_enrichment=None, 
                 libraries_supercandidate=None, conda_env="ndd_analysis"):
    
    # Load YAML config if provided
    if config_file:
        print(f"Loading configuration from: {config_file}", file=sys.stderr)
        config = load_yaml_config(config_file)
        
        # YAML config overrides defaults but command-line args override YAML
        if not data_dir and 'data_dir' in config:
            data_dir = config['data_dir']
        if not lib_dir and 'lib_dir' in config:
            lib_dir = config['lib_dir']
        if not output_dir and 'output_dir' in config:
            output_dir = config['output_dir']
        if not set0 and 'set0' in config:
            set0 = config['set0']
        if not set1 and 'set1' in config:
            set1 = config['set1']
        if not set2 and 'set2' in config:
            set2 = config['set2']
        if not background and 'background' in config:
            background = config['background']
        if not libraries_enrichment and 'libraries_enrichment' in config:
            libraries_enrichment = config['libraries_enrichment']
        if not libraries_supercandidate and 'libraries_supercandidate' in config:
            libraries_supercandidate = config['libraries_supercandidate']
        if 'conda_env' in config:
            conda_env = config['conda_env']
    
    # Check conda environment
    print(f"=== Checking conda environment: {conda_env} ===", file=sys.stderr)
    if check_conda_env(conda_env):
        print(f"  ✓ Conda environment '{conda_env}' found", file=sys.stderr)
    else:
        print(f"  ⚠ Warning: Conda environment '{conda_env}' not found", file=sys.stderr)
        print(f"  Will attempt to use system python (this may fail if dependencies are missing)", file=sys.stderr)
    print(file=sys.stderr)
    
    # Set default libraries if not provided
    if not libraries_enrichment:
        libraries_enrichment = [
            "GO_Biological_Process_2026",
            "GO_Biological_Process_Cancer_2026",
            "GO_Cellular_Component_2026",
            "GO_Molecular_Function_2026",
            "KEGG_2021_Human",
            "Reactome_Pathways_2024",
            "SynGO_2024",
            "SynGO_BP_2024",
            "SynGO_CC_2024",
            "MONDO_2026",
            "MONDO_GROUPS_2026"
        ]
    
    if not libraries_supercandidate:
        libraries_supercandidate = [
            'GO_Biological_Process_2026', 
            'GO_Cellular_Component_2026',
            'GO_Molecular_Function_2026', 
            'KEGG_2021_Human', 
            'Reactome_Pathways_2024'
        ]
    
    # Set up absolute paths
    base_dir = Path(__file__).parent
    
    # Use specified directories or defaults
    if data_dir:
        data_dir = Path(data_dir).resolve()
    else:
        data_dir = get_absolute_path("data")
    
    if lib_dir:
        lib_dir = Path(lib_dir).resolve()
    else:
        lib_dir = get_absolute_path("libs")
    
    if output_dir:
        results_dir = Path(output_dir).resolve()
    else:
        results_dir = get_absolute_path("results")
    
    # Create results directory if it doesn't exist
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Define file paths for gene sets
    gene_files = {
        "set0": Path(set0).resolve() if set0 else data_dir / "gene_set0.txt",
        "set1": Path(set1).resolve() if set1 else data_dir / "gene_set1.txt",
        "set2": Path(set2).resolve() if set2 else data_dir / "gene_set2.txt",
        "background": Path(background).resolve() if background else data_dir / "gene_all.txt"
    }
    
    gmt_file = lib_dir / "MONDO_GROUPS_2026.gmt"
    
    # Check if all input files exist
    print("=== Configuration ===", file=sys.stderr)
    print(f"Conda environment: {conda_env}", file=sys.stderr)
    print(f"Data directory: {data_dir}", file=sys.stderr)
    print(f"Library directory: {lib_dir}", file=sys.stderr)
    print(f"Output directory: {results_dir}", file=sys.stderr)
    if config_file:
        print(f"Config file: {config_file}", file=sys.stderr)
    print(file=sys.stderr)
    
    print("=== Libraries Configuration ===", file=sys.stderr)
    print(f"Enrichment analysis libraries ({len(libraries_enrichment)}):", file=sys.stderr)
    for lib in libraries_enrichment:
        print(f"  • {lib}", file=sys.stderr)
    print(f"\nSupercandidate libraries ({len(libraries_supercandidate)}):", file=sys.stderr)
    for lib in libraries_supercandidate:
        print(f"  • {lib}", file=sys.stderr)
    print(file=sys.stderr)
    
    print("=== Checking input files ===", file=sys.stderr)
    
    for name, path in gene_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing input file: {path}")
        print(f"  ✓ Found {name}: {path}", file=sys.stderr)
    
    if not gmt_file.exists():
        raise FileNotFoundError(f"Missing GMT file: {gmt_file}")
    print(f"  ✓ Found GMT file: {gmt_file}", file=sys.stderr)
    
    print(f"\n=== Starting NDD Analysis Pipeline ===", file=sys.stderr)
    
    # Get absolute paths to script files (now in scripts directory)
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    script2 = get_absolute_path("scripts/2_supercandidate.py")
    script3 = get_absolute_path("scripts/3_count_supercandidate.py")
    script4 = get_absolute_path("scripts/4_compare_subsets.py")
    
    # Check if scripts exist
    for script in [script1, script2, script3, script4]:
        if not script.exists():
            raise FileNotFoundError(f"Missing script: {script}")

    # Convert library lists to comma-separated strings
    enrichment_libs_str = ",".join(libraries_enrichment)
    supercandidate_libs_str = ",".join(libraries_supercandidate)
    
    try:
        # Step 1: Enrichment analysis
        run_command_conda([
            "python", str(script1),
            str(gene_files["set2"]),
            str(gene_files["set1"]),
            str(gene_files["set0"]),
            str(gene_files["background"]),
            "--output_dir", str(results_dir),
            "--libraries", enrichment_libs_str
        ], step_name="Running enrichment analysis", env_name=conda_env)
        
        # Step 2: Supercandidate identification
        run_command_conda([
            "python", str(script2),
            str(gene_files["set1"]),
            str(gene_files["set2"]),
            "--output_dir", str(results_dir),
            "--libraries", supercandidate_libs_str
        ], step_name="Identifying supercandidate genes", env_name=conda_env)
        
        # Step 3: Score distribution analysis
        mondo_enrichment = results_dir / "gene_set2_MONDO_GROUPS_2026.tsv"
        run_command_conda([
            "python", str(script3),
            str(results_dir / "supercandidate.tsv"),
            str(gene_files["set1"]),
            str(mondo_enrichment),
            str(gmt_file),
            "-p", "1.00",
            "--plot", str(results_dir / "dist_mondo_supercandidate.png"),
            "--output", str(results_dir / "dist_mondo_supercandidate.txt")
        ], step_name="Analyzing score distributions", env_name=conda_env)
        
        # Step 4: Subset comparison
        run_command_conda([
            "python", str(script4),
            str(results_dir / "dist_mondo_supercandidate.txt"),
            "2", "3",
            "--label_col", "1",
            "--matrix", str(results_dir / "mondo_supercandidate_matrix.txt"),
            "--dendrogram", str(results_dir / "mondo_supercandidate_fisher.png"),
            "--output_file", str(results_dir / "mondo_supercandidate_fisher.txt")
        ], step_name="Comparing subsets", env_name=conda_env)
        
        print("\n=== Pipeline completed successfully ===", file=sys.stderr)
        print(f"Results saved to: {results_dir}", file=sys.stderr)
        print("\nOutput files:", file=sys.stderr)
        for file in sorted(results_dir.glob("*")):
            if file.is_file():
                print(f"  • {file.name}", file=sys.stderr)
    
    except subprocess.CalledProcessError:
        print("\n=== Pipeline failed ===", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n=== Pipeline failed with unexpected error: {e} ===", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="NDD Analysis Pipeline - Run enrichment analysis and supercandidate identification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default directories and libraries with ndd_analysis conda environment
  %(prog)s
  
  # Use a YAML configuration file
  %(prog)s -c config.yaml
  
  # Specify a different conda environment
  %(prog)s --conda_env my_env
  
  # Specify custom libraries via command line
  %(prog)s --libraries_enrichment "GO_BP_2026,KEGG_2021" --libraries_supercandidate "GO_BP_2026"

YAML Configuration File Format:
  data_dir: /Users/emidio/projects/ndd-analysis/data
  lib_dir: /Users/emidio/projects/ndd-analysis/libs
  output_dir: /Users/emidio/projects/ndd-analysis/results
  conda_env: ndd_analysis
  set0: gene_set0.txt
  set1: gene_set1.txt
  set2: gene_set2.txt
  background: gene_all.txt
  libraries_enrichment:
    - GO_Biological_Process_2026
    - GO_Cellular_Component_2026
    - GO_Molecular_Function_2026
    - KEGG_2021_Human
    - Reactome_Pathways_2024
  libraries_supercandidate:
    - GO_Biological_Process_2026
    - GO_Cellular_Component_2026
    - GO_Molecular_Function_2026
    - KEGG_2021_Human
    - Reactome_Pathways_2024
        """
    )
    
    # Directory options
    parser.add_argument(
        "-d", "--data_dir",
        type=str,
        default=None,
        help="Directory containing input gene files (default: 'data' in script directory)"
    )
    parser.add_argument(
        "-l", "--lib_dir",
        type=str,
        default=None,
        help="Directory containing GMT library files (default: 'libs' in script directory)"
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default=None,
        help="Directory for output results (default: 'results' in script directory)"
    )
    
    # Gene set options
    parser.add_argument(
        "-s0", "--set0",
        type=str,
        default=None,
        help="Path to gene set 0 file (default: data/gene_set0.txt)"
    )
    parser.add_argument(
        "-s1", "--set1",
        type=str,
        default=None,
        help="Path to gene set 1 file (default: data/gene_set1.txt)"
    )
    parser.add_argument(
        "-s2", "--set2",
        type=str,
        default=None,
        help="Path to gene set 2 file (default: data/gene_set2.txt)"
    )
    parser.add_argument(
        "-b", "--background",
        type=str,
        default=None,
        help="Path to background gene file (default: data/gene_all.txt)"
    )
    
    # Library options
    parser.add_argument(
        "--libraries_enrichment",
        type=str,
        default=None,
        help="Comma-separated list of libraries for enrichment analysis"
    )
    parser.add_argument(
        "--libraries_supercandidate",
        type=str,
        default=None,
        help="Comma-separated list of libraries for supercandidate identification"
    )
    
    # Conda environment option
    parser.add_argument(
        "--conda_env",
        type=str,
        default="ndd_analysis",
        help="Name of the conda environment to use (default: ndd_analysis)"
    )
    
    # YAML configuration file
    parser.add_argument(
        "-c", "--config",
        type=str,
        default=None,
        help="YAML configuration file (command-line args override YAML values)"
    )
    
    args = parser.parse_args()
    
    # Parse comma-separated library strings if provided via command line
    libraries_enrichment = None
    libraries_supercandidate = None
    
    if args.libraries_enrichment:
        libraries_enrichment = [lib.strip() for lib in args.libraries_enrichment.split(",")]
    if args.libraries_supercandidate:
        libraries_supercandidate = [lib.strip() for lib in args.libraries_supercandidate.split(",")]
    
    run_pipeline(
        data_dir=args.data_dir,
        lib_dir=args.lib_dir,
        output_dir=args.output_dir,
        set0=args.set0,
        set1=args.set1,
        set2=args.set2,
        background=args.background,
        config_file=args.config,
        libraries_enrichment=libraries_enrichment,
        libraries_supercandidate=libraries_supercandidate,
        conda_env=args.conda_env
    )
