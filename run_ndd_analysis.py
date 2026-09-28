#!/usr/bin/env python3
import os
import argparse
import subprocess
import sys
import yaml
import shutil
from pathlib import Path

def get_absolute_path(relative_path):
    """Convert relative path to absolute path based on script location"""
    return (Path(__file__).parent / relative_path).resolve()

def get_conda_python(env_name="ndd_analysis"):
    """Get the Python executable path for a conda environment."""
    try:
        result = subprocess.run(
            ["conda", "run", "-n", env_name, "which", "python"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        python_path = result.stdout.strip()
        return python_path
    except (subprocess.CalledProcessError, FileNotFoundError):
        # `conda` may be absent from PATH (e.g. the pipeline is launched from an
        # already-activated environment, or from a container that only exposes the
        # env's bin/). Fall back to locating the environment prefix by hand.
        candidate_bases = []
        conda_prefix = os.environ.get("CONDA_PREFIX", None)
        if conda_prefix:
            # CONDA_PREFIX is either the base install (.../miniconda3) or a named
            # environment (.../miniconda3/envs/<name>); cover both.
            candidate_bases.append(conda_prefix)
            parent = os.path.dirname(conda_prefix)
            if os.path.basename(parent) == "envs":
                candidate_bases.append(os.path.dirname(parent))
        candidate_bases.extend([
            os.path.expanduser("~/miniconda3"),
            os.path.expanduser("~/anaconda3"),
            "/opt/miniconda3",
            "/opt/anaconda3"
        ])

        for base in candidate_bases:
            python_path = os.path.join(base, "envs", env_name, "bin", "python")
            if os.path.exists(python_path):
                return python_path

        # Last resort: if this script is itself running inside an interpreter that
        # has the dependencies, reuse it rather than an arbitrary `python` on PATH.
        sys.stderr.write(
            f"Warning: Could not find conda environment '{env_name}', "
            f"falling back to the current interpreter ({sys.executable})\n"
        )
        return sys.executable

def run_command(command, step_name=None, env_name="ndd_analysis"):
    """Run a command with proper error handling and output capture."""
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    if command[0] in ["python", "python3"]:
        command[0] = get_conda_python(env_name)
    
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

def run_command_conda(command, step_name=None, env_name="ndd_analysis"):
    """Run a command using conda run."""
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    if shutil.which("conda"):
        conda_command = ["conda", "run", "-n", env_name] + command
    else:
        # No `conda` launcher on PATH: run the command directly with the
        # environment's interpreter (resolved by get_conda_python).
        conda_command = list(command)
        if conda_command and conda_command[0] in ("python", "python3"):
            conda_command[0] = get_conda_python(env_name)

    print(f"Running: {' '.join(conda_command)}", file=sys.stderr)
    
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

def run_bash_command(command, step_name=None):
    """Run a bash command."""
    if step_name:
        print(f"\n{step_name}...", file=sys.stderr)
    
    print(f"Running: {command}", file=sys.stderr)
    
    try:
        result = subprocess.run(
            command,
            shell=True,
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

def load_yaml_config(config_path):
    """Load configuration from YAML file"""
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config

def check_conda_env(env_name="ndd_analysis"):
    """Check if the conda environment exists."""
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        for line in result.stdout.split('\n'):
            if line.strip() and not line.startswith('#'):
                if line.split()[0] == env_name:
                    return True
        return False
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.stderr.write("Warning: conda not found or not accessible\n")
        return False

def generate_gene_sets_from_file(combined_file, output_dir):
    """Generate gene set files from a combined gene file."""
    print(f"\nGenerating gene sets from combined file: {combined_file}", file=sys.stderr)
    
    combined_path = Path(combined_file).resolve()
    if not combined_path.exists():
        raise FileNotFoundError(f"Combined gene file not found: {combined_path}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    set0_file = output_dir / "gene_set0.txt"
    set1_file = output_dir / "gene_set1.txt"
    set2_file = output_dir / "gene_set2.txt"
    background_file = output_dir / "gene_all.txt"
    
    counts = {"no_evidence": 0, "candidate": 0, "curated": 0, "total": 0}
    
    with open(set0_file, 'w') as f0, \
         open(set1_file, 'w') as f1, \
         open(set2_file, 'w') as f2, \
         open(background_file, 'w') as fb:
        
        with open(combined_path, 'r') as fin:
            for line in fin:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 2:
                    continue
                
                gene_id = fields[0].strip()
                classification = fields[1].strip().lower()
                
                if not gene_id:
                    continue
                
                fb.write(f"{gene_id}\n")
                counts["total"] += 1
                
                if classification == "no_evidence":
                    f0.write(f"{gene_id}\n")
                    counts["no_evidence"] += 1
                elif classification == "candidate":
                    f1.write(f"{gene_id}\n")
                    counts["candidate"] += 1
                elif classification == "curated":
                    f2.write(f"{gene_id}\n")
                    counts["curated"] += 1
                else:
                    sys.stderr.write(f"Warning: Unknown classification '{classification}' for gene {gene_id}\n")
    
    print(f"  Generated gene sets:", file=sys.stderr)
    print(f"    • set0 (no_evidence): {counts['no_evidence']} genes", file=sys.stderr)
    print(f"    • set1 (candidate): {counts['candidate']} genes", file=sys.stderr)
    print(f"    • set2 (curated): {counts['curated']} genes", file=sys.stderr)
    print(f"    • background (all): {counts['total']} genes", file=sys.stderr)
    
    return {
        "set0": set0_file,
        "set1": set1_file,
        "set2": set2_file,
        "background": background_file
    }

def move_file(src, dst_dir, new_name=None):
    """Move a file to destination directory with optional new name."""
    src_path = Path(src)
    if not src_path.exists():
        sys.stderr.write(f"Warning: Source file not found: {src}\n")
        return False
    
    dst_path = Path(dst_dir)
    dst_path.mkdir(parents=True, exist_ok=True)
    
    if new_name:
        dst_file = dst_path / new_name
    else:
        dst_file = dst_path / src_path.name
    
    shutil.move(str(src_path), str(dst_file))
    return True

def copy_file(src, dst_dir, new_name=None):
    """Copy a file to destination directory with optional new name."""
    src_path = Path(src)
    if not src_path.exists():
        sys.stderr.write(f"Warning: Source file not found: {src}\n")
        return False

    dst_path = Path(dst_dir)
    dst_path.mkdir(parents=True, exist_ok=True)

    if new_name:
        dst_file = dst_path / new_name
    else:
        dst_file = dst_path / src_path.name

    shutil.copy2(str(src_path), str(dst_file))  # Changed from move to copy2
    return True

def step_initial_calculations(gene_files, output_dirs, libraries_enrichment, libraries_supercandidate, conda_env):
    """Step 1: Initial calculations"""
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 1: INITIAL CALCULATIONS", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    results_dir = output_dirs['main']
    figures_dir = output_dirs['docs_figures']
    enrichment_figures_dir = output_dirs['docs_enrichment_figures']
    enrichment_figures_dir.mkdir(parents=True, exist_ok=True)
    
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    script2 = get_absolute_path("scripts/2_supercandidate.py")
    script3 = get_absolute_path("scripts/3_count_supercandidate.py")
    script4 = get_absolute_path("scripts/4_compare_subsets.py")
    
    enrichment_libs_str = ",".join(libraries_enrichment)
    supercandidate_libs_str = ",".join(libraries_supercandidate)
    
    # Enrichment analysis
    run_command_conda([
        "python", str(script1),
        str(gene_files["set2"]),
        str(gene_files["set1"]),
        str(gene_files["set0"]),
        str(gene_files["background"]),
        "--output_dir", str(results_dir),
        "--libraries", enrichment_libs_str
    ], step_name="Running enrichment analysis", env_name=conda_env)
    
    # Move enrichment dot plots
    print("\nMoving enrichment dot plots...", file=sys.stderr)
    for png_file in results_dir.glob("*dotplot*.png"):
        if move_file(png_file, enrichment_figures_dir):
            print(f"  ✓ Moved: {png_file.name}", file=sys.stderr)
    
    # Supercandidate identification
    run_command_conda([
        "python", str(script2),
        str(gene_files["set1"]),
        str(gene_files["set2"]),
        "--output_dir", str(results_dir),
        "--libraries", supercandidate_libs_str
    ], step_name="Identifying supercandidate genes", env_name=conda_env)
    
    # Score distribution analysis
    gmt_file = gene_files.get("gmt_file")
    mondo_enrichment = results_dir / "gene_set2_MONDO_GROUPS_2026.tsv"
    run_command_conda([
        "python", str(script3),
        str(results_dir / "supercandidate.tsv"),
        str(gene_files["set1"]),
        str(mondo_enrichment),
        str(gmt_file),
        "-p", "1.00",
        "--plot", str(figures_dir / "dist_mondo_supercandidate.png"),
        "--output", str(results_dir / "dist_mondo_supercandidate.txt")
    ], step_name="Analyzing score distributions", env_name=conda_env)
    
    # Subset comparison
    run_command_conda([
        "python", str(script4),
        str(results_dir / "dist_mondo_supercandidate.txt"),
        "2", "3",
        "--label_col", "1",
        "--matrix", str(results_dir / "mondo_supercandidate_matrix.txt"),
        "--dendrogram", str(figures_dir / "mondo_supercandidate_fisher.png"),
        "--output_file", str(results_dir / "mondo_supercandidate_fisher.txt")
    ], step_name="Comparing subsets", env_name=conda_env)

def step_moe_analysis(gene_files, output_dirs, conda_env):
    """Step 2: MOE analysis"""
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 2: MOE ANALYSIS", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    results_dir = output_dirs['main']
    moe_dir = output_dirs['moe']
    enrichment_figures_dir = output_dirs['docs_enrichment_figures']
    moe_dir.mkdir(parents=True, exist_ok=True)
    enrichment_figures_dir.mkdir(parents=True, exist_ok=True)
    
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    supercandidate_file = results_dir / "supercandidate.tsv"
    
    score45_file = moe_dir / "score45.txt"
    with open(score45_file, 'w') as f:
        subprocess.run(
            ["awk", '-F', '\\t', 'NR>1 && !/^#/ && ($2==4 || $2==5) {print $1}', str(supercandidate_file)],
            stdout=f, check=True
        )
    
    score13_file = moe_dir / "score13.txt"
    with open(score13_file, 'w') as f:
        subprocess.run(
            ["awk", '-F', '\\t', 'NR>1 && !/^#/ && ($2==1 || $2==2 || $2==3) {print $1}', str(supercandidate_file)],
            stdout=f, check=True
        )
    
    # Run candidate enrichment
    run_command_conda([
        "python", str(script1),
        str(score45_file),
        str(score13_file),
        str(gene_files["set1"]),
        str(gene_files["background"]),
        "--output_dir", str(moe_dir),
        "--summary_file", "summary_file_candidate.tsv"
    ], step_name="Running candidate enrichment analysis", env_name=conda_env)
    
    # Move enrichment dot plots
    print("\nMoving enrichment dot plots...", file=sys.stderr)
    for png_file in moe_dir.glob("*dotplot*.png"):
        if move_file(png_file, enrichment_figures_dir):
            print(f"  ✓ Moved: {png_file.name}", file=sys.stderr)

def step_cancer_analysis(gene_files, output_dirs, conda_env):
    """Step 3: Cancer analysis (without scatter plots)"""
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 3: CANCER ANALYSIS", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    results_dir = output_dirs['main']
    cancer_dir = output_dirs['cancer_data']
    cancer_out = output_dirs['cancer_results']
    enrichment_figures_dir = output_dirs['docs_enrichment_figures']
    
    cancer_dir.mkdir(parents=True, exist_ok=True)
    cancer_out.mkdir(parents=True, exist_ok=True)
    enrichment_figures_dir.mkdir(parents=True, exist_ok=True)
    
    set_operations = get_absolute_path("utils/file_ops/set_operations.sh")
    run_cancer_comparison = get_absolute_path("utils/cancer_comparison/run_cancer_comparison.sh")
    script1 = get_absolute_path("scripts/1_enrichr_all.py")
    
    cancer_file = gene_files.get("cancer_file")
    
    # Extract cancer genes
    run_bash_command(
        f"{set_operations} {gene_files['set2']} {cancer_file} 1 1 -o intersect -w FALSE > {cancer_dir}/cancer_gs2.txt",
        step_name="Extracting cancer genes (set2)"
    )
    run_bash_command(
        f"{set_operations} {gene_files['set1']} {cancer_file} 1 1 -o intersect -w FALSE > {cancer_dir}/cancer_gs1.txt",
        step_name="Extracting cancer genes (set1)"
    )
    run_bash_command(
        f"{set_operations} {gene_files['set0']} {cancer_file} 1 1 -o intersect -w FALSE > {cancer_dir}/cancer_gs0.txt",
        step_name="Extracting cancer genes (set0)"
    )
    
    # Non-cancer genes
    run_bash_command(
        f"{set_operations} {gene_files['set1']} {cancer_dir}/cancer_gs1.txt 1 1 -o diff12 > {cancer_dir}/noncancer_gs1.txt",
        step_name="Extracting non-cancer genes (set1)"
    )
    run_bash_command(
        f"{set_operations} {gene_files['set2']} {cancer_dir}/cancer_gs2.txt 1 1 -o diff12 > {cancer_dir}/noncancer_gs2.txt",
        step_name="Extracting non-cancer genes (set2)"
    )
    run_bash_command(
        f"{set_operations} {gene_files['set0']} {cancer_dir}/cancer_gs0.txt 1 1 -o diff12 > {cancer_dir}/noncancer_gs0.txt",
        step_name="Extracting non-cancer genes (set0)"
    )
    
    # Cancer gene enrichment
    run_command_conda([
        "python", str(script1),
        str(cancer_dir / "cancer_gs2.txt"),
        str(cancer_dir / "noncancer_gs2.txt"),
        str(gene_files["set2"]),
        str(gene_files["background"]),
        "--output_dir", str(cancer_out),
        "--summary_file", "cancer_gs2_summary_file.tsv"
    ], step_name="Enrichment for cancer genes (set2)", env_name=conda_env)
    
    run_command_conda([
        "python", str(script1),
        str(cancer_dir / "cancer_gs1.txt"),
        str(cancer_dir / "noncancer_gs1.txt"),
        str(gene_files["set1"]),
        str(gene_files["background"]),
        "--output_dir", str(cancer_out),
        "--summary_file", "cancer_gs1_summary_file.tsv"
    ], step_name="Enrichment for cancer genes (set1)", env_name=conda_env)
    
    run_command_conda([
        "python", str(script1),
        str(cancer_dir / "cancer_gs0.txt"),
        str(cancer_dir / "noncancer_gs0.txt"),
        str(gene_files["set0"]),
        str(gene_files["background"]),
        "--output_dir", str(cancer_out),
        "--summary_file", "cancer_gs0_summary_file.tsv"
    ], step_name="Enrichment for cancer genes (set0)", env_name=conda_env)
    
    # Move enrichment dot plots
    print("\nMoving enrichment dot plots...", file=sys.stderr)
    for png_file in cancer_out.glob("*dotplot*.png"):
        if move_file(png_file, enrichment_figures_dir):
            print(f"  ✓ Moved: {png_file.name}", file=sys.stderr)
    
    # Aggregate cancer results and run Fisher test with BH correction
    cs1 = sum(1 for _ in open(cancer_dir / "cancer_gs1.txt"))
    ncs1 = sum(1 for _ in open(cancer_dir / "noncancer_gs1.txt"))
    s1 = sum(1 for _ in open(gene_files["set1"]))

    cs2 = sum(1 for _ in open(cancer_dir / "cancer_gs2.txt"))
    ncs2 = sum(1 for _ in open(cancer_dir / "noncancer_gs2.txt"))
    s2 = sum(1 for _ in open(gene_files["set2"]))

    # Create processed output directory
    processed_dir = cancer_out
    processed_dir.mkdir(parents=True, exist_ok=True)

    # Run cancer comparison for set 1
    cmd_set1 = (
        f"{run_cancer_comparison} "
        f"--cancer-file {cancer_out / 'cancer_gs1_GO_Biological_Process_2026.tsv'} "
        f"--noncancer-file {cancer_out / 'noncancer_gs1_GO_Biological_Process_2026.tsv'} "
        f"--set-file {results_dir / 'gene_set1_GO_Biological_Process_2026.tsv'} "
        f"--n-cancer {cs1} "
        f"--n-noncancer {ncs1} "
        f"--n-set {s1} "
        f"--output-file {processed_dir / 'cancer_compara_set1_adj.txt'}"
    )
    print(f"Running cancer comparison set1: {cmd_set1}", file=sys.stderr)
    run_bash_command(cmd_set1, step_name="Cancer comparison analysis (set1)")

    # Run cancer comparison for set 2
    cmd_set2 = (
        f"{run_cancer_comparison} "
        f"--cancer-file {cancer_out / 'cancer_gs2_GO_Biological_Process_2026.tsv'} "
        f"--noncancer-file {cancer_out / 'noncancer_gs2_GO_Biological_Process_2026.tsv'} "
        f"--set-file {results_dir / 'gene_set2_GO_Biological_Process_2026.tsv'} "
        f"--n-cancer {cs2} "
        f"--n-noncancer {ncs2} "
        f"--n-set {s2} "
        f"--output-file {processed_dir / 'cancer_compara_set2_adj.txt'}"
    )
    print(f"Running cancer comparison set2: {cmd_set2}", file=sys.stderr)
    run_bash_command(cmd_set2, step_name="Cancer comparison analysis (set2)")


def step_documentation(gene_files, output_dirs, conda_env):
    """Step 4: Generate documentation - scatter plots, tables, figures, and organize files"""
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 4: GENERATING DOCUMENTATION", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    docs_dir = output_dirs['docs']
    tables_dir = output_dirs['docs_tables']
    figures_dir = output_dirs['docs_figures']
    summaries_dir = output_dirs['docs_summaries']
    matrices_dir = output_dirs['docs_matrices']
    cancer_out = output_dirs['cancer_results']
    
    # Get paths to utility scripts
    script_scatter = get_absolute_path("utils/plotting/scatter_cancer.py")
    calculate_table_2 = get_absolute_path("utils/table_generation/calculate_table_2.sh")
    calculate_table_s1 = get_absolute_path("utils/table_generation/calculate_table_s1.sh")
    calculate_table_s2 = get_absolute_path("utils/table_generation/calculate_table_s2.sh")
    calculate_table_s3 = get_absolute_path("utils/table_generation/calculate_table_s3.sh")
    calculate_table_compara = get_absolute_path("utils/table_generation/calculate_table_compara.sh")
    venn_plot = get_absolute_path("utils/plotting/venn-plot.py")
    upset_plot = get_absolute_path("utils/plotting/upset-plot.py")
    violin_script = get_absolute_path("utils/plotting/violin.py")
    ks_test_script = get_absolute_path("utils/plotting/ks-test.py")
    set_operations = get_absolute_path("utils/file_ops/set_operations.sh")
    
    # Get file paths
    combined_file = gene_files.get("combined_file")
    cancer_file = gene_files.get("cancer_file")
    mondo_file = gene_files.get("mondo_file")
    supercandidate_file = output_dirs['main'] / "supercandidate.tsv"
    
    # Part A: Generate cancer scatter plots
    print("\n--- Part A: Generating cancer scatter plots ---", file=sys.stderr)
    
    if Path(script_scatter).exists():
        run_command_conda([
            "python", str(script_scatter),
            str(cancer_out / "cancer_compara_set1_adj.txt"), "9", "8", "16",
            "--output", str(figures_dir / "plot-compara-set1.png")
        ], step_name="Creating scatter plot (set1)", env_name=conda_env)
        
        run_command_conda([
            "python", str(script_scatter),
            str(cancer_out / "cancer_compara_set2_adj.txt"), "9", "8", "16",
            "--output", str(figures_dir / "plot-compara-set2.png")
        ], step_name="Creating scatter plot (set2)", env_name=conda_env)
        
        if Path(calculate_table_compara).exists():
            print("\nGenerating Table S7 ...", file=sys.stderr)
            table_s7_output = tables_dir / "table_s7.tsv"
            with open(table_s7_output, 'w') as f:
                subprocess.run(
                    ["bash", str(calculate_table_compara), str(cancer_out / "cancer_compara_set2_adj.txt")],
                    stdout=f, check=True
                )
            print(f"  ✓ Table S7 saved to: {table_s7_output}", file=sys.stderr)

            print("\nGenerating Table S8 ...", file=sys.stderr)
            table_s8_output = tables_dir / "table_s8.tsv"
            with open(table_s8_output, 'w') as f:
                subprocess.run(
                    ["bash", str(calculate_table_compara), str(cancer_out / "cancer_compara_set1_adj.txt")],
                    stdout=f, check=True
                )
            print(f"  ✓ Table S8 saved to: {table_s8_output}", file=sys.stderr)

    
    # Part B: Generate tables
    print("\n--- Part B: Generating tables ---", file=sys.stderr)
    
    if Path(calculate_table_2).exists() and combined_file:
        print("\nGenerating Table 2...", file=sys.stderr)
        table2_output = tables_dir / "table_2.tsv"
        with open(table2_output, 'w') as f:
            subprocess.run(
                ["bash", str(calculate_table_2), str(combined_file)],
                stdout=f, check=True
            )
        print(f"  ✓ Table 2 saved to: {table2_output}", file=sys.stderr)
    
    if Path(calculate_table_s1).exists():
        print("\nGenerating Table S1...", file=sys.stderr)
        table_s1_output = tables_dir / "table_s1.tsv"
        with open(table_s1_output, 'w') as f:
            subprocess.run(
                ["bash", str(calculate_table_s1), str(cancer_file), str(mondo_file)],
                stdout=f, check=True
            )
        print(f"  ✓ Table S1 saved to: {table_s1_output}", file=sys.stderr)
    
    if Path(calculate_table_s2).exists() and combined_file:
        print("\nGenerating Table S2...", file=sys.stderr)
        table_s2_output = tables_dir / "table_s2.tsv"
        with open(table_s2_output, 'w') as f:
            subprocess.run(
                ["bash", str(calculate_table_s2), str(combined_file),
                 str(cancer_file), str(mondo_file)],
                stdout=f, check=True
            )
        print(f"  ✓ Table S2 saved to: {table_s2_output}", file=sys.stderr)
    
    if Path(calculate_table_s3).exists() and supercandidate_file.exists():
        print("\nGenerating Table S3...", file=sys.stderr)
        table_s3_output = tables_dir / "table_s3.tsv"
        with open(table_s3_output, 'w') as f:
            subprocess.run(
                ["bash", str(calculate_table_s3), str(supercandidate_file),
                 str(cancer_file), str(mondo_file)],
                stdout=f, check=True
            )
        print(f"  ✓ Table S3 saved to: {table_s3_output}", file=sys.stderr)


    
    # Part C: Generate Venn diagrams
    print("\n--- Part C: Generating Venn diagrams ---", file=sys.stderr)
    
    if Path(venn_plot).exists() and combined_file:
        print("\nGenerating Venn diagram (candidate)...", file=sys.stderr)
        venn_candidate_output = figures_dir / "venn_candidate.png"
        run_command_conda([
            "python", str(venn_plot),
            str(combined_file), "candidate", "sanchis,candidate,orpha",
            str(venn_candidate_output)
        ], step_name="Generating Venn diagram (candidate)", env_name=conda_env)
        print(f"  ✓ Venn diagram saved to: {venn_candidate_output}", file=sys.stderr)
        
        print("\nGenerating Venn diagram (curated)...", file=sys.stderr)
        venn_curated_output = figures_dir / "venn_curated.png"
        run_command_conda([
            "python", str(venn_plot),
            str(combined_file), "curated", "sanchis,hc,orpha",
            str(venn_curated_output)
        ], step_name="Generating Venn diagram (curated)", env_name=conda_env)
        print(f"  ✓ Venn diagram saved to: {venn_curated_output}", file=sys.stderr)

    # Part C.2: Generate UpSet plots, kept alongside the Venn diagrams above. Unlike the
    # 3-set-capped Venn diagrams, these show every source-intersection combination
    # explicitly (SFARI included, plus the Orphanet neuro/develop AND-relationship that
    # defines curated eligibility, instead of pre-collapsing them into one "orpha" set).
    print("\n--- Part C.2: Generating UpSet plots ---", file=sys.stderr)

    if Path(upset_plot).exists() and combined_file:
        print("\nGenerating UpSet plot (candidate)...", file=sys.stderr)
        upset_candidate_output = figures_dir / "upset_candidate.png"
        run_command_conda([
            "python", str(upset_plot),
            str(combined_file), "candidate", "sanchis,candidate,neuro,develop,sfari",
            str(upset_candidate_output)
        ], step_name="Generating UpSet plot (candidate)", env_name=conda_env)
        print(f"  ✓ UpSet plot saved to: {upset_candidate_output}", file=sys.stderr)

        print("\nGenerating UpSet plot (curated)...", file=sys.stderr)
        upset_curated_output = figures_dir / "upset_curated.png"
        run_command_conda([
            "python", str(upset_plot),
            str(combined_file), "curated", "sanchis,hc,neuro,develop,sfari+",
            str(upset_curated_output)
        ], step_name="Generating UpSet plot (curated)", env_name=conda_env)
        print(f"  ✓ UpSet plot saved to: {upset_curated_output}", file=sys.stderr)

    # Part D: COPY summary files (CHANGED FROM MOVE)
    print("\n--- Part D: Copying summary files ---", file=sys.stderr)
    for dir_key in ['main', 'moe', 'cancer_results']:
        for summary_file in output_dirs[dir_key].glob("*summary*"):
            if copy_file(summary_file, summaries_dir):
                print(f"  ✓ Copied summary: {summary_file.name}", file=sys.stderr)
    
    # Part E: COPY matrix files (CHANGED FROM MOVE)
    print("\n--- Part E: Copying matrix files ---", file=sys.stderr)

    # Main matrices (from Step 1)
    for matrix_file in output_dirs['main'].glob("*aggregated_matrix*"):
        if copy_file(matrix_file, matrices_dir / "main"):
            print(f"  ✓ Copied: {matrix_file.name}", file=sys.stderr)

    # Cancer matrices (from Step 3)
    for matrix_file in output_dirs['cancer_results'].glob("*aggregated_matrix*"):
        if copy_file(matrix_file, matrices_dir / "cancer"):
            print(f"  ✓ Copied: {matrix_file.name}", file=sys.stderr)

    # Part F: COPY key files (CHANGED FROM MOVE)
    print("\n--- Part F: Copying key files ---", file=sys.stderr)
    
    key_files = [
        (output_dirs['main'] / "supercandidate.tsv",
         docs_dir,
         "supercandidate.tsv"),
        
        (output_dirs['main'] / "supercandidate_filtered.tsv",
         docs_dir,
         "supercandidate_filtered.tsv"),
        
        (output_dirs['main'] / "dist_mondo_supercandidate.txt",
         docs_dir,
         "dist_mondo_supercandidate.txt"),
        
        (output_dirs['cancer_results'] / "cancer_compara_set1_adj.txt",
         docs_dir,
         "cancer_compara_set1_adj.txt"),
        
        (output_dirs['cancer_results'] / "cancer_compara_set2_adj.txt",
         docs_dir,
         "cancer_compara_set2_adj.txt")
    ]
    
    for src_file, dst_dir, new_name in key_files:
        if src_file.exists():
            # File in original location
            if copy_file(src_file, dst_dir, new_name):
                print(f"  ✓ Copied: {src_file.name} → {dst_dir / new_name}", file=sys.stderr)
        elif (dst_dir / new_name).exists():
            # Already in docs/ from previous run
            print(f"  • Already in docs/: {new_name}", file=sys.stderr)
        else:
            sys.stderr.write(f"  ✗ Missing: {new_name}\n")
    
    # Part G: Copy files to tables (ALREADY USES copy_file)
    print("\n--- Part G: Copy key files to tables ---", file=sys.stderr)
    key_files = [
        (output_dirs['docs'] / "summaries/summary_table.tsv", 
         output_dirs['docs'] / "tables", 
         "table_3.tsv"),
        # Written by step_initial_calculations into results/main/, not into docs/.
        (output_dirs['main'] / "mondo_supercandidate_matrix.txt",
         matrices_dir,
         None),
        (matrices_dir / "mondo_supercandidate_matrix.txt",
         output_dirs['docs'] / "tables",
         "table_s4.tsv")
     
    ]
    
    for src_file, dst_dir, new_name in key_files:
        if copy_file(src_file, dst_dir, new_name):
            print(f"  ✓ Copied key file: {src_file.name}", file=sys.stderr)

    # Part H: Generate SysNDD comparison tables and figures
    print("\n--- Part H: Generating SysNDD comparison ---", file=sys.stderr)
    
    # Get SysNDD data file
    sysndd_file = gene_files.get("sysndd_file")
    
    # Check if SysNDD file exists
    if sysndd_file and Path(sysndd_file).exists() and Path(violin_script).exists() and Path(ks_test_script).exists():
        
        # Part H1: SysNDD distribution across gene classes
        print("\n--- Part H1: SysNDD distribution across gene classes ---", file=sys.stderr)
        
        violin_all_output = figures_dir / "violin_all.txt"
        
        # Run violin.py for all three gene classes
        # Column 5 = SysNDD score, Column 2 = gene class (1=no_evidence, 2=candidate, 3=curated)
        run_command_conda([
            "python", str(violin_script),
            str(sysndd_file), "5", "2",
            str(violin_all_output)
        ], step_name="Generating violin plot (all gene classes)", env_name=conda_env)
        
        # Move results to tables
        if violin_all_output.exists():
            if copy_file(violin_all_output, tables_dir, "table_s5a.tsv"):
                print(f"  ✓ Table S5a saved to: {tables_dir / 'table_s5a.tsv'}", file=sys.stderr)
        
        # Part H2: SysNDD distribution across MOE score tiers
        print("\n--- Part H2: SysNDD distribution across MOE score tiers ---", file=sys.stderr)
        
        # Combine SysNDD data with supercandidate.tsv (MOE scores)
        supercandidate_file = docs_dir / "supercandidate.tsv"
        sysndd_moe_file = docs_dir / "sysndd_moe.tsv"
        
        if supercandidate_file.exists():
            # Combine SysNDD scores with MOE scores
            # Extract columns: gene, MOE_score, SysNDD_score
            combined_cmd = (
                f"{set_operations} {sysndd_file} {supercandidate_file} 1 1 -o intersect "
                f"| awk -v OFS='\\t' '{{print $1,$2,$3,$7}}' > {sysndd_moe_file}"
            )
            run_bash_command(combined_cmd, step_name="Combining SysNDD scores with MOE scores")
            
            if sysndd_moe_file.exists():
                # Run violin.py for MOE subsets
                # Column 4 = SysNDD score, Column 2 = MOE score
                violin_moe_output = figures_dir / "violin_moe.txt"
                run_command_conda([
                    "python", str(violin_script),
                    str(sysndd_moe_file), "4", "2",
                    str(violin_moe_output)
                ], step_name="Generating violin plot (MOE subsets)", env_name=conda_env)
                
                # Move results to tables
                if violin_moe_output.exists():
                    if copy_file(violin_moe_output, tables_dir, "table_s5b.tsv"):
                        print(f"  ✓ Table S5b saved to: {tables_dir / 'table_s5b.tsv'}", file=sys.stderr)
                
                # Run ks-test.py for MOE subsets
                ks_moe_output = figures_dir / "ks_moe.txt"
                run_command_conda([
                    "python", str(ks_test_script),
                    str(sysndd_moe_file), "4", "2",
                    str(ks_moe_output)
                ], step_name="Generating KS test matrix (MOE subsets)", env_name=conda_env)
                
                # Move results to tables
                if ks_moe_output.exists():
                    if copy_file(ks_moe_output, tables_dir, "table_s6.tsv"):
                        print(f"  ✓ Table S6 saved to: {tables_dir / 'table_s6.tsv'}", file=sys.stderr)
                        
                # Move sysndd_moe.tsv to docs for reference
                #if copy_file(sysndd_moe_file, docs_dir):
                #    print(f"  ✓ sysndd_moe.tsv copied to docs/", file=sys.stderr)
    else:
        if not sysndd_file or not Path(sysndd_file).exists():
            sys.stderr.write("Warning: SysNDD file not found, skipping SysNDD comparison\n")
        if not Path(violin_script).exists():
            sys.stderr.write(f"Warning: violin.py not found: {violin_script}\n")
        if not Path(ks_test_script).exists():
            sys.stderr.write(f"Warning: ks-test.py not found: {ks_test_script}\n")

    # Create README
    readme_path = docs_dir / "README.md"
    with open(readme_path, 'w') as f:
        f.write("# NDD Analysis Results\n\n")
        f.write("## Directory Structure\n\n")
        f.write("### Figures\n- Venn diagrams\n- Score distribution plots\n- Cancer comparison scatter plots\n- SysNDD violin plots\n\n")
        f.write("### Figures/Enrichment\n- Enrichment dot plots for all libraries\n\n")
        f.write("### Tables\n")
        f.write("- Table 2: Gene classification summary\n")
        f.write("- Table 3: Summary enrichment analysis\n")
        f.write("- Table S1: MONDO and Cancer gene statistics\n")
        f.write("- Table S2: Combined gene set statistics\n")
        f.write("- Table S3: Supercandidate gene statistics\n")
        f.write("- Table S4: Supercandidate fisher table\n")
        f.write("- Table S5a: SysNDD score distribution across gene classes\n")
        f.write("- Table S5b: SysNDD score distribution across MOE score tiers\n")
        f.write("- Table S6: KS test results comparing SysNDD scores across MOE subsets\n")
        f.write("- Table S7: Curated set cancer compara summary\n")
        f.write("- Table S8: Candidate set cancer compara summary\n")
        f.write("- Table S9: MOE threshold operating characteristics against ClinVar P/LP\n")
        f.write("- Table S10: MOE vs ClinVar P/LP adjusted for gnomAD pLI (logistic regression)\n")
        f.write("- Table S11: ROC/PR AUC of MOE vs annotation-density baselines\n")
        f.write("- Table S12: MOE vs ClinVar P/LP adjusted for GO and PubMed counts\n")
        f.write("- Table S13: Spearman correlation of MOE with GO and PubMed annotation density\n")
        f.write("- Table S14: MONDO and HPO term coverage by MOE tier\n\n")
        f.write("### Validation (generated in Step 6)\n")
        f.write("- Independent validation of the MOE score against ClinVar P/LP status,\n")
        f.write("  gnomAD pLI constraint and annotation-density baselines\n\n")
        f.write("### Summaries\n- Summary tables for each analysis\n\n")
        f.write("### Matrices\n- Overlap matrices\n- Aggregated matrices\n\n")
        f.write("### Excel Files (generated in Step 5)\n")
        f.write("- `supfile_1_ndd_gene_lists.xlsx`: Complete gene list, MOE list, Cancer list, MONDO list, SysNDD MOE list\n")
        f.write("- `supfile_2_curated_enrichment.xlsx`: Curated set enrichment results\n")
        f.write("- `supfile_3_candidate_enrichment.xlsx`: Candidate set enrichment results\n")
        f.write("- `supfile_4_no_evidence_enrichment.xlsx`: No evidence set enrichment results\n")
        f.write("- `supfile_5_cancer_compara.xlsx`: Comparison of significance of cancer and non cancer gene functions \n")
        f.write("- `ndd-report.xlsx`: Complete report with all tables\n\n")
        f.write("## Key Files\n")
        f.write("- `supercandidate.tsv`: Supercandidate genes\n")
        f.write("- `cancer_compara_set1_adj.txt`: Cancer comparison results (set1)\n")
        f.write("- `cancer_compara_set2_adj.txt`: Cancer comparison results (set2)\n")
        f.write("- `sysndd_moe.tsv`: SysNDD scores combined with MOE scores\n")
    
    print(f"\n  ✓ Created README: {readme_path}", file=sys.stderr)

def step_generate_excel(gene_files, output_dirs, conda_env):
    """Step 5: Generate Excel files"""
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 5: GENERATING EXCEL FILES", file=sys.stderr)
    print("="*60, file=sys.stderr)
    
    tsv2excel_script = get_absolute_path("utils/table_generation/tsv2excel.py")
    
    if not Path(tsv2excel_script).exists():
        sys.stderr.write(f"Error: tsv2excel.py not found: {tsv2excel_script}\n")
        return
    
    docs_dir = output_dirs['docs']
    main_dir = output_dirs['main']
    moe_dir = output_dirs['moe']
    tables_dir = output_dirs['docs_tables']
    summaries_dir = output_dirs['docs_summaries']
    matrices_dir = output_dirs['docs_matrices']
    
    # ============================================
    # Excel 1: Gene Lists
    # ============================================
    print("\n--- Excel 1: Gene Lists ---", file=sys.stderr)
    
    excel1_files = []
    
    # Complete gene list
    combined_file = gene_files.get("combined_file")
    if combined_file and Path(combined_file).exists():
        excel1_files.append(f"{combined_file}:Complete_Gene_List:1")
    
    # MOE list (supercandidate.tsv)
    if docs_dir.exists():
        supercandidate_file = docs_dir / "supercandidate.tsv"
        if supercandidate_file.exists():
            excel1_files.append(f"{supercandidate_file}:Supercandidates:1")

    # Cancer list
    cancer_file = gene_files.get("cancer_file")
    if cancer_file and Path(cancer_file).exists():
        excel1_files.append(f"{cancer_file}:Cancer_List:0")
    
    # MONDO list
    mondo_file = gene_files.get("mondo_file")
    if mondo_file and Path(mondo_file).exists():
        excel1_files.append(f"{mondo_file}:MONDO_List:0")
    
    # SysNDD MOE file
    sysndd_moe_file = docs_dir / "sysndd_moe.tsv"
    if sysndd_moe_file.exists():
        excel1_files.append(f"{sysndd_moe_file}:SysNDD_MOE:1")
    
    if excel1_files:
        excel1_output = docs_dir / "supfile_1_ndd_gene_lists.xlsx"
        cmd = [get_conda_python(conda_env), str(tsv2excel_script), str(excel1_output)]
        cmd.extend(excel1_files)
        cmd.extend(['--separator','\\t'])
        
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, check=True)
            print(f"  ✓ Excel 1 saved to: {excel1_output}", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Warning: Could not generate Excel 1: {e}\n")
    
    # ============================================
    # Excel 2-4: Enrichment results for each gene set
    # ============================================
    
    enrichment_libraries = [
        ("GO_Biological_Process_2026", "GO_BP", "1"),
        ("GO_Cellular_Component_2026", "GO_CC", "1"),
        ("GO_Molecular_Function_2026", "GO_MF", "1"),
        # GO Slim: the redundancy-reduced counterpart to the three full-GO
        # sheets above (reviewer point R1.2 asks for a redundancy-reduction
        # step; GOslim collapses GO_BP's ~5,200 curated-set terms to 66).
        ("GOslim_Biological_Process_2026", "GOslim_BP", "1"),
        ("GOslim_Cellular_Component_2026", "GOslim_CC", "1"),
        ("GOslim_Molecular_Function_2026", "GOslim_MF", "1"),
        ("KEGG_2021_Human", "KEGG", "1"),
        ("Reactome_Pathways_2024", "Reactome", "1"),
        ("SynGO_2024", "SynGO", "1"),
        ("MONDO_GROUPS_2026", "MONDO_Groups", "1")
    ]
    
    gene_sets = [
        {"name": "curated", "prefix": "gene_set2", "excel_name": "supfile_2_curated_enrichment.xlsx"},
        {"name": "candidate", "prefix": "gene_set1", "excel_name": "supfile_3_candidate_enrichment.xlsx"},
        {"name": "no_evidence", "prefix": "gene_set0", "excel_name": "supfile_4_no_evidence_enrichment.xlsx"}
    ]
    
    for gene_set in gene_sets:
        print(f"\n--- Excel: {gene_set['name'].capitalize()} Enrichment ---", file=sys.stderr)
        
        excel_files = []
        
        for library, sheet_name , header in enrichment_libraries:
            file_path = main_dir / f"{gene_set['prefix']}_{library}.tsv"
            
            if file_path.exists():
                excel_files.append(f"{file_path}:{sheet_name}:{header}")
            else:
                # Silently omitting a missing sheet here is exactly the failure
                # mode that hides a stale results/ directory (step 1 not yet
                # rerun after a libraries_enrichment/pipeline change) -- make
                # it loud instead.
                sys.stderr.write(
                    f"  ⚠ Warning: {file_path} not found; sheet '{sheet_name}' "
                    f"will be MISSING from {gene_set['excel_name']}. Re-run step 1 "
                    f"(enrichment) if this library was recently added.\n"
                )
        
        if excel_files:
            excel_output = docs_dir / gene_set['excel_name']
            cmd = [get_conda_python(conda_env), str(tsv2excel_script), str(excel_output), 
                      "--skip_first_column"]
            cmd.extend(excel_files)
            cmd.extend(['--separator','\\t'])
            
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       text=True, check=True)
                print(f"  ✓ {gene_set['name'].capitalize()} enrichment Excel saved to: {excel_output}", file=sys.stderr)
                print(f"    Sheets included: {len(excel_files)}", file=sys.stderr)
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"Warning: Could not generate {gene_set['name']} Excel: {e}\n")
        else:
            sys.stderr.write(f"Warning: No enrichment files found for {gene_set['name']}\n")

    # ============================================
    # Excel 5: Cancer Compara - All tables
    # ============================================

    print("\n--- Excel 5: Cancer compara ---", file=sys.stderr)

    excel5_files = []   
 
    if docs_dir.exists():
        cancer_compara1 = docs_dir / "cancer_compara_set1_adj.txt"
        if cancer_compara1.exists():
            excel5_files.append(f"{cancer_compara1}:Cancer_Comparison_Set1:1")

        cancer_compara2 = docs_dir / "cancer_compara_set2_adj.txt"
        if cancer_compara2.exists():
            excel5_files.append(f"{cancer_compara2}:Cancer_Comparison_Set2:1")

    if excel5_files:
        excel5_output = docs_dir / "supfile_5_cancer_compara.xlsx"
        cmd = [get_conda_python(conda_env), str(tsv2excel_script), str(excel5_output)]
        cmd.extend(excel5_files)
        cmd.extend(['--separator','\\t'])

        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, check=True)
            print(f"  ✓ Cancer Compara Excel saved to: {excel5_output}", file=sys.stderr)
            print(f"    Total sheets included: {len(excel5_files)}", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Warning: Could not generate Cancer Compara Excel: {e}\n")
    else:
        sys.stderr.write("Warning: No files found to include in Cancer Compara\n")

    # ============================================
    # Excel 6: NDD Report - All tables
    # ============================================
    print("\n--- Excel 6: NDD Report ---", file=sys.stderr)
    
    report_files = []
    
    # Tables from docs/tables/
    if tables_dir.exists():
        for tsv_file in sorted(tables_dir.glob("*.tsv")):
            sheet_name = tsv_file.stem.replace("table_", "Table_").replace("table", "Table")
            report_files.append(f"{tsv_file}:{sheet_name}:1")
    
    if report_files:
        report_output = docs_dir / "ndd-report.xlsx"
        cmd = [get_conda_python(conda_env), str(tsv2excel_script), str(report_output)]
        cmd.extend(report_files)
        cmd.extend(['--separator','\\t'])
        
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, check=True)
            print(f"  ✓ NDD Report Excel saved to: {report_output}", file=sys.stderr)
            print(f"    Total sheets included: {len(report_files)}", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Warning: Could not generate NDD Report Excel: {e}\n")
    else:
        sys.stderr.write("Warning: No files found to include in NDD Report\n")

def step_reviewer_validation(gene_files, output_dirs, conda_env, data_dir):
    """Step 6: Independent validation of the MOE score.

    Addresses the reviewer requests that the MOE >= 4 threshold be validated
    against data the score was not derived from (R1.3), that the score be
    shown not to be a proxy for annotation density (R2.9), and that the
    behaviour of the individual score tiers be characterised (R2.10).

    Everything here runs off a single merged per-gene table built by
    utils/moe_validation/build_master_table.py, so the figures and tables quoted in the
    response to the reviewers are regenerated by the same command that
    regenerates the manuscript figures.
    """
    print("\n" + "="*60, file=sys.stderr)
    print("STEP 6: MOE VALIDATION (reviewer points R1.3, R2.9, R2.10)", file=sys.stderr)
    print("="*60, file=sys.stderr)

    main_dir = output_dirs['main']
    figures_dir = output_dirs['docs_figures']
    tables_dir = output_dirs['docs_tables']
    validation_dir = output_dirs['docs_validation']
    validation_dir.mkdir(parents=True, exist_ok=True)

    data_dir = Path(data_dir)
    supercandidate_file = main_dir / "supercandidate.tsv"
    if not supercandidate_file.exists():
        sys.stderr.write(
            f"Warning: {supercandidate_file} not found (run step 1 first); "
            f"skipping MOE validation\n"
        )
        return

    build_master_table = get_absolute_path("utils/moe_validation/build_master_table.py")
    clinvar_validation = get_absolute_path("utils/moe_validation/moe_clinvar_validation.py")
    pli_validation = get_absolute_path("utils/moe_validation/moe_pli_validation.py")
    bias_baseline = get_absolute_path("utils/moe_validation/moe_annotation_bias_baseline.py")
    violin_script = get_absolute_path("utils/plotting/violin.py")
    ks_test_script = get_absolute_path("utils/plotting/ks-test.py")
    bar_script = get_absolute_path("utils/plotting/bar.py")
    fisher_test_script = get_absolute_path("utils/plotting/fisher-test.py")
    roc_pr_script = get_absolute_path("utils/plotting/roc-pr.py")

    sfari_file = data_dir / "SFARI-Gene_genes_03-28-2024release_05-17-2024.csv"
    go_counts = data_dir / "gene_go_annotation_counts.tsv"
    pubmed_counts = data_dir / "gene_pubmed_counts.tsv"

    # Part A: merge the per-gene master table and the two validation inputs
    print("\n--- Part A: Building the per-gene master table ---", file=sys.stderr)
    run_command_conda([
        "python", str(build_master_table),
        "--data-dir", str(data_dir),
        "--supercandidate", str(supercandidate_file),
        "--sfari-file", str(sfari_file),
        "--output-dir", str(main_dir)
    ], step_name="Building master table", env_name=conda_env)

    master_table = main_dir / "ndd_master_table.tsv"
    pli_input = main_dir / "moe_pli_input.tsv"
    clinvar_input = main_dir / "moe_clinvar_input.tsv"

    for required in (master_table, pli_input, clinvar_input):
        if not required.exists():
            sys.stderr.write(f"Warning: {required} was not produced; skipping MOE validation\n")
            return

    # Part B: threshold operating characteristics against ClinVar P/LP (R1.3)
    print("\n--- Part B: ClinVar operating characteristics (R1.3) ---", file=sys.stderr)
    run_command_conda([
        "python", str(clinvar_validation),
        str(master_table),
        "--output-dir", str(validation_dir)
    ], step_name="MOE vs ClinVar P/LP validation", env_name=conda_env)

    # Part C: adjustment for gnomAD constraint (R1.3)
    print("\n--- Part C: gnomAD pLI adjusted model (R1.3) ---", file=sys.stderr)
    run_command_conda([
        "python", str(pli_validation),
        str(master_table),
        "--output-dir", str(validation_dir)
    ], step_name="MOE vs gnomAD pLI validation", env_name=conda_env)

    # Part D: annotation-density baselines (R2.9)
    print("\n--- Part D: Annotation-bias baselines (R2.9) ---", file=sys.stderr)
    if go_counts.exists() and pubmed_counts.exists():
        run_command_conda([
            "python", str(bias_baseline),
            str(master_table),
            "--go-counts", str(go_counts),
            "--pubmed-counts", str(pubmed_counts),
            "--output-dir", str(validation_dir)
        ], step_name="MOE vs annotation-density baselines", env_name=conda_env)
    else:
        sys.stderr.write(
            "Warning: GO/PubMed count files not found; skipping annotation-bias baselines\n"
        )

    # Part E: per-tier distributions and pairwise tests (R2.10)
    print("\n--- Part E: Per-tier distributions and pairwise tests (R2.10) ---", file=sys.stderr)
    # The MOE tiers are ordered, so the ClinVar comparison is run one-sided in
    # the direction the odds ratio is reported (higher tier enriched); the
    # continuous pLI comparison uses the two-sided KS test as before.
    tier_analyses = [
        (violin_script, pli_input, "moe_pli_violin", "pLI violin by MOE tier", []),
        (ks_test_script, pli_input, "moe_pli_ks", "pLI KS matrix by MOE tier", []),
        (bar_script, clinvar_input, "moe_clinvar_bar", "ClinVar proportion by MOE tier", []),
        (fisher_test_script, clinvar_input, "moe_clinvar_fisher",
         "ClinVar Fisher matrix by MOE tier", ["--alternative", "greater"]),
        (roc_pr_script, clinvar_input, "moe_clinvar_rocpr",
         "MOE ROC/PR against ClinVar P/LP", []),
    ]
    for script, input_file, stem, label, extra_args in tier_analyses:
        if not Path(script).exists():
            sys.stderr.write(f"Warning: {script} not found; skipping {label}\n")
            continue
        run_command_conda([
            "python", str(script),
            str(input_file), "2", "3",
            str(figures_dir / f"{stem}.txt"),
            "--png"
        ] + extra_args, step_name=f"Generating {label}", env_name=conda_env)

    # Part E.2: MONDO vs HPO term coverage across MOE tiers (R2.10)
    print("\n--- Part E.2: MONDO / HPO coverage by MOE tier (R2.10) ---", file=sys.stderr)
    mondo_hpo_script = get_absolute_path("utils/moe_validation/mondo_hpo_by_moe.py")
    if Path(mondo_hpo_script).exists():
        run_command_conda([
            "python", str(mondo_hpo_script),
            str(master_table),
            "--output-prefix", str(figures_dir / "mondo_hpo_by_moe"),
            "--png"
        ], step_name="Generating MONDO/HPO coverage by MOE tier", env_name=conda_env)
    else:
        sys.stderr.write(f"Warning: {mondo_hpo_script} not found; skipping MONDO/HPO coverage\n")

    # Part E.3: MOE component redundancy, HPO-augmented score, threshold choice,
    # ClinGen haploinsufficiency, annotation-density strata (R1.3, R2.9, R2.10)
    print("\n--- Part E.3: MOE component analyses (R1.3/R2.9/R2.10) ---", file=sys.stderr)
    component_script = get_absolute_path("utils/moe_validation/moe_component_analysis.py")
    if Path(component_script).exists():
        run_command_conda([
            "python", str(component_script), str(master_table),
            "--hpo", str(Path(data_dir) / "gene_hpo_terms.tsv"),
            "--libs", str(get_absolute_path("libs")),
            "--output-prefix", str(validation_dir / "moe_component")
        ], step_name="Running MOE component analyses", env_name=conda_env)
    else:
        sys.stderr.write(f"Warning: {component_script} not found; skipping component analyses\n")

    # Part E.4: Table 1 rule check and per-source unique contribution (R1.1, R2.3)
    print("\n--- Part E.4: Classification rule check and source contribution (R1.1/R2.3) ---", file=sys.stderr)
    source_script = get_absolute_path("utils/table_generation/source_contribution.py")
    if Path(source_script).exists():
        run_command_conda([
            "python", str(source_script), str(master_table),
            "--output-prefix", str(validation_dir / "source")
        ], step_name="Checking classification rule and source contribution", env_name=conda_env)
    else:
        sys.stderr.write(f"Warning: {source_script} not found; skipping source contribution\n")

    # Part F: promote the validation tables into docs/tables with table numbers
    print("\n--- Part F: Copying validation tables ---", file=sys.stderr)
    validation_tables = [
        (validation_dir / "moe_clinvar_operating_characteristics.tsv", "table_s9.tsv"),
        (validation_dir / "moe_pli_logistic_regression.tsv", "table_s10.tsv"),
        (validation_dir / "moe_annotation_bias_baseline_auc.tsv", "table_s11.tsv"),
        (validation_dir / "moe_annotation_bias_adjusted_model.tsv", "table_s12.tsv"),
        (validation_dir / "moe_annotation_density_correlations.tsv", "table_s13.tsv"),
        (figures_dir / "mondo_hpo_by_moe.tsv", "table_s14.tsv"),
        (validation_dir / "moe_component_redundancy.tsv", "table_s15.tsv"),
        (validation_dir / "moe_component_pca.tsv", "table_s16.tsv"),
        (validation_dir / "moe_component_score_variants.tsv", "table_s17.tsv"),
        (validation_dir / "moe_component_thresholds.tsv", "table_s18.tsv"),
        (validation_dir / "moe_component_clingen_hi.tsv", "table_s19.tsv"),
        (validation_dir / "moe_component_annotation_strata.tsv", "table_s20.tsv"),
        (validation_dir / "source_contribution.tsv", "table_s21.tsv"),
        (validation_dir / "source_sfari_only.tsv", "table_s22.tsv"),
        (validation_dir / "source_rule_check.tsv", "table_s23.tsv"),
        (validation_dir / "source_rule_counts.tsv", "table_s24.tsv"),
        (validation_dir / "source_decision_tree.tsv", "table_s25.tsv"),
    ]
    for src_file, new_name in validation_tables:
        if copy_file(src_file, tables_dir, new_name):
            print(f"  ✓ Copied validation table: {src_file.name} → {new_name}", file=sys.stderr)


def run_pipeline(data_dir=None, lib_dir=None, output_dir=None,
                 set0=None, set1=None, set2=None, background=None,
                 combined_file=None,
                 cancer_file=None, mondo_file=None, sysndd_file=None,
                 config_file=None, libraries_enrichment=None,
                 libraries_supercandidate=None, conda_env="ndd_analysis",
                 run_steps=None):
    """Main pipeline runner."""
    
    # Load YAML config if provided
    if config_file:
        print(f"Loading configuration from: {config_file}", file=sys.stderr)
        config = load_yaml_config(config_file)
        
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
        if not combined_file and 'combined_file' in config:
            combined_file = config['combined_file']
        if not cancer_file and 'cancer_file' in config:
            cancer_file = config['cancer_file']
        if not mondo_file and 'mondo_file' in config:
            mondo_file = config['mondo_file']
        if not sysndd_file and 'sysndd_file' in config:
            sysndd_file = config['sysndd_file']
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
    
    # Set default libraries
    if not libraries_enrichment:
        libraries_enrichment = [
            "GO_Biological_Process_2026",
            "GO_Biological_Process_Cancer_2026",
            "GO_Cellular_Component_2026",
            "GO_Molecular_Function_2026",
            "GOslim_Biological_Process_2026",
            "GOslim_Cellular_Component_2026",
            "GOslim_Molecular_Function_2026",
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
    
    # Set up paths
    base_dir = Path(__file__).parent
    
    if data_dir:
        data_dir = Path(data_dir).resolve()
    else:
        data_dir = get_absolute_path("data")
    
    if lib_dir:
        lib_dir = Path(lib_dir).resolve()
    else:
        lib_dir = get_absolute_path("libs")
    
    if output_dir:
        base_output_dir = Path(output_dir).resolve()
    else:
        base_output_dir = get_absolute_path("results")
    
    # Create output directory structure
    base_output_dir.mkdir(parents=True, exist_ok=True)
    
    output_dirs = {
        'main': base_output_dir / "main",
        'moe': base_output_dir / "moe",
        'cancer_data': base_output_dir / "cancer" / "data",
        'cancer_results': base_output_dir / "cancer" / "results",
        'docs': base_output_dir / "docs",
        'docs_figures': base_output_dir / "docs" / "figures",
        'docs_enrichment_figures': base_output_dir / "docs" / "figures" / "enrichment",
        'docs_tables': base_output_dir / "docs" / "tables",
        'docs_summaries': base_output_dir / "docs" / "summaries",
        'docs_matrices': base_output_dir / "docs" / "matrices",
        'docs_validation': base_output_dir / "docs" / "validation",
        'generated_sets': base_output_dir / "generated_sets"
    }
    
    for dir_path in output_dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    # If combined_file is provided, generate gene sets from it
    if combined_file:
        print(f"\n=== Using combined gene file: {combined_file} ===", file=sys.stderr)
        generated_files = generate_gene_sets_from_file(combined_file, output_dirs['generated_sets'])
        
        if not set0:
            set0 = generated_files["set0"]
        if not set1:
            set1 = generated_files["set1"]
        if not set2:
            set2 = generated_files["set2"]
        if not background:
            background = generated_files["background"]
    
    # Define file paths
    gene_files = {
        "set0": Path(set0).resolve() if set0 else data_dir / "gene_set0.txt",
        "set1": Path(set1).resolve() if set1 else data_dir / "gene_set1.txt",
        "set2": Path(set2).resolve() if set2 else data_dir / "gene_set2.txt",
        "background": Path(background).resolve() if background else data_dir / "gene_all.txt",
        "combined_file": Path(combined_file).resolve() if combined_file else None,
        "cancer_file": Path(cancer_file).resolve() if cancer_file else data_dir / "cancer.txt",
        "mondo_file": Path(mondo_file).resolve() if mondo_file else data_dir / "mondo.txt",
        "gmt_file": lib_dir / "MONDO_GROUPS_2026.gmt",
        "sysndd_file": Path(sysndd_file).resolve() if sysndd_file else data_dir / "SysNDD_all.txt"
    }
    
    # Check input files
    print("=== Checking input files ===", file=sys.stderr)
    for name, path in gene_files.items():
        if path is not None and Path(path).exists():
            print(f"  ✓ Found {name}: {path}", file=sys.stderr)
        elif name == "combined_file" and path is None:
            print(f"  • {name}: Not provided (using individual gene sets)", file=sys.stderr)
        elif name == "sysndd_file" and path is not None and not Path(path).exists():
            print(f"  • {name}: Not found at {path} (will skip SysNDD comparison)", file=sys.stderr)
        else:
            raise FileNotFoundError(f"Missing input file: {path}")
    
    # Print output directory structure
    print("\n=== Output Directory Structure ===", file=sys.stderr)
    for name, path in output_dirs.items():
        print(f"  • {name}: {path}", file=sys.stderr)
    
    # Determine which steps to run
    if run_steps is None:
        run_steps = [1, 2, 3, 4, 5, 6]
    
    print(f"\n=== Steps to run: {run_steps} ===", file=sys.stderr)
    
    # Run selected steps
    if 1 in run_steps:
        step_initial_calculations(gene_files, output_dirs, libraries_enrichment, libraries_supercandidate, conda_env)
    
    if 2 in run_steps:
        step_moe_analysis(gene_files, output_dirs, conda_env)
    
    if 3 in run_steps:
        step_cancer_analysis(gene_files, output_dirs, conda_env)
    
    if 4 in run_steps:
        step_documentation(gene_files, output_dirs, conda_env)
    
    if 5 in run_steps:
        step_generate_excel(gene_files, output_dirs, conda_env)
    
    if 6 in run_steps:
        step_reviewer_validation(gene_files, output_dirs, conda_env, data_dir)
    
    print("\n=== Pipeline completed successfully ===", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="NDD Analysis Pipeline - Complete analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all steps
  %(prog)s -c config.yml
  
  # Run with combined file
  %(prog)s -c config.yml --combined_file /path/to/combined_genes.tsv
  
  # Run only documentation generation
  %(prog)s -c config.yml --steps 4
  
  # Run only Excel generation
  %(prog)s -c config.yml --steps 5
  
  # Run documentation and Excel
  %(prog)s -c config.yml --steps 4 5
  
  # Run only the MOE validation analyses (reviewer points R1.3, R2.9, R2.10)
  %(prog)s -c config.yml --steps 6
        """
    )
    
    parser.add_argument("--steps", "-s", type=int, nargs='+', choices=[1, 2, 3, 4, 5, 6], default=None)
    parser.add_argument("-d", "--data_dir", type=str, default=None)
    parser.add_argument("-l", "--lib_dir", type=str, default=None)
    parser.add_argument("-o", "--output_dir", type=str, default=None)
    parser.add_argument("-s0", "--set0", type=str, default=None)
    parser.add_argument("-s1", "--set1", type=str, default=None)
    parser.add_argument("-s2", "--set2", type=str, default=None)
    parser.add_argument("-b", "--background", type=str, default=None)
    parser.add_argument("--combined_file", type=str, default=None)
    parser.add_argument("--cancer_file", type=str, default=None)
    parser.add_argument("--mondo_file", type=str, default=None)
    parser.add_argument("--sysndd_file", type=str, default=None)
    parser.add_argument("--libraries_enrichment", type=str, default=None)
    parser.add_argument("--libraries_supercandidate", type=str, default=None)
    parser.add_argument("--conda_env", type=str, default="ndd_analysis")
    parser.add_argument("-c", "--config", type=str, default=None)
    
    args = parser.parse_args()
    
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
        combined_file=args.combined_file,
        cancer_file=args.cancer_file,
        mondo_file=args.mondo_file,
        sysndd_file=args.sysndd_file,
        config_file=args.config,
        libraries_enrichment=libraries_enrichment,
        libraries_supercandidate=libraries_supercandidate,
        conda_env=args.conda_env,
        run_steps=args.steps
    )
