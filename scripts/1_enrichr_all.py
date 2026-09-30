#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path
import gseapy as gp
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns  # For heatmap visualization

# utils/enrichment/enrichr_odds_ratio_ci.py appends a 95% CI to the Odds Ratio
# column (reviewer point R1.2: "Report effect sizes, gene counts and
# confidence intervals"). It lives under utils/, a sibling of scripts/, so add
# it to sys.path relative to this file rather than assuming a fixed layout.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils" / "enrichment"))
from enrichr_odds_ratio_ci import add_ci  # noqa: E402


def read_gene_list(file_path):
    """
    Reads a gene list from a file and returns it as a list.
    Handles empty files by returning an empty list.
    Skips comment lines starting with #.
    """
    genes = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                # Keep the symbol as written (HGNC case, e.g. C9orf72). Case is
                # reconciled per library by resolve_symbols(): libraries do not
                # agree (GO/KEGG/Reactome/SynGO write C9ORF72, GOslim/MONDO write
                # C9orf72), so neither blanket upper- nor lower-casing is safe.
                genes.append(line)
    except FileNotFoundError:
        sys.stderr.write(f"Warning: File not found: {file_path}\n")
        return []
    return genes


def library_symbols(library):
    """Return (symbols, index) for a library: the set of gene symbols exactly as
    written in libs/<library>.gmt, and a map UPPER(symbol) -> {symbols}."""
    library_file = os.path.join(f"{ndd_path}/libs", f"{library}.gmt")
    symbols = set()
    with open(library_file) as fh:
        for line in fh:
            symbols.update(g for g in line.rstrip("\n").split("\t")[2:] if g)
    index = {}
    for g in symbols:
        index.setdefault(g.upper(), set()).add(g)
    return symbols, index


def resolve_symbols(genes, symbols, index):
    """Case exception for matching gene symbols to one library.

    gseapy matches symbols case-sensitively. A gene is used as written when the
    library contains it exactly; otherwise, if the library contains exactly one
    symbol that equals it ignoring case (e.g. C9orf72 vs C9ORF72, typically the
    open-reading-frame genes), the library's spelling is used and the case is
    reported. Genes absent from the library, or whose case-insensitive match is
    ambiguous (several spellings in the library), are upper-cased: they match
    no term either way, but gseapy (enrichr.py, _local_enrichment) upper-cases
    the whole query and background when <90% of the query is upper case and
    the first library terms are upper case -- which would undo the resolution
    above (e.g. C18orf32 -> C18ORF32, absent from MONDO). Upper-casing the
    unmatched genes keeps that heuristic from firing (checked by
    assert_no_gseapy_recase()).
    Returns (resolved_genes, [(input, library_symbol_or_candidates, status)]).
    """
    resolved, exceptions = [], []
    for g in genes:
        if g in symbols:
            resolved.append(g)
            continue
        hits = index.get(g.upper(), set())
        if len(hits) == 1:
            lib_g = next(iter(hits))
            resolved.append(lib_g)
            exceptions.append((g, lib_g, "case_match"))
        else:
            resolved.append(g.upper())
            if len(hits) > 1:
                exceptions.append((g, ",".join(sorted(hits)), "ambiguous_case"))
    return resolved, exceptions


def _mostly_upper(genes):
    """gseapy's check_uppercase(): >=90% of symbols are upper case."""
    genes = [str(g) for g in genes]
    return bool(genes) and sum(g.isupper() for g in genes) / len(genes) >= 0.9


def assert_no_gseapy_recase(query, library):
    """Fail if gseapy would silently upper-case this query (see resolve_symbols)."""
    gmt = gp.get_library(os.path.join(f"{ndd_path}/libs", f"{library}.gmt"))
    top = list(gmt.keys())[:min(len(gmt), 10)]
    if all(_mostly_upper(gmt[k]) for k in top) and not _mostly_upper(query):
        sys.exit(f"Error: {library}: gseapy would upper-case the query and undo the "
                 f"symbol case resolution; aborting.")


def perform_enrichment(gene_lists, background, library, output_dir, input_files):
    """
    Perform enrichment analysis for multiple gene lists using gseapy for a single library.
    Skips empty gene lists.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Perform enrichment analysis for each gene list
    for i, (gene_list, input_file) in enumerate(zip(gene_lists, input_files)):
        # Get the basename of the input file (without extension)
        basename = os.path.splitext(os.path.basename(input_file))[0]
        
        # Skip empty gene lists
        if len(gene_list) == 0:
            sys.stderr.write(f"Warning: Gene list {basename} is empty. Skipping enrichment for {library}.\n")
            # Create an empty output file to maintain consistency
            output_file = os.path.join(output_dir, f"{basename}_{library}.tsv")
            empty_df = pd.DataFrame(columns=["Gene_set", "Term", "Overlap", "P-value",
                                             "Adjusted P-value", "Odds Ratio",
                                             "OR_CI_low", "OR_CI_high", "Combined Score"])
            empty_df.to_csv(output_file, sep="\t", index=False)
            continue
        
        sys.stderr.write(f"Processing gene list {basename} for library {library}...\n")

        # Run enrichment analysis
        library_file = os.path.join(f"{ndd_path}/libs", f"{library}.gmt")
        
        try:
            enr_results = gp.enrichr(
                gene_list=gene_list,
                gene_sets=library_file,
                background=background,
                outdir=None,  # Do not save output automatically
            )
        except Exception as e:
            sys.stderr.write(f"Error during enrichment for {basename} with {library}: {e}\n")
            # Create an empty output file
            output_file = os.path.join(output_dir, f"{basename}_{library}.tsv")
            empty_df = pd.DataFrame(columns=["Gene_set", "Term", "Overlap", "P-value",
                                             "Adjusted P-value", "Odds Ratio",
                                             "OR_CI_low", "OR_CI_high", "Combined Score"])
            empty_df.to_csv(output_file, sep="\t", index=False)
            continue

        # Save results to a file
        output_file = os.path.join(output_dir, f"{basename}_{library}.tsv")
        if enr_results.res2d is None or len(enr_results.res2d) == 0:
            sys.stderr.write(f"No matching annotation for gene list {basename}.\n")
            # Create an empty output file
            empty_df = pd.DataFrame(columns=["Gene_set", "Term", "Overlap", "P-value",
                                             "Adjusted P-value", "Odds Ratio",
                                             "OR_CI_low", "OR_CI_high", "Combined Score"])
            empty_df.to_csv(output_file, sep="\t", index=False)
            continue
        
        enr_results.res2d = enr_results.res2d.sort_values(by=["Adjusted P-value", "Odds Ratio"], ascending=[True, False])
        # list_size/background_size are exactly the sizes gseapy used for this
        # call, so the 2x2 table reconstructed from Overlap is exact, not an
        # approximation from a size passed in separately.
        enr_results.res2d = add_ci(
            enr_results.res2d, list_size=len(gene_list), background_size=len(background),
            on_bad_row="nan",
        )
        enr_results.res2d.to_csv(output_file, sep="\t", index=False)
        sys.stderr.write(f"  Results saved to {output_file}\n")

        # Generate a dot plot for the top 20 enriched terms
        generate_dotplot(enr_results.res2d, basename, library, output_dir)


def generate_dotplot(enrichment_results, basename, library, output_dir):
    """
    Generates a dot plot for the top 20 enriched terms with a p-value threshold of 0.01.
    """
    # Filter results based on p-value threshold
    min_non_zero_pvalue = 1e-324
    filtered_results = enrichment_results[enrichment_results["Adjusted P-value"] < 0.01]
    
    if len(filtered_results) == 0:
        sys.stderr.write(f"  No significant terms found for {basename} ({library}) to generate a dot plot.\n")
        return
    
    filtered_results = filtered_results.copy()
    filtered_results.loc[filtered_results['Adjusted P-value'] == 0.0, 'Adjusted P-value'] = min_non_zero_pvalue
    filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Combined Score'] = 1e4

    # Select the top 20 terms
    filtered_results_sorted = filtered_results.sort_values(by=["Adjusted P-value", "Odds Ratio"], ascending=[True, False])
    top_terms = filtered_results_sorted.head(20)

    # Save the figure
    output_file = os.path.join(output_dir, f"{basename}_{library}_dotplot.png")
    try:
        gp.dotplot(
            top_terms,
            cutoff=0.01,
            top_term=20,
            size=5,
            ofname=output_file,
        )
        sys.stderr.write(f"  Dot plot saved to {output_file}\n")
    except Exception as e:
        sys.stderr.write(f"  Warning: Could not generate dot plot: {e}\n")


def plot_overlap_heatmap(overlap_matrix, gene_set_names, output_dir, library, aggregated_matrix):
    """
    Plots the overlap matrix as a heatmap with color and saves it as a PNG file.
    """
    # Create a DataFrame for the overlap matrix
    overlap_df = pd.DataFrame(overlap_matrix, index=gene_set_names, columns=gene_set_names)

    # Convert the aggregated matrix to a DataFrame
    annotation_matrix = pd.DataFrame(aggregated_matrix[1:, 1:], index=gene_set_names, columns=gene_set_names)

    # Plot the heatmap
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        overlap_df,
        annot=annotation_matrix,
        fmt="",
        cmap='OrRd',
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="black",
    )

    # Customize the colorbar
    if ax.collections[0].colorbar is not None:
        cbar = ax.collections[0].colorbar
        cbar.outline.set_edgecolor('black')
        cbar.outline.set_linewidth(0.5)

    plt.title(f"Overlap Matrix for {library}", fontsize=14, pad=20)
    ax.tick_params(axis='x', which='both', pad=5)
    ax.tick_params(axis='y', which='both', pad=5)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1)

    output_file = os.path.join(output_dir, f"overlap_heatmap_{library}.png")
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    sys.stderr.write(f"  Overlap heatmap saved to {output_file}\n")


def aggregate_matrices(p_value_matrix, overlap_matrix, output_file, gene_set_names):
    """
    Aggregates the p-value and overlap matrices into a single matrix.
    """
    assert p_value_matrix.shape == overlap_matrix.shape, "Matrices must be of the same shape"
    
    n = p_value_matrix.shape[0]
    aggregated_matrix = np.empty((n + 1, n + 1), dtype=object)
    
    aggregated_matrix[0, 0] = "set"
    for i in range(n):
        aggregated_matrix[0, i + 1] = gene_set_names[i]
        aggregated_matrix[i + 1, 0] = gene_set_names[i]
    
    for i in range(n):
        for j in range(n):
            if i == j:
                aggregated_matrix[i + 1, j + 1] = '-'
            elif i < j:
                aggregated_matrix[i + 1, j + 1] = f"{p_value_matrix[i, j]:.1e}"
            else:
                aggregated_matrix[i + 1, j + 1] = f"{overlap_matrix[i, j]:.3f}"
    
    np.savetxt(output_file, aggregated_matrix, delimiter='\t', fmt='%s')
    return aggregated_matrix


def calculate_overlap_matrix(input_files, library, output_dir, adjusted_p_threshold, background_terms):
    """
    Calculates the overlap matrix for multiple gene sets (dynamic size).
    """
    n_sets = len(input_files)

    if n_sets < 2:
        sys.stderr.write(f"Skipping overlap matrix for {library}: only 1 gene set provided.\n")
        return
 
    # Initialize matrices
    overlap_fraction_matrix = np.zeros((n_sets, n_sets), dtype=float)
    p_value_matrix = np.zeros((n_sets, n_sets), dtype=float)

    # Read significant terms
    significant_terms = {}
    for i, input_file in enumerate(input_files):
        basename = os.path.splitext(os.path.basename(input_file))[0]
        file_name = f"{basename}_{library}.tsv"
        file_path = os.path.join(output_dir, file_name)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep="\t")
            df_filtered = df[df["Adjusted P-value"] < adjusted_p_threshold]
            significant_terms[f"gene_list_{i + 1}"] = set(df_filtered["Term"])
        else:
            sys.stderr.write(f"Warning: File {file_path} not found.\n")
            significant_terms[f"gene_list_{i + 1}"] = set()

    # Calculate overlap and p-values
    for i in range(n_sets):
        for j in range(n_sets):
            set1 = significant_terms[f"gene_list_{i + 1}"]
            set2 = significant_terms[f"gene_list_{j + 1}"]
            overlap = len(set1.intersection(set2))
            min_size = min(len(set1), len(set2))
            overlap_fraction = overlap / min_size if min_size > 0 else 0.0
            overlap_fraction_matrix[i, j] = overlap_fraction

            M = len(background_terms)
            n = len(set1)
            N = len(set2)
            k = overlap
            p_value = hypergeom.sf(k - 1, M, n, N)
            p_value_matrix[i, j] = p_value

    gene_set_names = [os.path.splitext(os.path.basename(file))[0] for file in input_files]
    set_names_joined = "_vs_".join(os.path.splitext(os.path.basename(file))[0] for file in input_files)
    aggregated_output_file = os.path.join(output_dir, f"{set_names_joined}_aggregated_matrix_{library}.tsv")
    aggregated_matrix = aggregate_matrices(p_value_matrix, overlap_fraction_matrix, aggregated_output_file, gene_set_names)
    sys.stderr.write(f"Aggregated matrix saved to {aggregated_output_file}\n")

    plot_overlap_heatmap(overlap_fraction_matrix, gene_set_names, output_dir, library, aggregated_matrix)


def get_background_terms(library, background_genes):
    """
    Returns the set of terms associated with the background genes.
    """
    library_file = os.path.join(f"{ndd_path}/libs", f"{library}.gmt")
    gene_sets = gp.get_library(library_file)
    
    background_genes_set = set(background_genes)
    
    background_terms = set()
    for term, genes in gene_sets.items():
        if any(gene in background_genes_set for gene in genes):
            background_terms.add(term)

    return background_terms


def generate_summary_table(input_files, libraries, output_dir, adjusted_p_threshold):
    """
    Generates a summary table for any number of gene sets (including 1).
    """
    n_sets = len(input_files)
    set_names = [os.path.splitext(os.path.basename(file))[0] for file in input_files]
    
    summary_data = []

    for library in libraries:
        significant_terms = {}
        total_terms = {}
        all_terms = set()

        for i, input_file in enumerate(input_files):
            basename = os.path.splitext(os.path.basename(input_file))[0]
            file_name = f"{basename}_{library}.tsv"
            file_path = os.path.join(output_dir, file_name)
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep="\t")
                df_filtered = df[df["Adjusted P-value"] < adjusted_p_threshold]
                significant_terms[f"gene_list_{i + 1}"] = set(df_filtered["Term"])
                total_terms[f"gene_list_{i + 1}"] = len(df["Term"])
                all_terms.update(df["Term"])
            else:
                sys.stderr.write(f"Warning: File {file_path} not found.\n")
                significant_terms[f"gene_list_{i + 1}"] = set()
                total_terms[f"gene_list_{i + 1}"] = 0

        # Prepare row data
        row_data = [library]
        
        # Add significant terms counts
        for i in range(n_sets):
            row_data.append(f"{len(significant_terms[f'gene_list_{i + 1}'])} ({total_terms[f'gene_list_{i + 1}']})")
        
        # Add pairwise overlaps (ONLY IF n_sets >= 2)
        if n_sets >= 2:
            M = len(all_terms)
            for i in range(n_sets):
                for j in range(i + 1, n_sets):
                    set_i = significant_terms[f"gene_list_{i + 1}"]
                    set_j = significant_terms[f"gene_list_{j + 1}"]
                    intersection = len(set_i.intersection(set_j))
                    min_size = min(len(set_i), len(set_j))
                    fraction = intersection / min_size if min_size > 0 else 0
                    
                    p_value = hypergeom.sf(intersection - 1, M, len(set_i), len(set_j))
                    row_data.append(f"{fraction:.3f} ({p_value:.2e})")
        
        summary_data.append(row_data)

    # Build column names dynamically
    columns = ["Library"]
    for i in range(n_sets):
        columns.append(f"Significant Terms ({set_names[i]})")
    
    # Only add overlap columns if n_sets >= 2
    if n_sets >= 2:
        for i in range(n_sets):
            for j in range(i + 1, n_sets):
                columns.append(f"Fraction Overlap ({set_names[i]} & {set_names[j]})")

    summary_df = pd.DataFrame(summary_data, columns=columns)
    return summary_df


def main():
    # Get the directory of the script
    script_dir = Path(__file__).parent.parent
    global ndd_path
    ndd_path = script_dir

    # Set up argument parser with NARG for positional arguments
    parser = argparse.ArgumentParser(description='Perform enrichment analysis on gene lists')
    parser.add_argument('gene_lists', nargs='+', 
                       help='Gene list files (at least 2 required, last one is background)')
    parser.add_argument('--output_dir', default=f"{ndd_path}/results",
                       help='Output directory for results')
    parser.add_argument('--summary_file', default="summary_table.tsv",
                       help='Name of the summary output file')
    parser.add_argument('--libraries', type=str, default=None,
                       help='Comma-separated list of libraries to use (default: all libraries)')
    args = parser.parse_args()

    # Check minimum number of positional arguments
    if len(args.gene_lists) < 2:
        sys.stderr.write("Error: At least 2 positional arguments required: 1 gene sets + 1 background\n")
        sys.stderr.write("Usage: script.py gene_set1 gene_set2 [gene_set3 ...] background\n")
        sys.exit(1)

    # Separate gene lists and background
    input_files = args.gene_lists[:-1]
    background_file = args.gene_lists[-1]

    # Read gene lists
    gene_lists = []
    for file in input_files:
        gene_list = read_gene_list(file)
        gene_lists.append(gene_list)
        sys.stderr.write(f"Number of genes in {os.path.basename(file)}: {len(gene_list)}\n")

    background_genes = read_gene_list(background_file)
    sys.stderr.write(f"Number of genes in background {os.path.basename(background_file)}: {len(background_genes)}\n")

    # Define libraries
    if args.libraries:
        libraries = [lib.strip() for lib in args.libraries.split(",")]
    else:
        libraries = [
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

    adjusted_p_threshold = 0.01

    # Process each library
    os.makedirs(args.output_dir, exist_ok=True)
    case_log = []   # rows: library, gene_list, input_symbol, library_symbol, status
    for library in libraries:
        # Reconcile symbol case with this library (see resolve_symbols) for the
        # gene lists AND the background: gseapy restricts every term to the
        # background, so a case mismatch there also drops genes from terms.
        symbols, index = library_symbols(library)
        lib_gene_lists = []
        for gl, f in zip(gene_lists, input_files):
            resolved, exc = resolve_symbols(gl, symbols, index)
            lib_gene_lists.append(resolved)
            case_log += [(library, os.path.basename(f), *e) for e in exc]
        lib_background, exc = resolve_symbols(background_genes, symbols, index)
        case_log += [(library, os.path.basename(background_file), *e) for e in exc]
        n_case = sum(1 for e in exc if e[2] == "case_match")
        n_amb = sum(1 for e in exc if e[2] == "ambiguous_case")
        if n_case:
            sys.stderr.write(f"Warning: {library}: {n_case} background symbols matched only ignoring case "
                             f"(library spelling used; see symbol_case_matches.tsv)\n")
        if n_amb:
            sys.stderr.write(f"Warning: {library}: {n_amb} symbols have ambiguous case matches and were left unmatched\n")

        background_terms = get_background_terms(library, lib_background)
        sys.stderr.write(f"Number of terms associated with background genes in {library}: {len(background_terms)}\n")

        for q in lib_gene_lists:
            if q:
                assert_no_gseapy_recase(q, library)
        perform_enrichment(lib_gene_lists, lib_background, library, args.output_dir, input_files)
        
        # Calculate overlap matrix if we have at least 2 gene sets
        if len(input_files) >= 2:
            calculate_overlap_matrix(input_files, library, args.output_dir, adjusted_p_threshold, background_terms)

    case_log_path = os.path.join(args.output_dir, "symbol_case_matches.tsv")
    pd.DataFrame(case_log, columns=["library", "gene_list", "input_symbol", "library_symbol", "status"]
                 ).to_csv(case_log_path, sep="\t", index=False)
    sys.stderr.write(f"Symbol case exceptions ({len(case_log)} rows) saved to {case_log_path}\n")

    # Generate summary table
    summary_table = generate_summary_table(input_files, libraries, args.output_dir, adjusted_p_threshold)
    summary_output_path = os.path.join(args.output_dir, args.summary_file)
    summary_table.to_csv(summary_output_path, sep="\t", index=False)
    sys.stderr.write(f"Summary table saved to {summary_output_path}\n")


if __name__ == "__main__":
    main()
