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


def read_gene_list(file_path):
    """
    Reads a gene list from a file and returns it as a list.
    """
    with open(file_path, "r") as file:
        return [line.strip().upper() for line in file if line.strip()]


def perform_enrichment(gene_lists, background, library, output_dir, input_files):
    """
    Perform enrichment analysis for multiple gene lists using gseapy for a single library.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Perform enrichment analysis for each gene list
    for i, (gene_list, input_file) in enumerate(zip(gene_lists, input_files)):
        # Get the basename of the input file (without extension)
        basename = os.path.splitext(os.path.basename(input_file))[0]
        sys.stderr.write(f"Processing gene list {basename} for library {library}...\n")

        # Run enrichment analysis
        library_file=os.path.join(f"{ndd_path}/libs", f"{library}.gmt")
        enr_results = gp.enrichr(
            gene_list=gene_list,
	    #gene_sets=library,
            gene_sets=library_file,
            #organism='human',
            background=background,
            outdir=None,  # Do not save output automatically
        )

        # Save results to a file
        output_file = os.path.join(output_dir, f"{basename}_{library}.tsv")
        if enr_results.res2d is None:
            sys.stderr.write(f"No maching annotation for gene list {basename}.\n")
            return
        enr_results.res2d = enr_results.res2d.sort_values(by=["Adjusted P-value", "Odds Ratio"], ascending=[True, False])
        enr_results.res2d.to_csv(output_file, sep="\t", index=False)
        sys.stderr.write(f"  Results saved to {output_file}\n")

        # Generate a dot plot for the top 20 enriched terms
        generate_dotplot(enr_results.res2d, basename, library, output_dir)


def generate_dotplot(enrichment_results, basename, library, output_dir):
    """
    Generates a dot plot for the top 20 enriched terms with a p-value threshold of 0.01.
    The figure is saved with the basename of the gene set and the name of the library.
    If no terms meet the threshold, a message is printed, and no plot is generated.
    """
    # Filter results based on p-value threshold
    min_non_zero_pvalue = 1e-324
    filtered_results = enrichment_results[enrichment_results["Adjusted P-value"] < 0.01]
    filtered_results.loc[filtered_results['Adjusted P-value'] == 0.0, 'Adjusted P-value'] = min_non_zero_pvalue
    odds_ratio_std=filtered_results['Odds Ratio'].std()
    filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Combined Score'] = 1e4   
    #filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Combined Score'] = (-np.log10(min_non_zero_pvalue) * (np.log10(filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Odds Ratio']) / odds_ratio_std))
    #filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Combined Score'] = (324.0* (np.log10(filtered_results.loc[filtered_results['Combined Score'] == np.inf, 'Odds Ratio']) / odds_ratio_std))



    # Check if there are any significant terms
    if filtered_results.empty:
        sys.stderr.write(f"  No significant terms found for {basename} ({library}) to generate a dot plot.\n")
        return

    # Select the top 20 terms (or fewer if there are not enough)
    filtered_results_sorted = filtered_results.sort_values(by=["Adjusted P-value", "Odds Ratio"], ascending=[True, False])
    top_terms = filtered_results_sorted.head(20)

    #print(top_terms.to_string())

    # Save the figure in the output directory
    output_file = os.path.join(output_dir, f"{basename}_{library}_dotplot.png")
    gp.dotplot(
        top_terms,
        #title=f"Top {len(top_terms)} Enriched Terms for {basename} ({library})",  # Add a title
        cutoff=0.01,  # P-value cutoff
        top_term=20,  # Show only top 20 terms
        size=5,  # Size of the dots
        ofname=output_file,  # Save the plot directly to this file
    )
    sys.stderr.write(f"  Dot plot saved to {output_file}\n")


def plot_overlap_heatmap(overlap_matrix, gene_set_names, output_dir, library, aggregated_matrix):
    """
    Plots the overlap matrix as a heatmap with color and saves it as a PNG file.
    The aggregated matrix is used as the annotation matrix.
    """
    # Create a DataFrame for the overlap matrix with gene set names as row/column labels
    overlap_df = pd.DataFrame(overlap_matrix, index=gene_set_names, columns=gene_set_names)

    # Convert the aggregated matrix to a DataFrame (skip the header row and column)
    annotation_matrix = pd.DataFrame(aggregated_matrix[1:, 1:], index=gene_set_names, columns=gene_set_names)

    # Plot the heatmap
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        overlap_df,
        annot=annotation_matrix,  # Use the aggregated matrix as annotations
        fmt="",  # Disable default formatting since we're providing custom text
        cmap='OrRd',  # Color map
        vmin=0,  # Minimum value for color scale
        vmax=1,  # Maximum value for color scale
        linewidths=0.5,  # Add lines between cells
        linecolor="black",
    )

    # Customize the colorbar edge color
    if ax.collections[0].colorbar is not None:
        cbar = ax.collections[0].colorbar
        cbar.outline.set_edgecolor('black')  # Set edge color of the colorbar
        cbar.outline.set_linewidth(0.5)      # Set edge line width

    plt.title(f"Overlap Matrix for {library}", fontsize=14, pad=20)

    # Add padding to xticks and yticks
    ax.tick_params(axis='x', which='both', pad=5)  # Add padding to x-axis tick labels
    ax.tick_params(axis='y', which='both', pad=5)  # Add padding to y-axis tick labels

    # Add axis lines
    for spine in ax.spines.values():  # Enable all spines (borders)
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1)

    # Save the heatmap to a file
    output_file = os.path.join(output_dir, f"overlap_heatmap_{library}.png")
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()
    sys.stderr.write(f"  Overlap heatmap saved to {output_file}\n")


def aggregate_matrices(p_value_matrix, overlap_matrix, output_file, gene_set_names):
    """
    Aggregates the p-value and overlap matrices into a single matrix.
    The diagonal elements are replaced with dashes, the upper triangle contains p-values
    in exponential notation, and the lower triangle contains overlap values.
    The final matrix is saved as a TSV file with gene set names in the initial row and column.
    Returns the aggregated matrix for use in the heatmap.
    """
    # Ensure the matrices are of the same shape
    assert p_value_matrix.shape == overlap_matrix.shape, "Matrices must be of the same shape"
    
    # Create an empty matrix to store the aggregated results
    n = p_value_matrix.shape[0]
    aggregated_matrix = np.empty((n + 1, n + 1), dtype=object)  # +1 for headers
    
    # Add headers for the gene sets
    aggregated_matrix[0, 0] = "set"  # set top-left corner
    for i in range(n):
        aggregated_matrix[0, i + 1] = gene_set_names[i]  # Column headers
        aggregated_matrix[i + 1, 0] = gene_set_names[i]  # Row headers
    
    # Fill the matrix
    for i in range(n):
        for j in range(n):
            if i == j:
                aggregated_matrix[i + 1, j + 1] = '-'  # Replace diagonal with dash
            elif i < j:
                # Upper triangle: p-value in exponential notation
                aggregated_matrix[i + 1, j + 1] = f"{p_value_matrix[i, j]:.1e}"
            else:
                # Lower triangle: overlap with 3 significant digits
                aggregated_matrix[i + 1, j + 1] = f"{overlap_matrix[i, j]:.3f}"
    
    # Save the aggregated matrix to a TSV file
    np.savetxt(output_file, aggregated_matrix, delimiter='\t', fmt='%s')

    # Return the aggregated matrix for use in the heatmap
    return aggregated_matrix


def calculate_overlap_matrix(input_files, library, output_dir, adjusted_p_threshold, background_terms):
    """
    Calculates a 3x3 matrix with the fraction of overlapping significant terms and their statistical significance for a single library.
    """
    # Initialize the overlap fraction matrix and p-value matrix
    overlap_fraction_matrix = np.zeros((3, 3), dtype=float)
    p_value_matrix = np.zeros((3, 3), dtype=float)

    # Read the significant terms for each gene set
    significant_terms = {}
    for i, input_file in enumerate(input_files):
        # Get the basename of the input file (without extension)
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

    # Calculate the overlap fraction matrix and p-value matrix
    for i in range(3):
        for j in range(3):
            set1 = significant_terms[f"gene_list_{i + 1}"]
            set2 = significant_terms[f"gene_list_{j + 1}"]
            overlap = len(set1.intersection(set2))
            min_size = min(len(set1), len(set2))
            overlap_fraction = overlap / min_size if min_size > 0 else 0.0
            overlap_fraction_matrix[i, j] = overlap_fraction

            # Calculate the hypergeometric survival function (p-value)
            M = len(background_terms)  # Number of terms associated with the background genes
            n = len(set1)              # Number of terms in set1
            N = len(set2)              # Number of terms in set2
            k = overlap                # Number of overlapping terms
            p_value = hypergeom.sf(k - 1, M, n, N)  # Survival function (P(X >= k))
            p_value_matrix[i, j] = p_value

    # Get the basenames of the input files for gene set names
    gene_set_names = [os.path.splitext(os.path.basename(file))[0] for file in input_files]

    # Aggregate the matrices and save the result
    aggregated_output_file = os.path.join(output_dir, f"aggregated_matrix_{library}.tsv")
    aggregated_matrix = aggregate_matrices(p_value_matrix, overlap_fraction_matrix, aggregated_output_file, gene_set_names)
    sys.stderr.write(f"Aggregated matrix saved to {aggregated_output_file}\n")

    # Plot the overlap matrix as a heatmap
    plot_overlap_heatmap(overlap_fraction_matrix, gene_set_names, output_dir, library, aggregated_matrix)


def get_background_terms(library, background_genes):
    """
    Returns the set of terms associated with the genes in the background list.
    """
    # Use gseapy to get the library
    library_file=os.path.join(f"{ndd_path}/libs", f"{library}.gmt")
    gene_sets = gp.get_library(library_file)
    #gene_sets = gp.get_library(library)

    # Convert background_genes to a set for faster membership checking
    background_genes_set = set(background_genes)

    # Collect all terms associated with the background genes
    background_terms = set()
    for term, genes in gene_sets.items():
        # Check if any gene in the term's gene list is in the background genes
        if any(gene in background_genes_set for gene in genes):
            background_terms.add(term)

    return background_terms


def generate_summary_table(input_files, libraries, output_dir, adjusted_p_threshold):
    """
    Generates a summary table with significant terms for each library and the fraction of common terms for the intersections.
    The set names are derived using the basename of the input files.
    For individual sets, the total number of terms mapping to the set is added in parentheses.
    For intersections, the fraction of overlapping terms (divided by the minimum number of significant terms) is shown,
    along with the p-value of the overlap in parentheses. The background terms for each library are calculated as the union
    of terms from the three sets.
    
    Parameters:
    - input_files: List of input file paths for the three gene sets.
    - libraries: List of libraries to analyze.
    - output_dir: Directory where the enrichment results are stored.
    - adjusted_p_threshold: Adjusted p-value threshold for significance.
    
    Returns:
    - summary_df: A DataFrame containing the summary table.
    """
    # Extract set names from the input file paths using basename
    set_names = [os.path.splitext(os.path.basename(file))[0] for file in input_files]

    # Initialize a list to store the summary data
    summary_data = []

    # Iterate over each library
    for library in libraries:
        # Read the significant terms for each gene set
        significant_terms = {}
        total_terms = {}
        all_terms = set()  # To store the union of terms from all three sets

        for i, input_file in enumerate(input_files):
            basename = os.path.splitext(os.path.basename(input_file))[0]
            file_name = f"{basename}_{library}.tsv"
            file_path = os.path.join(output_dir, file_name)
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep="\t")
                df_filtered = df[df["Adjusted P-value"] < adjusted_p_threshold]
                significant_terms[f"gene_list_{i + 1}"] = set(df_filtered["Term"])
                total_terms[f"gene_list_{i + 1}"] = len(df["Term"])  # Total terms in the set
                all_terms.update(df["Term"])  # Add terms to the union
            else:
                sys.stderr.write(f"Warning: File {file_path} not found.\n")
                significant_terms[f"gene_list_{i + 1}"] = set()
                total_terms[f"gene_list_{i + 1}"] = 0

        # Get the significant terms and total terms for each gene set
        set1_terms = significant_terms["gene_list_1"]
        set2_terms = significant_terms["gene_list_2"]
        set3_terms = significant_terms["gene_list_3"]

        total_set1_terms = total_terms["gene_list_1"]
        total_set2_terms = total_terms["gene_list_2"]
        total_set3_terms = total_terms["gene_list_3"]

        # Calculate the fraction of overlapping terms and p-values for each intersection
        intersection_12 = len(set1_terms.intersection(set2_terms))
        intersection_13 = len(set1_terms.intersection(set3_terms))
        intersection_23 = len(set2_terms.intersection(set3_terms))

        # Calculate the fraction of overlapping terms (divided by the minimum number of significant terms)
        fraction_12 = intersection_12 / min(len(set1_terms), len(set2_terms)) if min(len(set1_terms), len(set2_terms)) > 0 else 0
        fraction_13 = intersection_13 / min(len(set1_terms), len(set3_terms)) if min(len(set1_terms), len(set3_terms)) > 0 else 0
        fraction_23 = intersection_23 / min(len(set2_terms), len(set3_terms)) if min(len(set2_terms), len(set3_terms)) > 0 else 0

        # Calculate p-values for the intersections using the hypergeometric test
        M = len(all_terms)  # Total number of terms in the union of all three sets
        n1 = len(set1_terms)  # Number of terms in set1
        n2 = len(set2_terms)  # Number of terms in set2
        n3 = len(set3_terms)  # Number of terms in set3

        p_value_12 = hypergeom.sf(intersection_12 - 1, M, n1, n2)
        p_value_13 = hypergeom.sf(intersection_13 - 1, M, n1, n3)
        p_value_23 = hypergeom.sf(intersection_23 - 1, M, n2, n3)

        # Append the data for the current library to the summary list
        summary_data.append([
            library,
            f"{len(set1_terms)} ({total_set1_terms})",
            f"{len(set2_terms)} ({total_set2_terms})",
            f"{len(set3_terms)} ({total_set3_terms})",
            f"{fraction_12:.3f} ({p_value_12:.2e})",
            f"{fraction_13:.3f} ({p_value_13:.2e})",
            f"{fraction_23:.3f} ({p_value_23:.2e})"
        ])

    # Create a DataFrame from the summary data
    summary_df = pd.DataFrame(summary_data, columns=[
        "Library",
        f"Significant Terms ({set_names[0]})",
        f"Significant Terms ({set_names[1]})",
        f"Significant Terms ({set_names[2]})",
        f"Fraction Overlap ({set_names[0]} & {set_names[1]})",
        f"Fraction Overlap ({set_names[0]} & {set_names[2]})",
        f"Fraction Overlap ({set_names[1]} & {set_names[2]})"
    ])

    return summary_df


def main():
    # Get the directory of the script
    script_dir=Path(__file__).parent.parent
    global ndd_path 
    ndd_path = script_dir

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Perform enrichment analysis on gene lists')
    parser.add_argument('gene_list1', help='First gene list file')
    parser.add_argument('gene_list2', help='Second gene list file')
    parser.add_argument('gene_list3', help='Third gene list file')
    parser.add_argument('background_list', help='Background gene list file')
    parser.add_argument('--output_dir', default=f"{ndd_path}/results",
                       help='Output directory for results')
    parser.add_argument('--summary_file', default="summary_table.tsv",
                       help='Name of the summary output file')
    parser.add_argument('--libraries', type=str, default=None,
                       help='Comma-separated list of libraries to use (default: all libraries)')
    args = parser.parse_args()

    # Read gene lists from files
    input_files = [args.gene_list1, args.gene_list2, args.gene_list3]
    gene_lists = [read_gene_list(file) for file in input_files]
    background_genes = read_gene_list(args.background_list)

    # Print the number of genes in each list
    for i, gene_list in enumerate(gene_lists):
        sys.stderr.write(f"Number of genes in {os.path.basename(input_files[i])}: {len(gene_list)}\n")

    # Define the gene set libraries
    if args.libraries:
        libraries = [lib.strip() for lib in args.libraries.split(",")] 
    else:
        libraries = [
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

    # Define the adjusted p-value threshold
    adjusted_p_threshold = 0.01

    # Process each library one at a time
    for library in libraries:
        # Get the terms associated with the background genes
        background_terms = get_background_terms(library, background_genes)
        sys.stderr.write(f"Number of terms associated with background genes in {library}: {len(background_terms)}\n")

        # Perform enrichment analysis for the current library
        perform_enrichment(gene_lists, background_genes, library, args.output_dir, input_files)

    # Generate summary table
    summary_table = generate_summary_table(input_files, libraries, args.output_dir, adjusted_p_threshold)
    summary_output_path = os.path.join(args.output_dir, args.summary_file)
    summary_table.to_csv(summary_output_path, sep="\t", index=False)
    sys.stderr.write(f"Summary table saved to {summary_output_path}\n")

if __name__ == "__main__":
    main()
