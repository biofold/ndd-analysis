#!/usr/bin/env python3
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import argparse

libs = ['GO_Biological_Process_2026', 
             'GO_Cellular_Component_2026',
             'GO_Molecular_Function_2026', 
             'KEGG_2021_Human', 
             'Reactome_Pathways_2024']

def get_significant_terms(enrichment_file, p_value_threshold=0.01):
    """Reads an enrichment file and returns a set of significant terms."""
    if os.path.exists(enrichment_file):
        df = pd.read_csv(enrichment_file, sep='\t')
        significant_terms = set(df[df['Adjusted P-value'] < p_value_threshold]['Term'])
        return significant_terms
    else:
        raise FileNotFoundError(f"File not found: {enrichment_file}")

def get_genes_with_terms(enrichment_file, terms, p_value_filter=None):
    """Reads an enrichment file and returns a dictionary of genes and their terms."""
    if os.path.exists(enrichment_file):
        df = pd.read_csv(enrichment_file, sep='\t')
        if p_value_filter is not None:
            df = df[df['Adjusted P-value'] < p_value_filter]
        genes_with_terms = {}
        for term in terms:
            if term in df['Term'].values:
                genes = df[df['Term'] == term]['Genes'].str.split(';').explode().tolist()
                for gene in genes:
                    if gene not in genes_with_terms:
                        genes_with_terms[gene] = set()
                    genes_with_terms[gene].add(term)
        return genes_with_terms
    else:
        raise FileNotFoundError(f"File not found: {enrichment_file}")

def get_basename(file_path):
    """Returns the basename of a file (without extension)."""
    return os.path.splitext(os.path.basename(file_path))[0]

def match_terms(set1_file, set2_file, output_dir, p_value_threshold=0.01, p_value_filter=None, libraries=libs):
    """Processes enrichment files and generates the output."""
    set1_basename = get_basename(set1_file)
    set2_basename = get_basename(set2_file)
    all_genes_with_terms = {}

    for library in libraries:
        set1_enrichment_file = os.path.join(output_dir, f"{set1_basename}_{library}.tsv")
        set2_enrichment_file = os.path.join(output_dir, f"{set2_basename}_{library}.tsv")
        significant_terms = get_significant_terms(set2_enrichment_file, p_value_threshold)
        genes_with_terms = get_genes_with_terms(set1_enrichment_file, significant_terms, p_value_filter)

        for gene, terms in genes_with_terms.items():
            if gene not in all_genes_with_terms:
                all_genes_with_terms[gene] = set()
            all_genes_with_terms[gene].add(library)

    return pd.DataFrame({
        'Gene': list(all_genes_with_terms.keys()),
        'Score': [len(terms) for terms in all_genes_with_terms.values()],
        'Terms': ['|'.join(terms) for terms in all_genes_with_terms.values()]
    })

def save_matches(result_table, output_file):
    """Saves the result table to a file."""
    result_table_sorted = result_table.sort_values(by=["Score", "Gene"], ascending=[False, True])
    result_table_sorted.to_csv(output_file, sep='\t', index=False)
    sys.stderr.write(f"Results saved to {output_file}\n")

def plot_scores(output_dir, plot_prefix):
    """Generates a bar plot of the scores."""
    supercandidate_file = os.path.join(output_dir, 'supercandidate.tsv')
    supercandidate_filtered_file = os.path.join(output_dir, 'supercandidate_filtered.tsv')

    supercandidate_df = pd.read_csv(supercandidate_file, sep='\t')
    supercandidate_filtered_df = pd.read_csv(supercandidate_filtered_file, sep='\t')

    sc_scores = supercandidate_df['Score'].value_counts().sort_index()
    scf_scores = supercandidate_filtered_df['Score'].value_counts().sort_index()

    plt.figure(figsize=(10, 6))
    plt.bar(sc_scores.index - 0.2, sc_scores.values, width=0.4, color='#1f77b4')
    plt.bar(scf_scores.index + 0.2, scf_scores.values, width=0.4, color='#ff7f0e')

    plt.xlabel('Score')
    plt.ylabel('Number of Genes')
    plt.title('Score Distribution of Supercandidate Genes')
    plt.xticks(range(1, max(sc_scores.index.max(), scf_scores.index.max()) + 1))
    plt.grid(False)

    plot_file = os.path.join(output_dir, f'{plot_prefix}_scores.png')
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()
    sys.stderr.write(f"Bar plot saved to {plot_file}\n")

def plot_cumulative_scores(output_dir, plot_prefix, threshold=1):
    """Generates a cumulative bar plot of the scores."""
    supercandidate_file = os.path.join(output_dir, 'supercandidate.tsv')
    supercandidate_filtered_file = os.path.join(output_dir, 'supercandidate_filtered.tsv')

    sc_df = pd.read_csv(supercandidate_file, sep='\t')
    scf_df = pd.read_csv(supercandidate_filtered_file, sep='\t')

    sc_scores = sc_df[sc_df['Score'] >= threshold]['Score'].value_counts().sort_index()
    scf_scores = scf_df[scf_df['Score'] >= threshold]['Score'].value_counts().sort_index()

    sc_cumulative = sc_scores.sort_index(ascending=False).cumsum().sort_index()
    scf_cumulative = scf_scores.sort_index(ascending=False).cumsum().sort_index()

    plt.figure(figsize=(10, 6))
    plt.bar(sc_cumulative.index - 0.2, sc_cumulative.values, width=0.4, color='#1f77b4')
    plt.bar(scf_cumulative.index + 0.2, scf_cumulative.values, width=0.4, color='#ff7f0e')

    plt.xlabel('Score')
    plt.ylabel('Cumulative Number of Genes')
    plt.title(f'Cumulative Score Distribution (Score ≥ {threshold})')
    plt.xticks(range(threshold, max(sc_cumulative.index.max(), scf_cumulative.index.max()) + 1))
    plt.grid(False)

    plot_file = os.path.join(output_dir, f'{plot_prefix}_cumulative.png')
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()
    sys.stderr.write(f"Cumulative plot saved to {plot_file}\n")

def main():
    parser = argparse.ArgumentParser(description='Generate supercandidate gene lists and plots')
    parser.add_argument('set1_file', help='First gene set file')
    parser.add_argument('set2_file', help='Second gene set file')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    parser.add_argument('--output', default='supercandidate', help='Base name for output files')
    parser.add_argument('--plot', help='Prefix for plot filenames (plots will be generated only if this is provided)')
    parser.add_argument('--libraries', type=str, default=None,
                       help='Comma-separated list of libraries to use (default: all libraries)')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Libraries
    if args.libraries:
        libs = [lib.strip() for lib in args.libraries.split(",")]
    else:
        libs = ['GO_Biological_Process_2026',
                     'GO_Cellular_Component_2026',
                     'GO_Molecular_Function_2026',
                     'KEGG_2021_Human',
                     'Reactome_Pathways_2024']

    # Generate and save results
    results_all = match_terms(args.set1_file, args.set2_file, args.output_dir, libraries=libs)
    results_filter = match_terms(args.set1_file, args.set2_file, args.output_dir, 
                                 p_value_threshold=0.01, libraries=libs)

    save_matches(results_all, os.path.join(args.output_dir, f'{args.output}.tsv'))
    save_matches(results_filter, os.path.join(args.output_dir, f'{args.output}_filtered.tsv'))

    # Generate plots only if --plot is provided
    if args.plot:
        plot_scores(args.output_dir, args.plot)
        plot_cumulative_scores(args.output_dir, args.plot)

if __name__ == '__main__':
    main()
