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

# Selection rule: a term is significant when its adjusted p-value, at the full
# precision written by gseapy, is strictly below the threshold. On unrounded
# values '<' and '<=' select the same terms (no adjusted p-value equals 0.01
# exactly), but on ROUNDED values the MOE score changes: rounding to 3 decimals
# drops terms just below 0.01 with '<' (e.g. 0.009994 -> 0.010) and admits terms
# just above it with '<=' (0.01008 -> 0.010) -- 26 / 29 candidates change.
# read_enrichment() therefore refuses tables whose p-values look rounded.
PVAL_COL = 'Adjusted P-value'
ROUNDED_MAX_SIGDIGITS = 6      # full-precision floats are written with 15-20
DECISION_ZONE = (0.005, 0.02)  # reported separately: rounding here can cross the cutoff


def _significant_digits(text):
    """Number of significant digits in a numeric string (e.g. '0.0100' -> 1)."""
    m = text.strip().lower().lstrip('+-').split('e')[0].replace('.', '').lstrip('0').rstrip('0')
    return len(m)


def check_pvalue_precision(raw, path):
    """Stop if any adjusted p-value looks rounded (policy: no rounding anywhere).

    raw: the adjusted p-value column read as text, exactly as stored. Exact 0
    (underflow) and 1 (the Benjamini-Hochberg cap) are legitimately short.
    """
    raw = raw.dropna().astype(str)
    vals = pd.to_numeric(raw, errors='coerce')
    ok = vals.notna() & (vals > 0) & (vals < 1)
    raw, vals = raw[ok], vals[ok]
    if raw.empty:
        return
    short = raw.map(_significant_digits) <= ROUNDED_MAX_SIGDIGITS
    if short.any():
        lo, hi = DECISION_ZONE
        n_zone = int((short & (vals >= lo) & (vals <= hi)).sum())
        ex = ", ".join(raw[short].head(5))
        sys.exit(f"Error: {path}: {int(short.sum())} adjusted p-values look rounded "
                 f"(<= {ROUNDED_MAX_SIGDIGITS} significant digits, e.g. {ex}; {n_zone} of them "
                 f"between {lo} and {hi}). Selecting terms on rounded values changes MOE "
                 f"scores; regenerate the table at full precision (scripts/1_enrichr_all.py).")


def read_enrichment(enrichment_file):
    """Read an enrichment table, checking the stored precision of its p-values."""
    df = pd.read_csv(enrichment_file, sep='\t', dtype={PVAL_COL: str})
    check_pvalue_precision(df[PVAL_COL], enrichment_file)
    # float() is the exact decimal->double conversion; pd.to_numeric and read_csv's
    # default parser can be one ulp off (9,238 of 16,114 HC values), see policy above.
    df[PVAL_COL] = df[PVAL_COL].map(lambda x: float(x) if isinstance(x, str) else x)
    return df


def get_significant_terms(enrichment_file, p_value_threshold=0.01):
    """Reads an enrichment file and returns a set of significant terms."""
    if os.path.exists(enrichment_file):
        df = read_enrichment(enrichment_file)
        significant_terms = set(df[df['Adjusted P-value'] < p_value_threshold]['Term'])
        return significant_terms
    else:
        raise FileNotFoundError(f"File not found: {enrichment_file}")

def get_genes_with_terms(enrichment_file, terms, p_value_filter=None):
    """Reads an enrichment file and returns a dictionary of genes (uppercase) and their terms."""
    if os.path.exists(enrichment_file):
        df = read_enrichment(enrichment_file)
        if p_value_filter is not None:
            df = df[df['Adjusted P-value'] < p_value_filter]
        genes_with_terms = {}
        for term in terms:
            if term in df['Term'].values:
                genes = df[df['Term'] == term]['Genes'].str.split(';').explode().tolist()
                for gene in genes:
                    gene_upper = gene.strip().upper()
                    if gene_upper not in genes_with_terms:
                        genes_with_terms[gene_upper] = set()
                    genes_with_terms[gene_upper].add(term)
        return genes_with_terms
    else:
        raise FileNotFoundError(f"File not found: {enrichment_file}")

def get_all_genes_from_file(gene_file):
    """Reads a gene list file and returns a dict mapping uppercase gene -> original gene."""
    if os.path.exists(gene_file):
        gene_map = {}
        with open(gene_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                gene = line.split('\t')[0].strip()
                if gene:
                    gene_map[gene.upper()] = gene   # store original as value
        return gene_map
    else:
        raise FileNotFoundError(f"File not found: {gene_file}")

def get_basename(file_path):
    """Returns the basename of a file (without extension)."""
    return os.path.splitext(os.path.basename(file_path))[0]

def match_terms(set1_file, set2_file, output_dir, p_value_threshold=0.01, p_value_filter=None, libraries=libs, include_zero_scores=True):
    """Processes enrichment files and generates the output."""
    set1_basename = get_basename(set1_file)
    set2_basename = get_basename(set2_file)
    all_genes_with_terms = {}   # keys are uppercase gene names
    
    # Get the mapping from uppercase to original gene name
    gene_map = get_all_genes_from_file(set1_file)
    
    for library in libraries:
        set1_enrichment_file = os.path.join(output_dir, f"{set1_basename}_{library}.tsv")
        set2_enrichment_file = os.path.join(output_dir, f"{set2_basename}_{library}.tsv")
        significant_terms = get_significant_terms(set2_enrichment_file, p_value_threshold)
        genes_with_terms = get_genes_with_terms(set1_enrichment_file, significant_terms, p_value_filter)
        
        for gene_upper, terms in genes_with_terms.items():
            # Only include genes that are present in the original input list
            if gene_upper in gene_map:
                if gene_upper not in all_genes_with_terms:
                    all_genes_with_terms[gene_upper] = set()
                all_genes_with_terms[gene_upper].add(library)
    
    # Include genes with score 0 if requested
    if include_zero_scores:
        for gene_upper in gene_map:   # gene_map keys are uppercase
            if gene_upper not in all_genes_with_terms:
                all_genes_with_terms[gene_upper] = set()
    
    # Build output DataFrame with original gene names
    genes_output = []
    for gene_upper in all_genes_with_terms:
        original_gene = gene_map.get(gene_upper, gene_upper)
        genes_output.append(original_gene)
    
    return pd.DataFrame({
        'Gene': genes_output,
        'Score': [len(all_genes_with_terms[gene_upper]) for gene_upper in all_genes_with_terms],
        'Terms': ['|'.join(sorted(all_genes_with_terms[gene_upper])) if all_genes_with_terms[gene_upper] else '' for gene_upper in all_genes_with_terms]
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
    plt.bar(sc_scores.index - 0.2, sc_scores.values, width=0.4, color='#1f77b4', label='All')
    plt.bar(scf_scores.index + 0.2, scf_scores.values, width=0.4, color='#ff7f0e', label='Filtered')

    plt.xlabel('Score')
    plt.ylabel('Number of Genes')
    plt.title('Score Distribution of Supercandidate Genes')
    plt.xticks(range(0, max(sc_scores.index.max(), scf_scores.index.max()) + 1))
    plt.legend()
    plt.grid(False)

    plot_file = os.path.join(output_dir, f'{plot_prefix}_scores.png')
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()
    sys.stderr.write(f"Bar plot saved to {plot_file}\n")

def plot_cumulative_scores(output_dir, plot_prefix, threshold=0):
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
    plt.bar(sc_cumulative.index - 0.2, sc_cumulative.values, width=0.4, color='#1f77b4', label='All')
    plt.bar(scf_cumulative.index + 0.2, scf_cumulative.values, width=0.4, color='#ff7f0e', label='Filtered')

    plt.xlabel('Score')
    plt.ylabel('Cumulative Number of Genes')
    plt.title(f'Cumulative Score Distribution (Score ≥ {threshold})')
    plt.xticks(range(threshold, max(sc_cumulative.index.max(), scf_cumulative.index.max()) + 1))
    plt.legend()
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
    parser.add_argument('--include-zero-scores', action='store_true', default=True,
                       help='Include genes with score 0 (default: True)')
    parser.add_argument('--exclude-zero-scores', action='store_true',
                       help='Exclude genes with score 0')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.libraries:
        libs = [lib.strip() for lib in args.libraries.split(",")]
    else:
        libs = ['GO_Biological_Process_2026',
                'GO_Cellular_Component_2026',
                'GO_Molecular_Function_2026',
                'KEGG_2021_Human',
                'Reactome_Pathways_2024']

    include_zero_scores = not args.exclude_zero_scores

    results_all = match_terms(args.set1_file, args.set2_file, args.output_dir, 
                              libraries=libs, include_zero_scores=include_zero_scores)
    results_filter = match_terms(args.set1_file, args.set2_file, args.output_dir, 
                                 p_value_threshold=0.01, p_value_filter=0.01,
                                 libraries=libs, include_zero_scores=include_zero_scores)

    save_matches(results_all, os.path.join(args.output_dir, f'{args.output}.tsv'))
    save_matches(results_filter, os.path.join(args.output_dir, f'{args.output}_filtered.tsv'))

    print(f"Total genes in {args.output}.tsv: {len(results_all)}")
    print(f"  Score distribution:")
    print(results_all['Score'].value_counts().sort_index())
    print(f"\nTotal genes in {args.output}_filtered.tsv: {len(results_filter)}")
    print(f"  Score distribution:")
    print(results_filter['Score'].value_counts().sort_index())

    if args.plot:
        plot_scores(args.output_dir, args.plot)
        plot_cumulative_scores(args.output_dir, args.plot)

if __name__ == '__main__':
    main()
