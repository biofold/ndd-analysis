#!/usr/bin/env python3
import argparse
import sys
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def read_gmt(gmt_file):
    """Read GMT file and return dictionary of {term: set(genes)}"""
    term_genes = {}
    with open(gmt_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                term = parts[0]
                genes = set(g.upper() for g in parts[2:])  # Skip description field
                term_genes[term] = genes
    return term_genes

def process_enrichment(enrichment_file, p_threshold=0.01):
    """Read enrichment results and return significant terms (using 5th column for adjusted p-value)"""
    significant_terms = set()
    with open(enrichment_file) as f:
        header = next(f, None)  # Skip header
        for line_num, line in enumerate(f, 2):
            if line.strip():
                parts = line.strip().split('\t')
                try:
                    if len(parts) >= 5:
                        adj_pvalue = float(parts[4])  # 5th column (0-indexed)
                        if adj_pvalue <= p_threshold:
                            significant_terms.add(parts[1])  # Term is in 2nd column
                except (IndexError, ValueError) as e:
                    print(f'WARNING: Line {line_num} - {str(e)}', file=sys.stderr)
    return significant_terms

def create_paired_barplot(x_values, y1_values, y2_values, save_path, color1='#1f77b4', color2='#ff7f0e', ymax=2000,  dpi=300):
    """Create a paired bar plot comparing two sets of values"""
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    plt.rcParams.update({
        "font.weight": "bold",
        "xtick.labelsize": 15,
    })

    # Set positions for the bars
    x = np.arange(len(x_values))
    width = 0.4  # Width of the bars

    # Create bars
    ax.bar(x - width/2, y1_values, width, color=color1) #label='Genes in terms'
    ax.bar(x + width/2, y2_values, width, color=color2) #label='Genes not in terms'

    # Add labels and title
    ax.set_xticks(x)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xticklabels(x_values)
    #ax.set_ylabel('Gene Count', fontsize=14, fontweight='bold')
    #ax.set_xlabel('Score Group', fontsize=14, fontweight='bold')
    #ax.legend(fontsize=12)

    if ymax is not None:
        ax.set_ylim(top=ymax)
    # Adjust layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
    #print(f"Figure saved to {save_path}", file=sys.stderr)
    plt.close()

def process_files(score_file, background_file, enrichment_file, gmt_file, p_threshold=0.01, output_file=None, plot_file=None, has_header=True, verbose=False):
    # Read scored genes (convert to uppercase)
    gene_scores = {}
    score_groups = defaultdict(set)
    
    with open(score_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        if has_header:  # Skip header if present
            header = next(reader,None)
        for line_num, line in enumerate(reader, 1):
            if line and len(line) >= 2:
                gene, score = line[0], line[1]
                gene_upper = gene.upper()
                try:
                    score_float = float(score)
                    gene_scores[gene_upper] = score_float
                    score_groups[score_float].add(gene_upper)
                except ValueError:
                    print(f'WARNING: Invalid score format in line {line_num} - {line}', file=sys.stderr)   
 
    # Process background genes (convert to uppercase)
    with open(background_file, 'r') as f:
        background_genes = set(line.strip().upper() for line in f if line.strip())
    
    # Assign score=0 to missing genes
    missing_genes = background_genes - set(gene_scores.keys())
    for gene in missing_genes:
        gene_scores[gene] = 0.0
        score_groups[0.0].add(gene)
    
    # Process enrichment results and GMT file
    significant_terms = process_enrichment(enrichment_file, p_threshold)
    term_genes = read_gmt(gmt_file)
    
    # Create sets to track unique genes per score group
    score_term_genes = defaultdict(set)  # {score: set(genes_in_terms)}
    term_score_genes = defaultdict(lambda: defaultdict(set))  # {term: {score: set(genes)}}
    
    # First collect all genes that appear in ANY significant term
    all_term_genes = set()
    for term in significant_terms:
        if term in term_genes:
            all_term_genes.update(term_genes[term])
    
    # Then map genes to their score groups
    for gene in all_term_genes:
        if gene in gene_scores:
            score = gene_scores[gene]
            score_term_genes[score].add(gene)
    
    # Now build term-specific counts (only counting each gene once)
    for term in significant_terms:
        if term in term_genes:
            for gene in term_genes[term]:
                if gene in gene_scores:
                    score = gene_scores[gene]
                    term_score_genes[term][score].add(gene)
    
    # Prepare cumulative counts for each column
    sorted_scores = sorted(score_groups.keys(), reverse=True)
    
    # Initialize cumulative counts
    cum_in_terms = 0
    cum_not_in_terms = 0
    cum_total = 0
    
    # Prepare data for plotting
    plot_scores = []
    plot_in_terms = []
    plot_not_in_terms = []
    plot_total = []
    
    # Prepare output
    output = []
    output.append("Score\tGenes_in_terms\tGenes_not_in_terms\tTotal_genes\tCum_in_terms\tCum_not_in_terms\tCum_total")
    
    for score in sorted_scores:
        total_genes = len(score_groups[score])
        term_genes_count = len(score_term_genes.get(score, set()))
        nterm_genes = total_genes - term_genes_count
        
        # Update cumulative counts
        cum_in_terms += term_genes_count
        cum_not_in_terms += nterm_genes
        cum_total += total_genes
        
        output.append(f"{score:.0f}\t{term_genes_count}\t{nterm_genes}\t{total_genes}\t{cum_in_terms}\t{cum_not_in_terms}\t{cum_total}")
        
        # Store data for plotting
        plot_scores.append(str(int(score)))
        plot_in_terms.append(term_genes_count)
        plot_not_in_terms.append(nterm_genes)
        plot_total.append(total_genes)

        # Detailed term breakdown per score group
        if verbose and score in score_term_genes:
            output.append("  Genes in terms:")
            for gene in sorted(score_term_genes[score]):
                output.append(f"    {gene}")
    
    # Generate plot if requested
    if plot_file:

        # Reverse the plotting order
        plot_scores.reverse()
        plot_in_terms.reverse()
        plot_not_in_terms.reverse()
        plot_total.reverse()

        create_paired_barplot(
            x_values=plot_scores,
            y1_values=plot_total,
            y2_values=plot_in_terms,
            save_path=plot_file
        )
        #print(f"Plot saved to {plot_file}", file=sys.stderr)
    
    # 2. Term-centric view
    if verbose:
        output.append("\n# Term Analysis (adj. p <= {p_threshold}):")
        for term in sorted(term_score_genes.keys()):
            total_genes = sum(len(genes) for genes in term_score_genes[term].values())
            output.append(f"\n{term}: {total_genes} unique genes")
            for score in sorted(term_score_genes[term].keys(), reverse=True):
                count = len(term_score_genes[term][score])
                output.append(f"  - Score {score}: {count} genes")
    
    # Write output
    if output_file:
        with open(output_file, 'w') as out:
            out.write('\n'.join(output) + '\n')
        #print(f"Results saved to {output_file}", file=sys.stderr)
    else:
        print('\n'.join(output))
    
    # Print summary to stderr
    print(f"\nAnalysis complete:", file=sys.stderr)
    print(f"- Processed {len(gene_scores)} total genes", file=sys.stderr)
    print(f"- Found {len(significant_terms)} significant terms at adj. p <= {p_threshold}", file=sys.stderr)
    print(f"- {len(all_term_genes & gene_scores.keys())} unique genes in significant terms", file=sys.stderr)
    print(f"- {len(score_groups)} score groups analyzed", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Analyze gene score groups with enrichment terms and cumulative counts',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('score_file', help='File with gene<tab>score pairs')
    parser.add_argument('background_file', help='File with full gene list (one per line)')
    parser.add_argument('enrichment_file', help='Enrichment analysis results file')
    parser.add_argument('gmt_file', help='GMT format gene set file')
    parser.add_argument('-p', '--p_threshold', type=float, default=0.01,
                       help='Adjusted p-value threshold (column 5)')
    parser.add_argument('-o', '--output', help='Output file for tabular data (default: stdout)', default=None)
    parser.add_argument('--plot', help='Output file for bar plot (optional)', default=None)
    parser.add_argument("--no_header", help="Input file has no header row", action="store_true")
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Show detailed gene/term breakdown')
    
    args = parser.parse_args()
    
    process_files(
        args.score_file,
        args.background_file,
        args.enrichment_file,
        args.gmt_file,
        args.p_threshold,
        args.output,
        args.plot,
        has_header=not args.no_header,
        verbose=args.verbose
    )
