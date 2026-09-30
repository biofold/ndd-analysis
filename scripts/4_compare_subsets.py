import sys
import csv
import numpy as np
import scipy.stats as stats
from itertools import combinations
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

def plot_horizontal_dendrogram(distance_matrix, labels, output_file, method='average'):
    """Plot and save a horizontal dendrogram using distance matrix"""
    condensed_dist = squareform(distance_matrix)
    Z = linkage(condensed_dist, method=method)

    plt.figure(figsize=(10, 8))
    dendrogram(
        Z,
        labels=labels,
        orientation='left',
        leaf_font_size=16,
    )
    plt.xlabel('Distance', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Dendrogram saved to {output_file}", file=sys.stderr)

def compare_annotations(
    file_path, 
    annotated_col, 
    non_annotated_col, 
    output_file=None,
    label_col=None, 
    delimiter='\t',
    matrix_file=None,  # Changed to accept filename
    use_chi2=False,
    two_tailed=False,
    dendrogram_file=None,
    has_header=True
):
    # Convert to 0-based indexing
    annotated_col -= 1
    non_annotated_col -= 1
    if label_col is not None:
        label_col -= 1

    # Read data
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        if has_header:
            header = next(reader, None)
        data = list(reader)
    
    # Extract columns
    labels = [row[label_col] for row in data] if label_col is not None else [f"Subgroup_{i}" for i in range(len(data))]
    annotated = [int(row[annotated_col]) for row in data]
    non_annotated = [int(row[non_annotated_col]) for row in data]
    n = len(labels)

    order = sorted(range(n), key=lambda idx: labels[idx])
    labels   = [labels[idx] for idx in order]
    annotated = [annotated[idx] for idx in order]
    non_annotated = [non_annotated[idx] for idx in order]  
 
    # Initialize matrices
    pairwise_results = []
    p_matrix = np.full((n, n), np.nan)
    or_matrix = np.full((n, n), np.nan)
    dist_matrix = np.zeros((n, n))

    # Perform all pairwise comparisons
    for i, j in combinations(range(n), 2):
        table = [
            [annotated[j], non_annotated[j]],
            [annotated[i], non_annotated[i]]
        ]
        odds_ratio, _ = stats.fisher_exact(table)
        
        if use_chi2:
            _, p_val, _, _ = stats.chi2_contingency(table)
            test_type = "Chi2"
        else:
            if two_tailed:
                _, p_val = stats.fisher_exact(table)
                test_type = "Fisher_two-tailed"
            else:
                if odds_ratio > 1:
                    _, p_val = stats.fisher_exact(table, alternative='greater')
                    test_type = "Fisher_greater"
                elif odds_ratio < 1:
                    _, p_val = stats.fisher_exact(table, alternative='less')
                    test_type = "Fisher_less"
                else:
                    p_val = 1.0
                    test_type = "Fisher_equal"
        
        pairwise_results.append([
            labels[i], labels[j], f"{p_val:.3e}", f"{odds_ratio:.3f}", test_type
        ])
        
        or_matrix[i, j] = odds_ratio
        p_matrix[j, i] = p_val
        distance = -np.log10(p_val) if p_val > 0 else 10
        dist_matrix[i,j] = dist_matrix[j,i] = distance
    
    np.fill_diagonal(dist_matrix, 0)
    
    if dendrogram_file:
        plot_horizontal_dendrogram(dist_matrix, labels, dendrogram_file)
    
    # Handle matrix output to file
    
    if matrix_file:
        with open(matrix_file, 'w') as f:
            # Write header
            f.write("Group\t" + "\t".join(labels) + "\n")
            # Write matrix rows
            for i in range(n):
                row = [labels[i]]
                for j in range(n):
                    if i < j:
                        val = or_matrix[i,j]
                        row.append(repr(float(val)) if not np.isnan(val) else "NA")
                    elif i > j:
                        val = p_matrix[i,j]
                        row.append(repr(float(val)) if not np.isnan(val) else "NA")
                    else:
                        row.append("NA")
                f.write("\t".join(row) + "\n")
        print(f"Matrix saved to {matrix_file}", file=sys.stderr)
    
    # Always print pairwise results to stdout
    if output_file:
        with open(output_file, 'w') as f:
            # Write header
            f.write("Group1\tGroup2\tP-value\tOddsRatio\tTest\n")
            for row in pairwise_results:
                f.write("\t".join(row) + '\n')
    else:
        print("Group1\tGroup2\tP-value\tOddsRatio\tTest")
        for row in pairwise_results:
            print("\t".join(row))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Compare annotation frequencies between groups",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("file_path", help="Input file path")
    parser.add_argument("annotated_col", type=int, help="1-based column for annotated counts")
    parser.add_argument("non_annotated_col", type=int, help="1-based column for non-annotated counts")
    parser.add_argument("-o", "--output_file", help="Output file for the results of the statistical analysis.")
    parser.add_argument("--label_col", type=int, help="1-based column for group labels")
    parser.add_argument("--delimiter", help="Field delimiter", default="\t")
    parser.add_argument("--matrix", help="Output file for matrix (will be created if specified)")
    parser.add_argument("--chi2", help="Use Chi-square test", action="store_true")
    parser.add_argument("--two_tailed", help="Use two-tailed Fisher test", action="store_true")
    parser.add_argument("--dendrogram", help="Output file for dendrogram (e.g., clusters.png)")
    parser.add_argument("--no_header", help="Input file has no header row", action="store_true")
    args = parser.parse_args()
    
    compare_annotations(
        args.file_path,
        args.annotated_col,
        args.non_annotated_col,
        args.output_file,
        args.label_col,
        args.delimiter,
        args.matrix,  # Now passing filename
        args.chi2,
        args.two_tailed,
        args.dendrogram,
        has_header=not args.no_header
    )
