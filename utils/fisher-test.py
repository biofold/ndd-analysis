import pandas as pd
import numpy as np
import sys
import os
from scipy.stats import fisher_exact
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

def pairwise_fisher_test(filename, category_col, score_col, output_file, alpha=0.05, save_heatmap=True, save_pdf=True, save_png=False):
    """
    Perform pairwise Fisher's exact tests between categories on a binary (0/1)
    outcome and visualize results, mirroring ks-test.py's combined-heatmap
    convention (p-values below the diagonal, odds ratios above it, same
    threshold-diverging blue/white/red color scale).

    Args:
        filename: Path to input file (CSV, Excel, or TSV)
        category_col: 1-based column index for categories
        score_col: 1-based column index for the binary (0/1) outcome
        alpha: Significance level for hypothesis testing (default: 0.05)
        save_heatmap: Whether to save the heatmap visualization (default: True)
        save_pdf: Whether to save as PDF (default: True)
        save_png: Whether to save as PNG (default: False)
    """

    file_ext = os.path.splitext(filename)[1].lower()

    try:
        if file_ext == '.csv':
            df = pd.read_csv(filename)
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(filename)
        elif file_ext in ['.tsv', '.txt']:
            df = pd.read_csv(filename, sep='\t')
        else:
            df = pd.read_csv(filename, sep=None, engine='python')
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    cat_idx = category_col - 1
    score_idx = score_col - 1

    if cat_idx >= len(df.columns) or score_idx >= len(df.columns):
        print(f"Error: Column index out of range. File has {len(df.columns)} columns.")
        return None

    cat_col_name = df.columns[cat_idx]
    score_col_name = df.columns[score_idx]

    print(f"Category column: '{cat_col_name}' (column {category_col})")
    print(f"Score column (binary 0/1): '{score_col_name}' (column {score_col})")
    print(f"Significance level (alpha): {alpha}")

    df_clean = df[[cat_col_name, score_col_name]].dropna()
    df_clean[cat_col_name] = pd.to_numeric(df_clean[cat_col_name], errors='coerce')
    df_clean[score_col_name] = pd.to_numeric(df_clean[score_col_name], errors='coerce')
    df_clean = df_clean.dropna()

    if not set(df_clean[score_col_name].unique()) <= {0, 1}:
        print("Error: Score column must be binary (0/1).")
        return None

    categories = sorted(df_clean[cat_col_name].unique())
    print(f"\nCategories found (sorted): {categories}")
    print(f"Total data points: {len(df_clean)}")

    category_data = {}
    for cat in categories:
        sub = df_clean[df_clean[cat_col_name] == cat][score_col_name]
        category_data[cat] = (int(sub.sum()), len(sub))
        print(f"Category {cat}: {len(sub)} samples, {int(sub.sum())} positive")

    print("\n" + "="*80)
    print("PAIRWISE FISHER'S EXACT TEST RESULTS")
    print("="*80)

    n_categories = len(categories)
    or_matrix = np.ones((n_categories, n_categories))
    pvalue_matrix = np.ones((n_categories, n_categories))
    significant_matrix = np.zeros((n_categories, n_categories), dtype=bool)
    neg_log_pvalue_matrix = np.zeros((n_categories, n_categories))

    pairwise_results = []

    for i, cat1 in enumerate(categories):
        for j, cat2 in enumerate(categories):
            if i < j:
                k1, n1 = category_data[cat1]
                k2, n2 = category_data[cat2]
                table = [[k1, n1 - k1], [k2, n2 - k2]]
                odds_ratio, p_value = fisher_exact(table)

                epsilon = 1e-300
                neg_log_pvalue = -np.log10(max(p_value, epsilon))

                or_matrix[i, j] = odds_ratio
                or_matrix[j, i] = (1.0 / odds_ratio) if odds_ratio > 0 else np.inf
                pvalue_matrix[i, j] = p_value
                pvalue_matrix[j, i] = p_value
                neg_log_pvalue_matrix[i, j] = neg_log_pvalue
                neg_log_pvalue_matrix[j, i] = neg_log_pvalue

                is_significant = p_value < alpha
                significant_matrix[i, j] = is_significant
                significant_matrix[j, i] = is_significant

                pairwise_results.append({
                    'Category 1': cat1,
                    'Category 2': cat2,
                    'Odds Ratio': odds_ratio,
                    'P-value': p_value,
                    '-log10(P-value)': neg_log_pvalue,
                    'Significant': is_significant,
                    'Significance': '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
                })

    print("\nOdds Ratio Matrix (row category vs column category):")
    print("-"*80)
    header = "Category".ljust(12)
    for cat in categories:
        header += f"{cat:>12}"
    print(header)
    print("-"*80)
    for i, cat in enumerate(categories):
        row = f"{cat:<12}"
        for j in range(n_categories):
            row += f"{'---':>12}" if i == j else f"{or_matrix[i, j]:>12.4f}"
        print(row)

    print("\nP-value Matrix:")
    print("-"*80)
    header = "Category".ljust(12)
    for cat in categories:
        header += f"{cat:>12}"
    print(header)
    print("-"*80)
    for i, cat in enumerate(categories):
        row = f"{cat:<12}"
        for j in range(n_categories):
            row += f"{'---':>12}" if i == j else f"{pvalue_matrix[i, j]:>12.4e}"
        print(row)

    print(f"\nSignificance Matrix (threshold: p < {alpha}):")
    print("-"*80)
    header = "Category".ljust(12)
    for cat in categories:
        header += f"{cat:>12}"
    print(header)
    print("-"*80)
    for i, cat in enumerate(categories):
        row = f"{cat:<12}"
        for j in range(n_categories):
            if i == j:
                row += f"{'---':>12}"
            elif pvalue_matrix[i, j] < alpha:
                row += f"{'*':>12}"
            else:
                row += f"{'ns':>12}"
        print(row)

    print("\nDetailed Pairwise Results:")
    print("-"*80)
    print(f"{'Cat 1':<8} {'Cat 2':<8} {'OR':<10} {'P-value':<12} {'Significant':<12}")
    print("-"*80)
    for result in pairwise_results:
        print(f"{result['Category 1']:<8} {result['Category 2']:<8} "
              f"{result['Odds Ratio']:<10.4f} {result['P-value']:<12.4e} "
              f"{'Yes' if result['Significant'] else 'No':<12}")

    n_tests = len(pairwise_results)
    n_significant = sum(1 for r in pairwise_results if r['Significant'])

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Significance threshold: p < {alpha}")
    print(f"Total pairwise comparisons: {n_tests}")
    print(f"Significant differences: {n_significant}")
    print(f"Non-significant differences: {n_tests - n_significant}")

    bonferroni_alpha = alpha
    n_bonf_significant = sum(1 for r in pairwise_results if r['P-value'] < bonferroni_alpha)
    print(f"\nBonferroni correction:")
    print(f"Adjusted alpha: {bonferroni_alpha:.6f}")
    print(f"Significant after correction: {n_bonf_significant}")

    results_df = pd.DataFrame(pairwise_results)

    if output_file:
        base_filename = os.path.splitext(output_file)[0]
    else:
        base_filename = os.path.splitext(filename)[0]
        output_file = f"{base_filename}_stats.txt"
    pdf_file = f"{base_filename}.pdf"

    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ Results saved to: {output_file}")

    if save_heatmap and n_categories > 1:
        create_combined_heatmap(or_matrix, pvalue_matrix, categories, cat_col_name, score_col_name, pdf_file, alpha, save_pdf, save_png)

    return results_df, or_matrix, pvalue_matrix

def create_combined_heatmap(or_matrix, pvalue_matrix, categories, cat_col_name, score_col_name, pdf_file, alpha=0.05, save_pdf=True, save_png=False):
    """
    Create a combined heatmap with p-values below diagonal and odds ratios
    above diagonal, using the same threshold-diverging colormap as
    ks-test.py's create_combined_heatmap.
    """
    print("\n" + "="*80)
    print("CREATING COMBINED HEATMAP VISUALIZATION")
    print("="*80)

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    fig, ax = plt.subplots(figsize=(11, 10), facecolor='white')
    ax.set_facecolor('white')

    labels = [str(cat) for cat in categories]
    n_categories = len(categories)

    threshold_neg_log = -np.log10(alpha)

    lower_values = []
    for i in range(n_categories):
        for j in range(i):
            p_val = pvalue_matrix[i, j]
            epsilon = 1e-300
            lower_values.append(-np.log10(max(p_val, epsilon)))

    if lower_values:
        min_val = 0
        max_val = max(lower_values)
    else:
        min_val = 0
        max_val = threshold_neg_log * 2

    if max_val > min_val:
        threshold_position = (threshold_neg_log - min_val) / (max_val - min_val)
    else:
        threshold_position = 0.5

    colors = []
    positions = []

    n_blue_steps = 100
    blue_colors = plt.cm.Blues(np.linspace(0.5, 0.0, n_blue_steps))
    blue_colors[-1] = [1, 1, 1, 1]
    blue_positions = np.linspace(0, threshold_position, n_blue_steps)
    for i in range(n_blue_steps):
        colors.append(blue_colors[i])
        positions.append(blue_positions[i])

    colors.append([1, 1, 1, 1])
    positions.append(threshold_position)

    n_red_steps = 100
    red_colors = plt.cm.Reds(np.linspace(0.0, 0.5, n_red_steps))
    red_colors[0] = [1, 1, 1, 1]
    red_positions = np.linspace(threshold_position, 1.0, n_red_steps)
    for i in range(n_red_steps):
        colors.append(red_colors[i])
        positions.append(red_positions[i])

    cmap_lower = mcolors.LinearSegmentedColormap.from_list('threshold_colormap',
                                                            list(zip(positions, colors)))
    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)

    # Lower triangle: p-values
    for i in range(n_categories):
        for j in range(i):
            p_val = pvalue_matrix[i, j]
            neg_log_p = -np.log10(max(p_val, 1e-300))
            color = cmap_lower(norm(neg_log_p))
            rect = Rectangle((j - 0.5, i - 0.5), 1, 1,
                              linewidth=0, edgecolor='none',
                              facecolor=color, alpha=0.9)
            ax.add_patch(rect)
            p_text = f'{p_val:.2e}'
            ax.text(j, i, p_text, ha='center', va='center',
                    fontsize=18, color='black', fontweight='normal')

    # Upper triangle: odds ratios, same color as mirrored lower-triangle cell
    for i in range(n_categories):
        for j in range(i+1, n_categories):
            or_val = or_matrix[i, j]
            p_val = pvalue_matrix[i, j]
            neg_log_p = -np.log10(max(p_val, 1e-300))
            color = cmap_lower(norm(neg_log_p))
            rect = Rectangle((j - 0.5, i - 0.5), 1, 1,
                              linewidth=0, edgecolor='none',
                              facecolor=color, alpha=0.9)
            ax.add_patch(rect)
            or_text = f'{or_val:.2f}'
            ax.text(j, i, or_text, ha='center', va='center',
                    fontsize=18, color='black', fontweight='normal')

    for i in range(n_categories):
        rect = Rectangle((i - 0.5, i - 0.5), 1, 1,
                          linewidth=0, edgecolor='none', facecolor='#E8E8E8')
        ax.add_patch(rect)

    sm = plt.cm.ScalarMappable(cmap=cmap_lower, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('P-value (-log10 scale)', rotation=270, labelpad=20, fontsize=18)

    ax.set_xticks(range(n_categories))
    ax.set_yticks(range(n_categories))
    ax.set_xticklabels(labels, fontsize=18)
    ax.set_yticklabels(labels, fontsize=18)
    ax.set_xlabel(cat_col_name, fontsize=18, labelpad=10)
    ax.set_ylabel(cat_col_name, fontsize=18, labelpad=10)
    ax.set_title(f'Pairwise Fisher\'s Exact Test Results\n{score_col_name} Proportion Comparison',
                 fontsize=18, fontweight='normal', pad=15)

    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')
        spine.set_linewidth(0.8)

    ax.set_xlim(-0.5, n_categories - 0.5)
    ax.set_ylim(n_categories - 0.5, -0.5)

    if n_categories > 6:
        plt.xticks(rotation=45, ha='right')

    plt.tight_layout()

    base_filename = os.path.splitext(pdf_file)[0]

    if save_pdf:
        plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"\n✓ PDF heatmap saved to: {pdf_file}")

    if save_png:
        png_file = f"{base_filename}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ PNG heatmap saved to: {png_file}")

def main():
    """Main function to handle command line arguments"""

    if len(sys.argv) >= 4:
        filename = sys.argv[1]
        category_col = int(sys.argv[2])
        score_col = int(sys.argv[3])
        output_file = sys.argv[4] if len(sys.argv) > 4 else None
        alpha = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
    else:
        print("Pairwise Fisher's Exact Test with Combined Heatmap")
        print("-" * 60)

        filename = input("Enter the input file path: ").strip()

        if not os.path.exists(filename):
            print(f"Error: File '{filename}' does not exist.")
            return

        try:
            category_col = int(input("Enter the category column number (1-based): "))
            score_col = int(input("Enter the binary score column number (1-based): "))

            alpha_input = input("Enter significance threshold (default 0.05): ").strip()
            alpha = float(alpha_input) if alpha_input else 0.05
        except ValueError:
            print("Error: Invalid input.")
            return

    pairwise_fisher_test(filename, category_col, score_col, output_file, alpha, save_heatmap=True, save_pdf=True, save_png=False)

if __name__ == "__main__":
    required_packages = ['pandas', 'numpy', 'scipy', 'matplotlib']

    import subprocess

    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

    main()
