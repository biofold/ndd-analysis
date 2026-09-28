import pandas as pd
import numpy as np
import sys
import os
from scipy import stats
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

def pairwise_ks_test(filename, category_col, score_col, output_file, alpha=0.05, save_heatmap=True, save_pdf=True, save_png=False):
    """
    Perform pairwise Kolmogorov-Smirnov tests between categories and visualize results.
    
    Args:
        filename: Path to input file (CSV, Excel, or TSV)
        category_col: 1-based column index for categories
        score_col: 1-based column index for scores
        alpha: Significance level for hypothesis testing (default: 0.05)
        save_heatmap: Whether to save the heatmap visualization (default: True)
        save_pdf: Whether to save as PDF (default: True)
        save_png: Whether to save as PNG (default: False)
    """
    
    # Read the file based on extension
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
    
    # Convert 1-based to 0-based indexing
    cat_idx = category_col - 1
    score_idx = score_col - 1
    
    # Check if indices are valid
    if cat_idx >= len(df.columns) or score_idx >= len(df.columns):
        print(f"Error: Column index out of range. File has {len(df.columns)} columns.")
        return None
    
    # Get column names
    cat_col_name = df.columns[cat_idx]
    score_col_name = df.columns[score_idx]
    
    print(f"Category column: '{cat_col_name}' (column {category_col})")
    print(f"Score column: '{score_col_name}' (column {score_col})")
    print(f"Significance level (alpha): {alpha}")
    
    # Clean the data
    df_clean = df[[cat_col_name, score_col_name]].dropna()
    df_clean[cat_col_name] = pd.to_numeric(df_clean[cat_col_name], errors='coerce')
    df_clean[score_col_name] = pd.to_numeric(df_clean[score_col_name], errors='coerce')
    df_clean = df_clean.dropna()
    
    # Get unique categories and sort them
    categories = sorted(df_clean[cat_col_name].unique())
    print(f"\nCategories found (sorted): {categories}")
    print(f"Total data points: {len(df_clean)}")
    
    # Store data for each category
    category_data = {}
    for cat in categories:
        category_data[cat] = df_clean[df_clean[cat_col_name] == cat][score_col_name].values
        print(f"Category {cat}: {len(category_data[cat])} samples")
    
    # Perform pairwise KS tests
    print("\n" + "="*80)
    print("PAIRWISE KOLMOGOROV-SMIRNOV TEST RESULTS")
    print("="*80)
    
    # Create results matrices
    n_categories = len(categories)
    ks_matrix = np.zeros((n_categories, n_categories))
    pvalue_matrix = np.ones((n_categories, n_categories))
    significant_matrix = np.zeros((n_categories, n_categories), dtype=bool)
    neg_log_pvalue_matrix = np.zeros((n_categories, n_categories))
    
    # Store all pairwise results
    pairwise_results = []
    
    # Perform all pairwise comparisons
    for i, cat1 in enumerate(categories):
        for j, cat2 in enumerate(categories):
            if i < j:  # Only unique pairs
                data1 = category_data[cat1]
                data2 = category_data[cat2]
                
                # Perform KS test
                ks_stat, p_value = stats.ks_2samp(data1, data2)
                
                # Calculate -log10(p-value) for heatmap
                epsilon = 1e-300
                neg_log_pvalue = -np.log10(max(p_value, epsilon))
                
                # Store in matrices
                ks_matrix[i, j] = ks_stat
                ks_matrix[j, i] = ks_stat
                pvalue_matrix[i, j] = p_value
                pvalue_matrix[j, i] = p_value
                neg_log_pvalue_matrix[i, j] = neg_log_pvalue
                neg_log_pvalue_matrix[j, i] = neg_log_pvalue
                
                # Check significance
                is_significant = p_value < alpha
                significant_matrix[i, j] = is_significant
                significant_matrix[j, i] = is_significant
                
                pairwise_results.append({
                    'Category 1': cat1,
                    'Category 2': cat2,
                    'KS Statistic': ks_stat,
                    'P-value': p_value,
                    '-log10(P-value)': neg_log_pvalue,
                    'Significant': is_significant,
                    'Significance': '***' if p_value < 0.001 else '**' if p_value < 0.01 else '*' if p_value < 0.05 else 'ns'
                })
    
    # Print KS statistic matrix
    print("\nKS Statistic Matrix:")
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
            else:
                row += f"{ks_matrix[i, j]:>12.4f}"
        print(row)
    
    # Print p-value matrix
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
            if i == j:
                row += f"{'---':>12}"
            else:
                row += f"{pvalue_matrix[i, j]:>12.4e}"
        print(row)
    
    # Print significance matrix
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
    
    # Print detailed pairwise results table
    print("\nDetailed Pairwise Results:")
    print("-"*80)
    print(f"{'Cat 1':<8} {'Cat 2':<8} {'KS Stat':<10} {'P-value':<12} {'Significant':<12}")
    print("-"*80)
    
    for result in pairwise_results:
        cat1 = result['Category 1']
        cat2 = result['Category 2']
        ks_stat = result['KS Statistic']
        p_value = result['P-value']
        is_sig = "Yes" if result['Significant'] else "No"
        
        print(f"{cat1:<8} {cat2:<8} {ks_stat:<10.4f} {p_value:<12.4e} {is_sig:<12}")
    
    # Summary statistics
    n_tests = len(pairwise_results)
    n_significant = sum(1 for r in pairwise_results if r['Significant'])
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Significance threshold: p < {alpha}")
    print(f"Total pairwise comparisons: {n_tests}")
    print(f"Significant differences: {n_significant}")
    print(f"Non-significant differences: {n_tests - n_significant}")
    
    # Bonferroni correction
    #bonferroni_alpha = alpha / n_tests
    bonferroni_alpha = alpha
    n_bonf_significant = sum(1 for r in pairwise_results if r['P-value'] < bonferroni_alpha)
    
    print(f"\nBonferroni correction:")
    print(f"Adjusted alpha: {bonferroni_alpha:.6f}")
    print(f"Significant after correction: {n_bonf_significant}")
    
    # Create a DataFrame for export
    results_df = pd.DataFrame(pairwise_results)
    
    # Save results to CSV
    if output_file:
        base_filename = os.path.splitext(output_file)[0]
    else:
        base_filename = os.path.splitext(filename)[0]
        output_file = f"{base_filename}_stats.txt"
    pdf_file = f"{base_filename}.pdf"

    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ Results saved to: {output_file}")
    
    # Create and save heatmap
    if save_heatmap and n_categories > 1:
        create_combined_heatmap(ks_matrix, pvalue_matrix, categories, cat_col_name, score_col_name, pdf_file, alpha, save_pdf, save_png)
    
    return results_df, ks_matrix, pvalue_matrix

def create_combined_heatmap(ks_matrix, pvalue_matrix, categories, cat_col_name, score_col_name, pdf_file, alpha=0.05, save_pdf=True, save_png=False):
    """
    Create a combined heatmap with p-values below diagonal and KS statistics above diagonal.
    
    Args:
        ks_matrix: Matrix of KS statistics
        pvalue_matrix: Matrix of p-values
        categories: List of category labels
        cat_col_name: Name of category column
        score_col_name: Name of score column
        alpha: Significance threshold
        save_pdf: Whether to save as PDF
        save_png: Whether to save as PNG
    """
    print("\n" + "="*80)
    print("CREATING COMBINED HEATMAP VISUALIZATION")
    print("="*80)
    
    # Set PDF-friendly font settings
    plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    
    # Create figure with white background
    fig, ax = plt.subplots(figsize=(11, 10), facecolor='white')
    
    # Set plot area background to white
    ax.set_facecolor('white')
    
    labels = [str(cat) for cat in categories]
    n_categories = len(categories)
    
    # Calculate the -log10 of the threshold
    threshold_neg_log = -np.log10(alpha)
    
    # Calculate min and max for color scale
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
    
    # Create custom colormap with lighter colors converging to white at threshold
    # Below threshold: Medium blue -> White (converging to white at threshold)
    # Above threshold: White -> Medium red (diverging from white at threshold)
    
    # Calculate the relative position of the threshold
    if max_val > min_val:
        threshold_position = (threshold_neg_log - min_val) / (max_val - min_val)
    else:
        threshold_position = 0.5
    
    # Create a list of colors that properly converge to white
    colors = []
    positions = []
    
    # Non-significant part (0 to threshold): Medium blue converging to white
    n_blue_steps = 100
    blue_colors = plt.cm.Blues(np.linspace(0.5, 0.0, n_blue_steps))  # Medium to very light
    blue_colors[-1] = [1, 1, 1, 1]  # Ensure last color is white
    blue_positions = np.linspace(0, threshold_position, n_blue_steps)
    
    for i in range(n_blue_steps):
        colors.append(blue_colors[i])
        positions.append(blue_positions[i])
    
    # Add explicit white at the threshold
    colors.append([1, 1, 1, 1])  # Pure white
    positions.append(threshold_position)
    
    # Significant part (threshold to max): White diverging to medium red
    n_red_steps = 100
    red_colors = plt.cm.Reds(np.linspace(0.0, 0.5, n_red_steps))  # Very light to medium
    red_colors[0] = [1, 1, 1, 1]  # Ensure first color is white
    red_positions = np.linspace(threshold_position, 1.0, n_red_steps)
    
    for i in range(n_red_steps):
        colors.append(red_colors[i])
        positions.append(red_positions[i])
    
    cmap_lower = mcolors.LinearSegmentedColormap.from_list('threshold_colormap', 
                                                          list(zip(positions, colors)))
    
    # Use normal normalization (not centered)
    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
    
    # Plot lower triangle (p-values)
    for i in range(n_categories):
        for j in range(i):
            p_val = pvalue_matrix[i, j]
            neg_log_p = -np.log10(max(p_val, 1e-300))
            
            # Color based on colormap
            color = cmap_lower(norm(neg_log_p))
            
            # Create rectangle for lower triangle cell (no border)
            rect = Rectangle((j - 0.5, i - 0.5), 1, 1, 
                           linewidth=0, edgecolor='none', 
                           facecolor=color, alpha=0.9)
            ax.add_patch(rect)
            
            # Add p-value text in exponential form
            p_text = f'{p_val:.2e}'
            
            # Use black text for all cells since colors are lighter
            text_color = 'black'
            
            ax.text(j, i, p_text, ha='center', va='center', 
                   fontsize=18, color=text_color, fontweight='normal')
    
    # Plot upper triangle (KS statistics) - symmetric colors to lower triangle
    for i in range(n_categories):
        for j in range(i+1, n_categories):
            ks_stat = ks_matrix[i, j]
            p_val = pvalue_matrix[i, j]
            neg_log_p = -np.log10(max(p_val, 1e-300))
            
            # Use the same color as the corresponding lower triangle cell
            color = cmap_lower(norm(neg_log_p))
            
            # Create rectangle for upper triangle cell (no border)
            rect = Rectangle((j - 0.5, i - 0.5), 1, 1, 
                           linewidth=0, edgecolor='none', 
                           facecolor=color, alpha=0.9)
            ax.add_patch(rect)
            
            # Add KS statistic value
            ks_text = f'{ks_stat:.3f}'
            ax.text(j, i, ks_text, ha='center', va='center', 
                   fontsize=18, color='black', fontweight='normal')
    
    # Plot diagonal (empty cells)
    for i in range(n_categories):
        rect = Rectangle((i - 0.5, i - 0.5), 1, 1, 
                       linewidth=0, edgecolor='none', facecolor='#E8E8E8')
        ax.add_patch(rect)
        # No text in diagonal cells
    
    # Add colorbar with threshold marked
    sm = plt.cm.ScalarMappable(cmap=cmap_lower, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=18) 
    cbar.set_label('P-value (-log10 scale)', rotation=270, labelpad=20, fontsize=18)
    
    # Add reference line at threshold in colorbar
    #cbar.ax.axhline(y=threshold_neg_log, color='black', linestyle='-', linewidth=1.5, alpha=0.8)
    #cbar.ax.text(1.8, threshold_neg_log, f'p={alpha}', rotation=0, va='center', ha='center', 
    #            fontsize=9, fontweight='bold', color='black')
    
    # Add significance labels on colorbar
    ylim = cbar.ax.get_ylim()
    #cbar.ax.text(2.5, (ylim[0] + threshold_neg_log) / 2, 'ns', 
    #            rotation=270, va='center', ha='center', fontsize=8, style='italic')
    #cbar.ax.text(2.5, (threshold_neg_log + ylim[1]) / 2, 'sig', 
    #            rotation=270, va='center', ha='center', fontsize=8, style='italic')
    
    # Add legend for heatmap content
    #legend_text = f'Lower triangle: p-values\nUpper triangle: KS statistics\n\nColor scale based on -log10(p-value)'
    legend_text = ""
    ax.text(1.25, 0.5, legend_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='center', bbox=dict(boxstyle='round', 
                                                facecolor='white', 
                                                alpha=0.8, 
                                                edgecolor='#cccccc'))
    
    # Customize axes
    ax.set_xticks(range(n_categories))
    ax.set_yticks(range(n_categories))
    ax.set_xticklabels(labels, fontsize=18)
    ax.set_yticklabels(labels, fontsize=18)
    ax.set_xlabel(cat_col_name, fontsize=18, labelpad=10)
    ax.set_ylabel(cat_col_name, fontsize=18, labelpad=10)
    ax.set_title(f'Pairwise KS Test Results\n{score_col_name} Distribution Comparison', 
                fontsize=18, fontweight='normal', pad=15)
    
    # Remove grid lines
    ax.grid(False)
    
    # Add border around plot area
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')
        spine.set_linewidth(0.8)
    
    # Set axis limits
    ax.set_xlim(-0.5, n_categories - 0.5)
    ax.set_ylim(n_categories - 0.5, -0.5)
    
    # Rotate x-axis labels if there are many categories
    if n_categories > 6:
        plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    
    base_filename = os.path.splitext(pdf_file)[0]

    saved_files = []
    
    # Save as PDF
    if save_pdf:
        plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"\n✓ PDF heatmap saved to: {pdf_file}")
        saved_files.append(pdf_file)
    
    # Save as PNG
    if save_png:
        png_file = f"{base_filename}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ PNG heatmap saved to: {png_file}")
        saved_files.append(png_file)
    
    #plt.show()

def main():
    """Main function to handle command line arguments"""

    # Optional --png flag: also write a PNG alongside the PDF. The manuscript
    # figures are vector PDFs, but the response-to-reviewers document embeds
    # raster copies, so both are needed from the same invocation.
    save_png = "--png" in sys.argv
    if save_png:
        sys.argv = [a for a in sys.argv if a != "--png"]

    
    # Check if arguments are provided
    if len(sys.argv) >= 4:
        filename = sys.argv[1]
        category_col = int(sys.argv[2])
        score_col = int(sys.argv[3])
        output_file = sys.argv[4] if len(sys.argv) > 4 else None
        alpha = float(sys.argv[5]) if len(sys.argv) > 5 else 0.05
    else:
        # Interactive mode
        print("Pairwise Kolmogorov-Smirnov Test with Combined Heatmap")
        print("-" * 60)
        
        # Get filename
        filename = input("Enter the input file path: ").strip()
        
        # Check if file exists
        if not os.path.exists(filename):
            print(f"Error: File '{filename}' does not exist.")
            return
        
        # Get column indices
        try:
            category_col = int(input("Enter the category column number (1-based): "))
            score_col = int(input("Enter the score column number (1-based): "))
            
            alpha_input = input("Enter significance threshold (default 0.05): ").strip()
            alpha = float(alpha_input) if alpha_input else 0.05
        except ValueError:
            print("Error: Invalid input.")
            return
    
    # Perform pairwise KS tests (saves both PDF and PNG by default)
    pairwise_ks_test(filename, category_col, score_col, output_file, alpha, save_heatmap=True, save_pdf=True, save_png=save_png)

if __name__ == "__main__":
    # Check if required packages are installed
    required_packages = ['pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn']
    
    import subprocess
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    
    main()
