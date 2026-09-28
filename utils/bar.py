import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import matplotlib.colors as mcolors

def create_bar_plot(filename, category_col, score_col, output_file=None, save_pdf=True, save_png=False, save_stats=True):
    """
    Create a proportion bar plot (with 95% Wilson confidence intervals) for a
    binary (0/1) outcome broken down by category, using the same light
    pastel-blue palette, white background, and plot border convention as
    violin.py.

    Args:
        filename: Path to input file (CSV, Excel, or TSV)
        category_col: 1-based column index for categories
        score_col: 1-based column index for the binary (0/1) outcome
        output_file: Output file name without extension (optional)
        save_pdf: Whether to save as PDF (default: True)
        save_png: Whether to save as PNG (default: False)
        save_stats: Whether to save statistics to file (default: True)
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
            # Try to guess the separator
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
    print(f"Score column (binary 0/1): '{score_col_name}' (column {score_col})")

    # Clean the data - remove NaN values
    df_clean = df[[cat_col_name, score_col_name]].dropna()

    # Convert columns to numeric
    df_clean[cat_col_name] = pd.to_numeric(df_clean[cat_col_name], errors='coerce')
    df_clean[score_col_name] = pd.to_numeric(df_clean[score_col_name], errors='coerce')

    # Remove any remaining NaN after conversion
    df_clean = df_clean.dropna()

    if not set(df_clean[score_col_name].unique()) <= {0, 1}:
        print("Error: Score column must be binary (0/1).")
        return None

    # Get unique categories and sort them
    categories = sorted(df_clean[cat_col_name].unique())
    print(f"\nCategories found (sorted): {categories}")
    print(f"Total data points: {len(df_clean)}")

    # Calculate proportion and Wilson 95% CI for each category
    z = 1.96
    stats_dict = {}
    for cat in categories:
        cat_data = df_clean[df_clean[cat_col_name] == cat][score_col_name]
        n = len(cat_data)
        if n == 0:
            continue
        k = cat_data.sum()
        p = k / n
        center = (p + z**2 / (2 * n)) / (1 + z**2 / n)
        halfwidth = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / (1 + z**2 / n)
        stats_dict[cat] = {
            'Category': cat,
            'Count': n,
            'Positive': int(k),
            'Proportion': p,
            'CI_lower': center - halfwidth,
            'CI_upper': center + halfwidth,
        }

    stats_df = pd.DataFrame.from_dict(stats_dict, orient='index')
    stats_df = stats_df.reset_index(drop=True)

    # Reset to default style and clear any previous settings
    plt.rcParams.update(plt.rcParamsDefault)

    # Set font to be PDF-friendly
    plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

    # Create figure with white background
    fig, ax = plt.subplots(figsize=(10, 7), facecolor='white')
    ax.set_facecolor('white')

    # Same light pastel-blue single color used by violin.py's mid tone
    bar_color = '#B8DCF0'

    x = np.arange(len(categories))
    props = stats_df['Proportion'].values
    ci_lo = stats_df['CI_lower'].values
    ci_hi = stats_df['CI_upper'].values

    ax.bar(x, props, width=0.6, color=bar_color, edgecolor='#999999', linewidth=0.8,
           yerr=[props - ci_lo, ci_hi - props],
           error_kw=dict(ecolor='#333333', elinewidth=1.2, capsize=4))

    # Customize x-axis
    ax.set_xticks(x)
    ax.set_xticklabels([str(cat) for cat in categories])

    # Set labels with clean style
    ax.set_xlabel(cat_col_name, fontsize=13, fontweight='normal', labelpad=10)
    ax.set_ylabel(f'Proportion of {score_col_name} = 1', fontsize=13, fontweight='normal', labelpad=10)

    # Set title
    ax.set_title(f'Proportion of {score_col_name} by {cat_col_name}',
                 fontsize=15, fontweight='normal', pad=15)

    # Customize grid (light gray grid on y-axis only)
    ax.grid(True, axis='y', alpha=0.3, linestyle='-', linewidth=0.5, color='#CCCCCC')
    ax.grid(False, axis='x')
    ax.set_axisbelow(True)

    # Add border around the plot area (thinner lines) -- all four spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')
        spine.set_linewidth(0.8)

    plt.tight_layout()

    # Save the plot in different formats
    if output_file is None:
        base_filename = f"bar_plot_{cat_col_name}_vs_{score_col_name}"
    else:
        base_filename = os.path.splitext(output_file)[0]

    saved_files = []

    if save_pdf:
        pdf_file = f"{base_filename}.pdf"
        plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"\n✓ PDF saved as: {pdf_file}")
        saved_files.append(pdf_file)

    if save_png:
        png_file = f"{base_filename}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ PNG saved as: {png_file}")
        saved_files.append(png_file)

    if save_stats:
        csv_stats_file = f"{output_file}"
        stats_df.to_csv(csv_stats_file, sep='\t', index=False)
        print(f"✓ Statistics CSV saved as: {csv_stats_file}")
        saved_files.append(csv_stats_file)

    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(stats_df.to_string(index=False))

    return df_clean, stats_df

def main():
    """Main function to handle command line arguments"""

    if len(sys.argv) >= 4:
        filename = sys.argv[1]
        category_col = int(sys.argv[2])
        score_col = int(sys.argv[3])
        output_file = sys.argv[4] if len(sys.argv) > 4 else None
    else:
        print("Bar Plot Generator for Binary Outcomes (White Background with Clear Borders)")
        print("-" * 65)

        filename = input("Enter the input file path: ").strip()

        if not os.path.exists(filename):
            print(f"Error: File '{filename}' does not exist.")
            return

        try:
            category_col = int(input("Enter the category column number (1-based): "))
            score_col = int(input("Enter the binary score column number (1-based): "))
        except ValueError:
            print("Error: Column numbers must be integers.")
            return

        output_file = input("Enter output file name (or press Enter for default): ").strip()
        if not output_file:
            output_file = None

    create_bar_plot(filename, category_col, score_col, output_file,
                     save_pdf=True, save_png=False, save_stats=True)

if __name__ == "__main__":
    required_packages = ['pandas', 'matplotlib', 'numpy']

    import subprocess

    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

    main()
