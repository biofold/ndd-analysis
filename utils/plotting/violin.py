import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import os
from scipy import stats
import matplotlib.colors as mcolors

def create_violin_plot(filename, category_col, score_col, output_file=None, save_pdf=True, save_png=False, save_stats=True):
    """
    Create violin plots with lighter, pastel colors, white background, and plot border.
    
    Args:
        filename: Path to input file (CSV, Excel, or TSV)
        category_col: 1-based column index for categories
        score_col: 1-based column index for scores
        output_file: Output file name without extension (optional)
        save_pdf: Whether to save as PDF (default: True)
        save_png: Whether to save as PNG (default: True)
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
    print(f"Score column: '{score_col_name}' (column {score_col})")
    
    # Clean the data - remove NaN values
    df_clean = df[[cat_col_name, score_col_name]].dropna()
    
    # Convert columns to numeric
    df_clean[cat_col_name] = pd.to_numeric(df_clean[cat_col_name], errors='coerce')
    df_clean[score_col_name] = pd.to_numeric(df_clean[score_col_name], errors='coerce')
    
    # Remove any remaining NaN after conversion
    df_clean = df_clean.dropna()
    
    # Get unique categories and sort them
    categories = sorted(df_clean[cat_col_name].unique())
    print(f"\nCategories found (sorted): {categories}")
    print(f"Total data points: {len(df_clean)}")
    
    # Calculate statistics for each category
    stats_dict = {}
    for cat in categories:
        cat_data = df_clean[df_clean[cat_col_name] == cat][score_col_name]
        if len(cat_data) > 0:
            stats_dict[cat] = {
                'Category': cat,
                'Count': len(cat_data),
                'Mean': cat_data.mean(),
                'Median': cat_data.median(),
                'Std Dev': cat_data.std(),
                'Min': cat_data.min(),
                'Max': cat_data.max(),
                'Q1 (25th)': cat_data.quantile(0.25),
                'Q2 (50th)': cat_data.quantile(0.50),
                'Q3 (75th)': cat_data.quantile(0.75),
                'IQR': cat_data.quantile(0.75) - cat_data.quantile(0.25),
                'Skewness': cat_data.skew(),
                'Kurtosis': cat_data.kurtosis(),
                'SEM': cat_data.sem(),
                'Variance': cat_data.var()
            }
    
    # Create statistics DataFrame
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
    
    # Explicitly set plot area background to white
    ax.set_facecolor('white')
    
    # Create a single base color with very light, pastel variations
    # Using a soft blue as base
    base_color = '#87CEEB'  # Sky blue
    
    # Generate very light pastel colors
    n_categories = len(categories)
    
    # Create colors with increasing darkness but all very light
    if n_categories > 1:
        # Start with very light color and gradually darken slightly
        light_color = mcolors.to_rgb('#E8F4FD')  # Very light blue
        mid_color = mcolors.to_rgb('#B8DCF0')    # Light blue
        dark_color = mcolors.to_rgb('#8EC5E8')   # Slightly darker light blue
        
        colors = []
        for i in range(n_categories):
            if i < n_categories / 2:
                # Interpolate from very light to light
                t = i / max(n_categories / 2 - 1, 1)
                color = tuple(
                    light_c * (1 - t) + mid_c * t
                    for light_c, mid_c in zip(light_color, mid_color)
                )
            else:
                # Interpolate from light to slightly darker
                t = (i - n_categories / 2) / max(n_categories / 2 - 1, 1)
                color = tuple(
                    mid_c * (1 - t) + dark_c * t
                    for mid_c, dark_c in zip(mid_color, dark_color)
                )
            colors.append(color)
    else:
        colors = [mcolors.to_rgb('#B8DCF0')]
    
    # Create violin plot with custom colors
    parts = ax.violinplot(
        [df_clean[df_clean[cat_col_name] == cat][score_col_name].values for cat in categories],
        positions=range(len(categories)),
        showmeans=False,
        showmedians=True,
        showextrema=True,
        widths=0.8
    )
    
    # Customize violin colors with light pastels
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('#999999')  # Light gray edge
        pc.set_alpha(0.8)  # Slightly more opaque for better visibility
        pc.set_linewidth(0.8)  # Thinner line
    
    # Customize other parts
    parts['cmedians'].set_color('#333333')  # Darker for visibility
    parts['cmedians'].set_linewidth(1.0)  # Thinner line
    parts['cmins'].set_color('#666666')
    parts['cmaxes'].set_color('#666666')
    parts['cbars'].set_color('#666666')
    parts['cbars'].set_linewidth(0.8)  # Thinner line
    
    # Add individual data points with jitter (very light). Fixed seed: the jitter is
    # display only, and a seeded generator makes the figure byte-reproducible.
    rng = np.random.default_rng(0)
    for i, cat in enumerate(categories):
        cat_data = df_clean[df_clean[cat_col_name] == cat][score_col_name].values
        # Add jittered points with very light colors
        jitter = rng.normal(0, 0.04, size=len(cat_data))
        # Make points slightly darker than the violin but still light
        point_color = tuple(c * 0.7 for c in colors[i])
        ax.scatter(i + jitter, cat_data, alpha=0.3, s=8, 
                  color=point_color, edgecolors='none', zorder=3)
    
    # Customize x-axis
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([str(cat) for cat in categories])
    
    # Set labels with clean style
    ax.set_xlabel(cat_col_name, fontsize=13, fontweight='normal', labelpad=10)
    ax.set_ylabel(score_col_name, fontsize=13, fontweight='normal', labelpad=10)
    
    # Set title
    ax.set_title(f'Distribution of {score_col_name} by {cat_col_name}', 
                fontsize=15, fontweight='normal', pad=15)
    
    # Customize grid (light gray grid on y-axis only)
    ax.grid(True, axis='y', alpha=0.3, linestyle='-', linewidth=0.5, color='#CCCCCC')
    ax.grid(False, axis='x')
    ax.set_axisbelow(True)
    
    # Add border around the plot area (thinner lines)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('#333333')  # Dark gray border (not too harsh)
        spine.set_linewidth(0.8)    # Thin border line
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the plot in different formats
    if output_file is None:
        base_filename = f"violin_plot_{cat_col_name}_vs_{score_col_name}"
    else:
        # Remove extension if provided
        base_filename = os.path.splitext(output_file)[0]
    
    saved_files = []
    
    # Save as PDF
    if save_pdf:
        pdf_file = f"{base_filename}.pdf"
        plt.savefig(pdf_file, format='pdf', bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"\n✓ PDF saved as: {pdf_file}")
        saved_files.append(pdf_file)
    
    # Save as PNG
    if save_png:
        png_file = f"{base_filename}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"✓ PNG saved as: {png_file}")
        saved_files.append(png_file)
    
    # Save statistics to file
    if save_stats:
        # Save as CSV
        csv_stats_file = f"{output_file}"
        stats_df.to_csv(csv_stats_file, sep='\t', index=False)
        print(f"✓ Statistics CSV saved as: {csv_stats_file}")
        saved_files.append(csv_stats_file)
        
    
    # Print color information
    print("\nColor scheme (light pastels):")
    for i, (cat, color) in enumerate(zip(categories, colors)):
        hex_color = mcolors.to_hex(color)
        print(f"  Category {cat}: {hex_color}")
    
    # Show the plot
    #plt.show()
    
    # Print summary statistics to console
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(stats_df.to_string(index=False))
    
    return df_clean, stats_df

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
    else:
        # Interactive mode
        print("Violin Plot Generator (White Background with Clear Borders)")
        print("-" * 65)
        
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
        except ValueError:
            print("Error: Column numbers must be integers.")
            return
        
        # Optional output file
        output_file = input("Enter output file name (or press Enter for default): ").strip()
        if not output_file:
            output_file = None
    
    # Create the violin plot and save statistics
    create_violin_plot(filename, category_col, score_col, output_file, 
                      save_pdf=True, save_png=save_png, save_stats=True)

if __name__ == "__main__":
    # Check if required packages are installed
    required_packages = ['pandas', 'matplotlib', 'seaborn', 'numpy', 'scipy', 'openpyxl']
    
    import subprocess
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            print(f"Installing {package}...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    
    main()
