import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable

def create_scatter_plot(file_path, x_col, y_col, color_col, output_file=None):
    """
    Create a square scatter plot with:
    - Color scale from 0 to 30 (white at 2)
    - x=2 and y=2 reference lines
    - Point contours
    - Numbered axis labels
    - Axis limits extended by -1/+1
    - Fixed bin size of 5 for ticks
    - Colorbar same height as plot
    """
    
    # Read the file
    try:
        df = pd.read_csv(file_path, sep=None, engine='python', header=None)
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Check columns exist
    if max(x_col, y_col, color_col) >= len(df.columns):
        print("Error: One or more specified columns don't exist in the file")
        return
    
    # Extract data
    x = df.iloc[:, x_col]
    y = df.iloc[:, y_col]
    colors = df.iloc[:, color_col]
    
    # Create custom red-white-blue colormap centered at 2
    cmap = LinearSegmentedColormap.from_list('rwb', [
        (0.0, 'red'),          # 0 = red
        (2/30, 'white'),       # 2 = white
        (1.0, 'blue')          # 30 = blue
    ])
    
    # Set color scale from 0 to 30
    norm = plt.Normalize(vmin=0, vmax=30)
    
    # Create square figure with adjusted width to accommodate colorbar
    fig, ax = plt.subplots(figsize=(8.5, 8))
    
    # Scatter plot with contours
    sc = ax.scatter(
        x, y, 
        c=colors, 
        cmap=cmap, 
        norm=norm,
        alpha=0.7,
        edgecolors='lightgrey',
        linewidths=0.5,
        s=60
    )
    
    # Add reference lines
    ax.axvline(x=2, color='gray', linestyle='--', linewidth=1)
    ax.axhline(y=2, color='gray', linestyle='--', linewidth=1)
    
    # Create axes divider with increased padding
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    
    # Add colorbar with same height as plot
    cbar = plt.colorbar(sc, cax=cax)
    cbar.set_ticks(np.linspace(0, 30, 7))
    cbar.ax.tick_params(labelsize=12, width=1.5)  # Increased colorbar tick label size
    #for label in cbar.ax.get_yticklabels():
    #    label.set_fontweight('bold')  # Make text bold
    
    # Set equal aspect ratio for square plot
    ax.set_aspect('equal', 'box')
    
    # Calculate axis limits with -1/+1 padding
    x_min, x_max = -1, 31
    y_min, y_max = -1, 31
    
    # Set axis limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # Set fixed bin size of 5 for ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
   
    # Set larger font size for axes tick labels
    ax.tick_params(axis='both', which='major', labelsize=12, width=1.5)
    #for label in ax.get_xticklabels() + ax.get_yticklabels():
    #    label.set_fontweight('bold')
 
    # Set axis labels
    #ax.set_xlabel('X', fontsize=14)
    #ax.set_ylabel('Y', fontsize=14)
    
    # Remove title
    ax.set_title('')
    
    # Adjust layout to account for colorbar padding
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
    else:
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create scatter plot with extended axes and fixed tick spacing')
    parser.add_argument('file_path', help='Path to data file')
    parser.add_argument('x_col', type=int, help='X-axis column (1-based)')
    parser.add_argument('y_col', type=int, help='Y-axis column (1-based)')
    parser.add_argument('color_col', type=int, help='Color column (1-based)')
    parser.add_argument('--output', '-o', help='Output file to save plot (optional)', default=None)
    
    args = parser.parse_args()
    create_scatter_plot(args.file_path, args.x_col-1, args.y_col-1, args.color_col-1, args.output)
