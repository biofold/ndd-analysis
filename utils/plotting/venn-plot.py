#!/usr/bin/env python3

import sys
import os
from collections import defaultdict

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input.tsv> <class_name> [subgroups]")
        print()
        print("Arguments:")
        print("  input.tsv    - File with 4 columns: gene, class, score, pipe-separated subgroups")
        print("  class_name   - Class to filter (e.g., curated, candidate, no_evidence)")
        print("  subgroups    - Comma-separated subgroups to compare (e.g., hc,sanchis,orpha)")
        print("                 Note: 'orpha' is a special subgroup = union of 'neuro' and 'develop'")
        print("                 If not provided, uses all subgroups found in the class")
        print()
        print("Note: Lines starting with # are treated as comments and skipped")
        sys.exit(1)
    
    input_file = sys.argv[1]
    class_name = sys.argv[2]
    subgroups_arg = sys.argv[3] if len(sys.argv) > 3 else None
    output_file = sys.argv[4] if len(sys.argv) > 4 else f"venn_{class_name}_{'_'.join(subgroups_to_plot)}.png"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: File not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Parse subgroups if provided
    specified_subgroups = []
    if subgroups_arg:
        specified_subgroups = [s.strip().lower() for s in subgroups_arg.split(',') if s.strip()]
    
    # Read gene list and organize by subgroup
    gene_subgroups = defaultdict(set)  # subgroup -> set of genes
    class_genes = set()  # all genes in the specified class
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Skip comments (lines starting with #) and empty lines
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 4:
                print(f"Warning: Line {line_num} has fewer than 4 columns, skipping", file=sys.stderr)
                continue
            
            gene = fields[0].strip()
            gene_class = fields[1].strip()
            score = fields[2].strip()
            evidence = fields[3].strip()
            
            # Skip if gene is empty
            if not gene:
                continue
            
            # Filter by class
            if gene_class != class_name:
                continue
            
            class_genes.add(gene)
            
            # Parse subgroups if evidence is not 'none'
            if evidence.lower() != 'none' and evidence:
                subgroups = [s.strip().lower() for s in evidence.split('|') if s.strip()]
                for subgroup in subgroups:
                    gene_subgroups[subgroup].add(gene)
    
    # Check if we found any genes for the specified class
    if not class_genes:
        print(f"Error: No genes found for class '{class_name}'", file=sys.stderr)
        sys.exit(1)
    
    # Handle "orpha" subgroup = union of "neuro" and "develop"
    if 'neuro' in gene_subgroups or 'develop' in gene_subgroups:
        orpha_genes = set()
        if 'neuro' in gene_subgroups:
            orpha_genes.update(gene_subgroups['neuro'])
        if 'develop' in gene_subgroups:
            orpha_genes.update(gene_subgroups['develop'])
        if orpha_genes:
            gene_subgroups['orpha'] = orpha_genes
    
    # Determine which subgroups to use
    if specified_subgroups:
        # Use specified subgroups
        subgroups_to_plot = []
        missing = []
        
        for subgroup in specified_subgroups:
            if subgroup == 'orpha':
                # Special handling for orpha
                if 'orpha' in gene_subgroups:
                    subgroups_to_plot.append('orpha')
                else:
                    missing.append('orpha (requires neuro and/or develop)')
            elif subgroup in gene_subgroups:
                subgroups_to_plot.append(subgroup)
            else:
                missing.append(subgroup)
        
        if missing:
            print(f"Warning: Subgroups not found: {', '.join(missing)}", file=sys.stderr)
    else:
        # Use all subgroups found
        subgroups_to_plot = sorted(gene_subgroups.keys())
    
    if len(subgroups_to_plot) < 2:
        print(f"Error: Need at least 2 subgroups to create a Venn diagram", file=sys.stderr)
        print(f"Available subgroups: {', '.join(sorted(gene_subgroups.keys()))}", file=sys.stderr)
        sys.exit(1)
    
    if len(subgroups_to_plot) > 3:
        print(f"Warning: Venn diagram supports max 3 subgroups. Using first 3: {', '.join(subgroups_to_plot[:3])}", file=sys.stderr)
        subgroups_to_plot = subgroups_to_plot[:3]
    
    # Try to import matplotlib
    try:
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        import matplotlib.pyplot as plt
        from matplotlib_venn import venn2, venn3
    except ImportError as e:
        print(f"Error: Required libraries not found: {e}", file=sys.stderr)
        print("Please install: pip install matplotlib matplotlib-venn", file=sys.stderr)
        sys.exit(1)
    
    # Set ggplot2 style
    plt.style.use('ggplot')
    
    # Define light, pastel colors similar to ggplot2
    # These are lighter versions of the default ggplot2 palette
    light_colors_2 = ['#FBB4AE', '#B3CDE3']  # Light red/pink and light blue
    light_colors_3 = ['#FBB4AE', '#CCEBC5', '#B3CDE3']  # Light red/pink, light green, light blue
    
    # Alternative: even lighter pastel colors
    # light_colors_2 = ['#FDDBC7', '#D1E5F0']  # Very light orange and very light blue
    # light_colors_3 = ['#FDDBC7', '#D9F0D3', '#D1E5F0']  # Very light pastels
    
    # Create Venn diagram
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if len(subgroups_to_plot) == 2:
        # 2-way Venn diagram
        set1 = gene_subgroups[subgroups_to_plot[0]]
        set2 = gene_subgroups[subgroups_to_plot[1]]
        
        venn = venn2([set1, set2], set_labels=subgroups_to_plot, ax=ax)
        
        # Apply light colors
        if venn:
            for i, patch in enumerate(venn.patches):
                if patch:
                    patch.set_facecolor(light_colors_2[i % len(light_colors_2)])
                    patch.set_edgecolor('#333333')  # Dark gray edges
                    patch.set_linewidth(1.5)
                    patch.set_alpha(0.6)  # Semi-transparent
        
        plt.title(f"Venn Diagram - {class_name.capitalize()} Genes\n{subgroups_to_plot[0]} vs {subgroups_to_plot[1]}", 
                  fontsize=14, fontweight='bold')
        
    elif len(subgroups_to_plot) == 3:
        # 3-way Venn diagram
        set1 = gene_subgroups[subgroups_to_plot[0]]
        set2 = gene_subgroups[subgroups_to_plot[1]]
        set3 = gene_subgroups[subgroups_to_plot[2]]
        
        venn = venn3([set1, set2, set3], set_labels=subgroups_to_plot, ax=ax)
        
        # Apply light colors
        if venn:
            for i, patch in enumerate(venn.patches):
                if patch:
                    patch.set_facecolor(light_colors_3[i % len(light_colors_3)])
                    patch.set_edgecolor('#333333')  # Dark gray edges
                    patch.set_linewidth(1.5)
                    patch.set_alpha(0.6)  # Semi-transparent
        
        plt.title(f"Venn Diagram - {class_name.capitalize()} Genes\n{subgroups_to_plot[0]} vs {subgroups_to_plot[1]} vs {subgroups_to_plot[2]}", 
                  fontsize=14, fontweight='bold')
    
    # Add total genes information
    total_genes = len(class_genes)
    plt.text(0.5, -0.1, f"Total {class_name} genes: {total_genes:,}", 
             horizontalalignment='center', verticalalignment='center', 
             transform=ax.transAxes, fontsize=10, style='italic')
    
    # Save the figure
    # output_file = f"venn_{class_name}_{'_'.join(subgroups_to_plot)}.png"
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Venn diagram saved to: {output_file}")
    
    # Print statistics
    print(f"\nStatistics for {class_name} genes:")
    print(f"Total genes in class: {total_genes:,}")
    print()
    
    for subgroup in subgroups_to_plot:
        genes = gene_subgroups[subgroup]
        pct = (len(genes) / total_genes * 100) if total_genes > 0 else 0
        if subgroup == 'orpha':
            print(f"  {subgroup} (union of neuro + develop): {len(genes):,} genes ({pct:.1f}%)")
        else:
            print(f"  {subgroup}: {len(genes):,} genes ({pct:.1f}%)")
    
    # Print overlap information
    if len(subgroups_to_plot) == 2:
        set1 = gene_subgroups[subgroups_to_plot[0]]
        set2 = gene_subgroups[subgroups_to_plot[1]]
        overlap = set1 & set2
        print(f"\n  Overlap ({subgroups_to_plot[0]} ∩ {subgroups_to_plot[1]}): {len(overlap):,} genes")
    elif len(subgroups_to_plot) == 3:
        set1 = gene_subgroups[subgroups_to_plot[0]]
        set2 = gene_subgroups[subgroups_to_plot[1]]
        set3 = gene_subgroups[subgroups_to_plot[2]]
        overlap_all = set1 & set2 & set3
        overlap_12 = (set1 & set2) - set3
        overlap_13 = (set1 & set3) - set2
        overlap_23 = (set2 & set3) - set1
        print(f"\n  Overlap all three: {len(overlap_all):,} genes")
        print(f"  Overlap {subgroups_to_plot[0]} ∩ {subgroups_to_plot[1]} (only): {len(overlap_12):,} genes")
        print(f"  Overlap {subgroups_to_plot[0]} ∩ {subgroups_to_plot[2]} (only): {len(overlap_13):,} genes")
        print(f"  Overlap {subgroups_to_plot[1]} ∩ {subgroups_to_plot[2]} (only): {len(overlap_23):,} genes")

if __name__ == "__main__":
    main()
