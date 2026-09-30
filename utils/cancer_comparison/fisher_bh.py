#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats

def akashi(mat):
    """Calculate Z-score using Akashi's method."""
    n = mat[0][0] + mat[0][1] + mat[1][0] + mat[1][1]
    p = (mat[0][0] + mat[0][1]) / float(n)
    q = (mat[0][0] + mat[1][0]) / float(n)
    E = n * p * q
    V = p * (1 - p) * q * (1 - q) * (n ** 2) / float(n - 1)
    Z = (mat[0][0] - E) / np.sqrt(V)
    return Z

def cal_fisher(ct):
    """Calculate Fisher exact test statistics."""
    Z = akashi(ct)
    oddsr, pvalueg = stats.fisher_exact(ct, alternative='greater')
    oddsr, pvaluel = stats.fisher_exact(ct, alternative='less')
    oddsr, pvalue2 = stats.fisher_exact(ct)
    
    # Determine the side (greater or less) based on smaller p-value
    if pvaluel < pvalueg:
        pvalue = pvaluel
        side = 'L'
    else:
        pvalue = pvalueg
        side = 'G'
    
    return {
        'odds_ratio': oddsr,
        'p_value_left': pvaluel,
        'p_value_right': pvalueg,
        'p_value': pvalue,
        'p_value_two_sided': pvalue2,
        'z_score': Z,
        'side': side
    }

def benjamini_hochberg(p_values):
    """Compute BH-adjusted p-values (no rejection logic)."""
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])
    
    # Sort p-values and compute ranks
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    ranks = np.arange(1, n + 1)
    
    # Calculate adjusted p-values
    adjusted_p = sorted_p * n / ranks
    
    # Ensure monotonicity
    for i in range(n - 2, -1, -1):
        adjusted_p[i] = min(adjusted_p[i], adjusted_p[i + 1])
    
    # Restore original order
    final_p = np.zeros_like(adjusted_p)
    final_p[sorted_indices] = adjusted_p
    return final_p

def process_fisher(filename, vpos1, vpos2, c=0):
    """Process file and calculate Fisher exact test for each line."""
    results = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            
            v = line.split()
            
            try:
                ct = np.array([
                    [float(v[vpos1[0]]) + c, float(v[vpos1[1]]) + c],
                    [float(v[vpos2[0]]) + c, float(v[vpos2[1]]) + c]
                ])
                
                fisher_results = cal_fisher(ct)
                
                # Store original line and calculated values
                results.append({
                    'original_line': line,
                    'original_fields': v,
                    'odds_ratio': fisher_results['odds_ratio'],
                    'p_value_left': fisher_results['p_value_left'],
                    'p_value_right': fisher_results['p_value_right'],
                    'p_value': fisher_results['p_value'],
                    'p_value_two_sided': fisher_results['p_value_two_sided'],
                    'z_score': fisher_results['z_score'],
                    'side': fisher_results['side']
                })
            except (ValueError, IndexError) as e:
                sys.stderr.write(f"Warning: Could not process line: {line[:50]}... Error: {e}\n")
                continue
    
    return results

def format_output(results, adjusted_p_values=None):
    """Format results for output."""
    output_lines = []
    
    for idx, result in enumerate(results):
        # Start with original line
        output_parts = [result['original_line']]
        
        # Add Fisher test results
        '''
        output_parts.extend([
            f"OR: {result['odds_ratio']:.3f}",
            f"p-valuel: {result['p_value_left']:.2e}",
            f"p-valueg: {result['p_value_right']:.2e}",
            f"p-value: {result['p_value']:.2e}",
            f"p-value2: {result['p_value_two_sided']:.2e}",
            f"Z: {result['z_score']:.2f}",
            result['side']
        ])
        '''
        # Full precision (shortest exact repr): calculate_table_compara.sh
        # thresholds these values, and rounding before a threshold can flip it.
        output_parts.extend([
            repr(float(result['odds_ratio'])),
            repr(float(result['p_value_left'])),
            repr(float(result['p_value_right'])),
            repr(float(result['p_value'])),
            repr(float(result['p_value_two_sided'])),
            repr(float(result['z_score'])),
            result['side']
        ])
        
        # Add BH corrected p-value if available
        if adjusted_p_values is not None:
            output_parts.append(repr(float(adjusted_p_values[idx])))
        
        output_lines.append('\t'.join(output_parts))
    
    return output_lines

def main():
    parser = argparse.ArgumentParser(
        description="Calculate Fisher exact test and apply Benjamini-Hochberg correction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: BH correction on p-valueg (right-sided, greater)
  %(prog)s input.txt 3,4 5,6
  
  # BH correction on two-sided p-value
  %(prog)s input.txt 3,4 5,6 --pvalue_type two_sided
  
  # BH correction on directional p-value (min of left/right)
  %(prog)s input.txt 3,4 5,6 --pvalue_type directional
  
  # BH correction on left-sided p-value
  %(prog)s input.txt 3,4 5,6 --pvalue_type left
  
  # Sort by odds ratio descending
  %(prog)s input.txt 3,4 5,6 --sort odds_ratio --sort_order desc
        """
    )
    parser.add_argument("input_file", help="Input file (space/tab-delimited)")
    parser.add_argument("pos1", help="Comma-separated column numbers for first row of contingency table (1-based)")
    parser.add_argument("pos2", help="Comma-separated column numbers for second row of contingency table (1-based)")
    parser.add_argument("--output", default="/dev/stdout", help="Output file (default: stdout)")
    parser.add_argument("--constant", type=float, default=0,
                       help="Constant to add to all cells (default: 0)")
    parser.add_argument("--no_bh", action="store_true",
                       help="Skip Benjamini-Hochberg correction")
    parser.add_argument("--pvalue_type", 
                       choices=['right', 'left', 'two_sided', 'directional'],
                       default='right',
                       help="Which p-value to use for BH correction (default: right = p-valueg)")
    parser.add_argument("--sort", type=str, default=None,
                       help="Sort key (default: BH-adjusted p-value if BH applied, else p-value)")
    parser.add_argument("--sort_order", choices=['asc', 'desc'], default='asc',
                       help="Sort order (default: asc)")
    parser.add_argument("--no_sort", action="store_true",
                       help="Do not sort results (preserve original order)")
    args = parser.parse_args()
    
    # Parse column positions
    vpos1 = [int(i) - 1 for i in args.pos1.split(',')]
    vpos2 = [int(i) - 1 for i in args.pos2.split(',')]
    
    # Validate positions
    if len(vpos1) != 2 or len(vpos2) != 2:
        sys.stderr.write("Error: Each position must have exactly 2 columns\n")
        sys.exit(1)
    
    # Process Fisher test
    sys.stderr.write("Calculating Fisher exact test...\n")
    results = process_fisher(args.input_file, vpos1, vpos2, args.constant)
    
    if not results:
        sys.stderr.write("Error: No valid results found\n")
        sys.exit(1)
    
    # Select p-values for BH correction based on user preference
    if args.pvalue_type == 'right':
        p_values = [r['p_value_right'] for r in results]
        pvalue_label = "right-sided (p-valueg)"
    elif args.pvalue_type == 'left':
        p_values = [r['p_value_left'] for r in results]
        pvalue_label = "left-sided (p-valuel)"
    elif args.pvalue_type == 'two_sided':
        p_values = [r['p_value_two_sided'] for r in results]
        pvalue_label = "two-sided (p-value2)"
    elif args.pvalue_type == 'directional':
        p_values = [r['p_value'] for r in results]
        pvalue_label = "directional (p-value)"
    
    # Apply BH correction
    if not args.no_bh:
        sys.stderr.write(f"Applying Benjamini-Hochberg correction using {pvalue_label} p-values...\n")
        adjusted_p_values = benjamini_hochberg(np.array(p_values))
        # Add adjusted p-values to results
        for idx, result in enumerate(results):
            result['adjusted_p_value'] = adjusted_p_values[idx]
    else:
        adjusted_p_values = None
    
    # Sort results (default: by BH-adjusted p-value ascending)
    if not args.no_sort:
        # Determine sort key
        if args.sort:
            sort_key = args.sort
        else:
            # Default: sort by BH-adjusted p-value if available, else p-value
            if not args.no_bh:
                sort_key = 'adjusted_p'
            else:
                sort_key = 'p_value_right'
        
        sys.stderr.write(f"Sorting by: {sort_key} ({args.sort_order})...\n")
        
        if sort_key.startswith('column:'):
            # Sort by original column
            col_num = int(sort_key.split(':')[1]) - 1  # Convert to 0-based
            results.sort(key=lambda x: float(x['original_fields'][col_num]), 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'adjusted_p' and not args.no_bh:
            results.sort(key=lambda x: x['adjusted_p_value'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'p_value':
            results.sort(key=lambda x: x['p_value'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'p_value_right':
            results.sort(key=lambda x: x['p_value_right'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'p_value_left':
            results.sort(key=lambda x: x['p_value_left'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'p_value_two_sided':
            results.sort(key=lambda x: x['p_value_two_sided'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'odds_ratio':
            results.sort(key=lambda x: x['odds_ratio'], 
                        reverse=(args.sort_order == 'desc'))
        elif sort_key == 'z_score':
            results.sort(key=lambda x: x['z_score'], 
                        reverse=(args.sort_order == 'desc'))
        else:
            sys.stderr.write(f"Warning: Unknown sort key: {sort_key}. Skipping sorting.\n")
    
    # Generate output
    output_lines = format_output(results, adjusted_p_values)
    
    # Write output
    if args.output == "/dev/stdout":
        for line in output_lines:
            print(line)
    else:
        with open(args.output, 'w') as f:
            for line in output_lines:
                f.write(line + '\n')
        sys.stderr.write(f"Results saved to: {args.output}\n")
    
    # Print summary
    sys.stderr.write(f"Processed {len(results)} lines\n")
    if not args.no_bh:
        sys.stderr.write(f"BH correction applied to {pvalue_label} p-values\n")
    if not args.no_sort:
        sys.stderr.write(f"Results sorted by {sort_key}\n")

if __name__ == "__main__":
    main()
