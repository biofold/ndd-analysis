#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pandas as pd

def benjamini_hochberg(p_values):
    """Compute BH-adjusted p-values (no rejection logic)."""
    p_values = np.asarray(p_values)
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

def main():
    parser = argparse.ArgumentParser(
        description="Apply Benjamini-Hochberg correction to p-values in a file."
    )
    parser.add_argument("input_file", help="Input file (CSV/TSV/space-delimited)")
    parser.add_argument("column", type=int, help="Column number (1-based) containing p-values")
    parser.add_argument("--output", default="/dev/stdout", help="Output file (default: stdout)")
    args = parser.parse_args()

    try:
        # Read input file
        df = pd.read_csv(args.input_file, sep=None, engine='python', header=None)
        p_values = df.iloc[:, args.column-1].astype(float)
        
        # Compute adjusted p-values
        df["Adjusted_p"] = benjamini_hochberg(p_values)
        
        # Write output
        df.to_csv(args.output, sep='\t', index=False, header=False)
        if args.output != "/dev/stdout":
            print(f"Adjusted p-values saved to: {args.output}", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
