#!/usr/bin/env python3
"""
Convert multiple TSV files to an Excel file with each TSV as a separate sheet.
Sheet names are derived from the TSV filenames (without extension).

Usage: python tsv_to_excel.py output.xlsx file1.tsv file2.tsv file3.tsv ...
"""

import sys
import pandas as pd
from pathlib import Path


def tsv_to_excel(output_file, input_files):
    """
    Convert multiple TSV files to an Excel file with separate sheets.
    
    Args:
        output_file (str): Path to the output Excel file
        input_files (list): List of paths to TSV files
    """
    # Check if there are input files
    if not input_files:
        print("Error: No input TSV files provided.")
        return False
    
    try:
        # Create Excel writer
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for tsv_file in input_files:
                file_path = Path(tsv_file)
                
                # Check if file exists
                if not file_path.exists():
                    print(f"Warning: File '{tsv_file}' does not exist. Skipping...")
                    continue
                
                # Get sheet name from filename (without extension)
                # Excel sheet names must be 31 characters or less
                sheet_name = file_path.stem #[:31]
                
                print(f"Processing: {tsv_file} -> Sheet: '{sheet_name}'")
                
                # Read TSV file
                df = pd.read_csv(tsv_file, sep='\t')

                # Remove the first column
                df = df.iloc[:, 1:]  # Keep all rows, all columns except the first
                
                # Write to Excel sheet
                sheet_name=sheet_name[10:]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"\nSuccess! Excel file created: {output_file}")
        print(f"Total sheets created: {len([f for f in input_files if Path(f).exists()])}")
        return True
        
    except Exception as e:
        print(f"Error creating Excel file: {e}")
        return False


def main():
    # Check command line arguments
    if len(sys.argv) < 3:
        print("Usage: python tsv_to_excel.py <output.xlsx> <file1.tsv> [file2.tsv ...]")
        print("\nExample: python tsv_to_excel.py combined.xlsx data1.tsv data2.tsv data3.tsv")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    # Validate output file extension
    if not output_file.endswith('.xlsx'):
        output_file += '.xlsx'
        print(f"Note: Adding .xlsx extension to output file: {output_file}")
    
    # Convert TSV files to Excel
    success = tsv_to_excel(output_file, input_files)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
