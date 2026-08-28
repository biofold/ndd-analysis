#!/usr/bin/env python3
"""
Convert multiple delimited files to an Excel file with each file as a separate sheet.
Sheet names can be derived from the filenames (without extension) or specified explicitly.

Usage:
  python tsv2excel.py output.xlsx file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --header file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --separator ',' file1.csv file2.csv
  python tsv2excel.py output.xlsx --quotechar '"' file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --header --skip_first_column file1.tsv file2.tsv
"""

import sys
import argparse
import pandas as pd
from pathlib import Path


def parse_input_spec(input_spec):
    """
    Parse input specification to extract file path and optional sheet name.
    """
    if ':' in input_spec:
        parts = input_spec.rsplit(':', 1)
        file_path = parts[0]
        sheet_name = parts[1].strip() if len(parts) > 1 else None
    else:
        file_path = input_spec
        sheet_name = None
    return file_path, sheet_name


def sanitize_sheet_name(name, max_length=31):
    """
    Sanitize sheet name for Excel compatibility.
    """
    if not name:
        return "Sheet"
    invalid_chars = ['\\', '/', '?', '*', '[', ']', ':']
    for char in invalid_chars:
        name = name.replace(char, '_')
    name = name.strip("'")
    if len(name) > max_length:
        name = name[:max_length]
    if not name:
        name = "Sheet"
    return name


def tsv_to_excel(output_file, input_specs, header=False, separator='\t',
                 quotechar='"', skip_first_column=False):
    """
    Convert multiple delimited files to an Excel file with separate sheets.

    Args:
        output_file (str): Path to the output Excel file
        input_specs (list): List of input specifications
        header (bool): If True, treat the first row of each file as column names.
                       Default is False (no header row in the output).
        separator (str): Field delimiter (default: tab).
        quotechar (str): Character used to quote fields (default: ").
        skip_first_column (bool): If True, remove the first column from each file.
    """
    if not input_specs:
        print("Error: No input files provided.")
        return False

    try:
        used_sheet_names = set()
        processed_count = 0

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for input_spec in input_specs:
                file_path_str, custom_sheet_name = parse_input_spec(input_spec)
                file_path = Path(file_path_str)

                if not file_path.exists():
                    print(f"Warning: File '{file_path_str}' does not exist. Skipping...")
                    continue
                if file_path.stat().st_size == 0:
                    print(f"Warning: File '{file_path_str}' is empty. Skipping...")
                    continue

                sheet_name = custom_sheet_name if custom_sheet_name else file_path.stem
                sheet_name = sanitize_sheet_name(sheet_name)

                # Ensure unique sheet name
                if sheet_name in used_sheet_names:
                    base = sheet_name
                    counter = 2
                    while f"{base}_{counter}" in used_sheet_names:
                        counter += 1
                    sheet_name = f"{base}_{counter}"
                used_sheet_names.add(sheet_name)

                # Read file with or without header
                try:
                    if header:
                        # First row as column names
                        df = pd.read_csv(file_path, sep=separator, quotechar=quotechar,
                                         on_bad_lines='skip', engine='python')
                    else:
                        # No header: all rows are data
                        df = pd.read_csv(file_path, sep=separator, quotechar=quotechar,
                                         header=None, on_bad_lines='skip', engine='python')
                except Exception as e:
                    print(f"Warning: Could not read '{file_path_str}': {e}. Skipping...")
                    continue

                if df.empty:
                    print(f"Warning: '{file_path_str}' contains no data. Skipping...")
                    continue

                # Optionally skip first column
                if skip_first_column and len(df.columns) > 1:
                    df = df.iloc[:, 1:]

                # Write to Excel; include header row only if header=True
                df.to_excel(writer, sheet_name=sheet_name, index=False, header=header)
                processed_count += 1
                print(f"  ✓ Added sheet: '{sheet_name}' ({len(df)} rows)")

        if processed_count == 0:
            print("Error: No valid files to process.")
            return False

        print(f"\nSuccess! Excel file created: {output_file}")
        print(f"Total sheets created: {processed_count}")
        return True

    except Exception as e:
        print(f"Error creating Excel file: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Convert multiple delimited files to an Excel file with separate sheets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # No header (default)
  python tsv2excel.py output.xlsx data1.tsv data2.tsv

  # First row is header
  python tsv2excel.py output.xlsx --header data1.tsv data2.tsv

  # Custom separator (comma)
  python tsv2excel.py output.xlsx --separator ',' data1.csv data2.csv

  # Remove first column
  python tsv2excel.py output.xlsx --skip_first_column data1.tsv data2.tsv

  # Custom sheet names
  python tsv2excel.py output.xlsx data1.tsv:Table1 data2.tsv:Table2
        """
    )

    parser.add_argument("output_file", help="Path to the output Excel file (.xlsx)")
    parser.add_argument("input_files", nargs='+',
                        help="Input files (format: file.tsv or file.tsv:sheet_name)")
    parser.add_argument("--header", action="store_true",
                        help="Treat first row as header (default: no header)")
    parser.add_argument("--separator", "--sep", dest="separator", default="\t",
                        help="Field delimiter (default: \\t)")
    parser.add_argument("--quotechar", dest="quotechar", default='"',
                        help="Character used to quote fields (default: \")")
    parser.add_argument("--skip_first_column", "-s", action="store_true",
                        help="Remove the first column from each input file")

    args = parser.parse_args()

    output_file = args.output_file
    if not output_file.endswith('.xlsx'):
        output_file += '.xlsx'
        print(f"Note: Adding .xlsx extension to output file: {output_file}")

    success = tsv_to_excel(
        output_file,
        args.input_files,
        header=args.header,
        separator=args.separator,
        quotechar=args.quotechar,
        skip_first_column=args.skip_first_column
    )

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
