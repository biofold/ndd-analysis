#!/usr/bin/env python3
"""
Convert multiple delimited files to an Excel file with each file as a separate sheet.
Sheet names can be derived from the filenames (without extension) or specified explicitly.

Input format: file.tsv[:sheet_name[:header_rows]]

Where:
  file.tsv     - Path to input file
  sheet_name   - Optional: Name for the Excel sheet (default: filename without extension)
  header_rows  - Optional: Number of header rows
                 0 = no header (default)
                 N = number of rows to treat as header (e.g., 1 = first row is header)

Usage:
  python tsv2excel.py output.xlsx file1.tsv
  python tsv2excel.py output.xlsx file1.tsv:MySheet
  python tsv2excel.py output.xlsx file1.tsv:MySheet:1
  python tsv2excel.py output.xlsx file1.tsv:MySheet:0
  python tsv2excel.py output.xlsx --separator tab file1.tsv:Sheet1:1
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from openpyxl.styles import Font, PatternFill, Alignment
import re


def parse_input_spec(input_spec):
    """
    Parse input specification to extract file path, sheet name, and header rows.
    Format: file.tsv[:sheet_name[:header_rows]]
    
    header_rows:
      0 = no header (default)
      N = number of rows to treat as header (e.g., 1 = first row is header)
    """
    parts = input_spec.split(':')
    file_path = parts[0]
    
    sheet_name = None
    header_rows = 0  # Default: no header
    
    if len(parts) > 1 and parts[1].strip():
        sheet_name = parts[1].strip()
    
    if len(parts) > 2 and parts[2].strip():
        try:
            header_rows = int(parts[2].strip())
        except ValueError:
            print(f"Warning: Invalid header_rows value '{parts[2]}'. Using 0 (no header).", file=sys.stderr)
            header_rows = 0
    
    return file_path, sheet_name, header_rows


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


def get_separator(separator_arg):
    """
    Convert separator argument to actual separator character.
    """
    if separator_arg is None:
        return None
    
    separator_map = {
        'tab': '\t',
        '\\t': '\t',
        '\t': '\t',
        'comma': ',',
        ',': ',',
        'semicolon': ';',
        ';': ';',
        'space': ' ',
        ' ': ' ',
        'pipe': '|',
        '|': '|'
    }
    
    # Try direct mapping
    if separator_arg in separator_map:
        return separator_map[separator_arg]
    
    # Handle escaped characters
    if separator_arg.startswith('\\') and len(separator_arg) == 2:
        escape_map = {
            '\\t': '\t',
            '\\n': '\n',
            '\\r': '\r',
            '\\s': ' '
        }
        return escape_map.get(separator_arg, separator_arg)
    
    # If it's a single character, use it directly
    if len(separator_arg) == 1:
        return separator_arg
    
    return separator_arg


def read_file_with_separator(file_path, separator):
    """
    Read file with a specific separator.
    """
    lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.strip() or separator in line:
                    lines.append(line)
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin1') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    
    if not lines:
        return pd.DataFrame()
    
    # Parse lines using the specified separator
    parsed_lines = []
    for line in lines:
        if separator == ' ':
            fields = [f for f in line.split(' ') if f != '']
        else:
            fields = line.split(separator)
        parsed_lines.append(fields)
    
    # Find max columns
    max_cols = max(len(fields) for fields in parsed_lines)
    
    # Pad shorter lines
    for i, fields in enumerate(parsed_lines):
        if len(fields) < max_cols:
            fields.extend([''] * (max_cols - len(fields)))
        parsed_lines[i] = fields
    
    # Create DataFrame
    df = pd.DataFrame(parsed_lines)
    
    return df


def read_file_auto(file_path):
    """
    Auto-detect separator and read file.
    """
    # Read first few lines to detect separator
    lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin1') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    
    if not lines:
        return pd.DataFrame()
    
    # Check for tabs
    if any('\t' in line for line in lines):
        return read_file_with_separator(file_path, '\t')
    
    # Check for commas
    if any(',' in line for line in lines):
        return read_file_with_separator(file_path, ',')
    
    # Check for semicolons
    if any(';' in line for line in lines):
        return read_file_with_separator(file_path, ';')
    
    # Default to space
    return read_file_with_separator(file_path, ' ')


def detect_column_type(series):
    """
    Detect the type of a column by checking if ALL non-empty values are of the same type.
    Returns: 'integer', 'numeric', 'boolean', 'text', 'mixed', or 'empty'
    """
    # Drop null/empty values for type detection
    non_null = series.dropna()
    non_empty = non_null[non_null.astype(str).str.strip() != '']
    
    if len(non_empty) == 0:
        return 'empty'
    
    str_values = non_empty.astype(str).str.strip()
    
    # Check if all values are numeric
    numeric_pattern = re.compile(r'^[-+]?[\d]*\.?[\d]+(?:[eE][-+]?\d+)?$')
    numeric_count = sum(1 for v in str_values if numeric_pattern.match(v))
    
    if numeric_count == len(str_values):
        # All values are numeric - check if integer or float
        try:
            numeric_values = pd.to_numeric(non_empty, errors='raise')
            if (numeric_values == numeric_values.astype(int)).all():
                return 'integer'
            else:
                return 'numeric'
        except:
            return 'numeric'
    
    # Check if all values are boolean
    bool_values = {'true', 'false', 'yes', 'no', 'y', 'n', '1', '0', 
                   'g', 'l', 'enriched', 'depleted', 'up', 'down'}
    lower_values = str_values.str.lower()
    if all(v in bool_values for v in lower_values):
        return 'boolean'
    
    # If no numeric values at all, all are text
    if numeric_count == 0:
        return 'text'
    
    # Mixed types (some numeric, some text)
    return 'mixed'


def auto_detect_and_convert(df, verbose=False):
    """
    Automatically detect and convert column types.
    Only assigns a type if ALL values in a column are consistent.
    Mixed columns are converted to string.
    """
    converted_df = df.copy()
    
    for col_idx, col in enumerate(converted_df.columns):
        col_type = detect_column_type(converted_df[col])
        
        if col_type == 'integer':
            # All values are integers
            converted_df[col] = pd.to_numeric(converted_df[col], errors='coerce').astype('Int64')
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: integer (consistent)")
                
        elif col_type == 'numeric':
            # All values are numeric (floats)
            converted_df[col] = pd.to_numeric(converted_df[col], errors='coerce')
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: numeric (consistent)")
                
        elif col_type == 'boolean':
            # All values are boolean
            bool_map = {
                'true': True, 'false': False,
                'yes': True, 'no': False,
                'y': True, 'n': False,
                '1': True, '0': False,
                'g': True, 'l': False,
                'enriched': True, 'depleted': False,
                'up': True, 'down': False
            }
            converted_df[col] = converted_df[col].astype(str).str.strip().str.lower().map(bool_map)
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: boolean (consistent)")
                
        elif col_type == 'text':
            # All values are text - ensure string type
            converted_df[col] = converted_df[col].astype(str)
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: text (string)")
                
        elif col_type == 'mixed':
            # Mixed types - convert to string
            converted_df[col] = converted_df[col].astype(str)
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: mixed types (converted to string)")
                
        elif col_type == 'empty':
            # Empty column - convert to string
            converted_df[col] = converted_df[col].astype(str)
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: empty (string)")
    
    return converted_df


def tsv_to_excel(output_file, input_specs, separator=None,
                 quotechar='"', skip_first_column=False, auto_convert_types=True,
                 global_header=None):
    """
    Convert multiple delimited files to an Excel file with separate sheets.
    
    Args:
        output_file: Path to output Excel file
        input_specs: List of input specifications (file[:sheet[:header_rows]])
        separator: Field separator (None for auto-detect)
        quotechar: Quote character
        skip_first_column: Skip first column
        auto_convert_types: Auto-detect and convert column types
        global_header: Global header setting (True/False/None) - overrides per-file setting
    """
    if not input_specs:
        print("Error: No input files provided.")
        return False

    try:
        used_sheet_names = set()
        processed_count = 0

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for input_spec in input_specs:
                file_path_str, custom_sheet_name, header_rows = parse_input_spec(input_spec)
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

                # Read file
                if separator is not None:
                    print(f"  Using separator: '{repr(separator)}'")
                    df = read_file_with_separator(file_path, separator)
                else:
                    print(f"  Auto-detecting separator...")
                    df = read_file_auto(file_path)
                
                # Determine final header_rows
                if global_header is not None:
                    # Global setting overrides everything
                    final_header_rows = 1 if global_header else 0
                else:
                    # Use per-file setting (default 0)
                    final_header_rows = header_rows
                
                print(f"  File: '{file_path.name}'")
                print(f"    Header rows: {final_header_rows}")
                print(f"    Columns: {len(df.columns)}")
                print(f"    Rows: {len(df)}")

                if df.empty:
                    print(f"Warning: '{file_path_str}' contains no data. Skipping...")
                    continue

                # Apply header rows
                if final_header_rows > 0:
                    # Use specified number of rows as header
                    if final_header_rows >= len(df):
                        print(f"Warning: Header rows ({final_header_rows}) >= data rows ({len(df)}). Using no header.")
                        has_header = False
                        df.columns = [f'Column_{i+1}' for i in range(len(df.columns))]
                    else:
                        # Combine multiple header rows if needed
                        if final_header_rows == 1:
                            df.columns = df.iloc[0].tolist()
                        else:
                            # For multiple header rows, join them
                            header_cols = []
                            for col_idx in range(len(df.columns)):
                                col_parts = []
                                for row_idx in range(final_header_rows):
                                    val = df.iloc[row_idx, col_idx]
                                    if str(val).strip():
                                        col_parts.append(str(val).strip())
                                header_cols.append('_'.join(col_parts) if col_parts else f'Column_{col_idx+1}')
                            df.columns = header_cols
                        
                        df = df.iloc[final_header_rows:].reset_index(drop=True)
                        has_header = True
                else:
                    # No header
                    has_header = False
                    df.columns = [f'Column_{i+1}' for i in range(len(df.columns))]

                # Optionally skip first column
                if skip_first_column and len(df.columns) > 1:
                    df = df.iloc[:, 1:]

                # Auto-detect and convert types
                if auto_convert_types:
                    print(f"    Detecting types...")
                    df = auto_detect_and_convert(df, verbose=True)
                    print(f"    ✓ Type conversion complete")

                # Write to Excel
                df.to_excel(writer, sheet_name=sheet_name, index=False, header=has_header)
                
                # Apply formatting if has header
                if has_header:
                    worksheet = writer.sheets[sheet_name]
                    
                    # Format header row
                    for cell in worksheet[1]:
                        cell.font = Font(bold=True)
                        cell.fill = PatternFill(start_color="DDDDDD", end_color="DDDDDD", fill_type="solid")
                        cell.alignment = Alignment(horizontal='center', vertical='center')
                    
                    # Auto-adjust column widths
                    for column in worksheet.columns:
                        max_length = 0
                        column_letter = column[0].column_letter
                        for cell in column:
                            if cell.value is not None:
                                max_length = max(max_length, len(str(cell.value)))
                        adjusted_width = min(max_length + 2, 50)
                        worksheet.column_dimensions[column_letter].width = adjusted_width
                    
                    # Freeze header row
                    worksheet.freeze_panes = 'A2'
                    
                    # Add auto-filter
                    worksheet.auto_filter.ref = worksheet.dimensions
                
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
  python tsv2excel.py output.xlsx data1.tsv

  # With header (1 row)
  python tsv2excel.py output.xlsx data1.tsv:Sheet1:1

  # No header explicitly
  python tsv2excel.py output.xlsx data1.tsv:Sheet1:0

  # Force tab separator
  python tsv2excel.py output.xlsx --separator tab data1.tsv:Sheet1:1

  # Multiple files with different header settings
  python tsv2excel.py output.xlsx file1.tsv:Sheet1:1 file2.tsv:Sheet2:0 file3.tsv:Sheet3:2

  # Disable automatic type detection
  python tsv2excel.py output.xlsx --no-type-detection data1.tsv
        """
    )

    parser.add_argument("output_file", help="Path to the output Excel file (.xlsx)")
    parser.add_argument("input_files", nargs='+',
                        help="Input files (format: file.tsv[:sheet_name[:header_rows]])")
    parser.add_argument("--separator", "--sep", dest="separator", default=None,
                        help="Field delimiter: 'tab', 'comma', 'semicolon', 'space', or custom character")
    parser.add_argument("--quotechar", dest="quotechar", default='"',
                        help="Character used to quote fields (default: \")")
    parser.add_argument("--skip_first_column", "-s", action="store_true",
                        help="Remove the first column from each input file")
    parser.add_argument("--no-type-detection", dest="no_type_detection", action="store_true",
                        help="Disable automatic type detection (all values stored as text)")
    parser.add_argument("--header", action="store_true",
                        help="Force header for all files (overrides per-file setting)")
    parser.add_argument("--no-header", action="store_true",
                        help="Force no header for all files (overrides per-file setting)")

    args = parser.parse_args()

    # Global header settings
    if args.header and args.no_header:
        print("Error: Cannot specify both --header and --no-header")
        sys.exit(1)
    
    global_header = None
    if args.header:
        global_header = True
    elif args.no_header:
        global_header = False

    # Get separator
    separator = get_separator(args.separator)

    output_file = args.output_file
    if not output_file.endswith('.xlsx'):
        output_file += '.xlsx'
        print(f"Note: Adding .xlsx extension to output file: {output_file}")

    success = tsv_to_excel(
        output_file,
        args.input_files,
        separator=separator,
        quotechar=args.quotechar,
        skip_first_column=args.skip_first_column,
        auto_convert_types=not args.no_type_detection,
        global_header=global_header
    )

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
