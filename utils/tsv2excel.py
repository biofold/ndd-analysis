#!/usr/bin/env python3
"""
Convert multiple delimited files to an Excel file with each file as a separate sheet.
Sheet names can be derived from the filenames (without extension) or specified explicitly.
Automatically detects and converts data types to avoid "number stored as text" issues.

Usage:
  python tsv2excel.py output.xlsx file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --header file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --separator ',' file1.csv file2.csv
  python tsv2excel.py output.xlsx --quotechar '"' file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --header --skip_first_column file1.tsv file2.tsv
  python tsv2excel.py output.xlsx --no-type-detection file1.tsv file2.tsv
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
    Parse input specification to extract file path and optional sheet name.
    Format: file.tsv or file.tsv:sheet_name or file.tsv:sheet_name:header
    """
    parts = input_spec.split(':')
    file_path = parts[0]
    sheet_name = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
    has_header = parts[2].strip().lower() in ['header', 'true', 'yes', '1'] if len(parts) > 2 else None
    
    return file_path, sheet_name, has_header


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


def get_value_type(value):
    """
    Determine the type of a value: 'numeric', 'text', or 'empty'
    """
    if value is None or (isinstance(value, str) and value.strip() == ''):
        return 'empty'
    
    value_str = str(value).strip()
    
    # Check if numeric
    numeric_pattern = re.compile(r'^[-+]?[\d]*\.?[\d]+(?:[eE][-+]?\d+)?$')
    if numeric_pattern.match(value_str):
        return 'numeric'
    
    # Check if boolean
    bool_values = {'true', 'false', 'yes', 'no', 'y', 'n', '1', '0'}
    if value_str.lower() in bool_values:
        return 'boolean'
    
    return 'text'


def looks_like_header(first_row, second_row):
    """
    Determine if the first row looks like a header by comparing types with the second row.
    Generic approach: if the type pattern of the first row differs from the second row,
    the first row is likely a header.
    """
    if not first_row or not second_row:
        return False
    
    # Get types for each field in both rows
    first_types = [get_value_type(v) for v in first_row]
    second_types = [get_value_type(v) for v in second_row]
    
    # Remove 'empty' types for comparison
    first_non_empty = [t for t in first_types if t != 'empty']
    second_non_empty = [t for t in second_types if t != 'empty']
    
    if not first_non_empty or not second_non_empty:
        return False
    
    # Check if type patterns are different
    type_pattern_differs = first_non_empty != second_non_empty
    
    # Check specific cases:
    
    # Case 1: First row all text, second row has some numeric/boolean
    first_all_text = all(t == 'text' for t in first_non_empty)
    second_has_numeric = any(t in ['numeric', 'boolean'] for t in second_non_empty)
    
    if first_all_text and second_has_numeric:
        return True
    
    # Case 2: First row has different types than second row
    if type_pattern_differs:
        # Count type occurrences
        first_text_count = sum(1 for t in first_non_empty if t == 'text')
        first_numeric_count = sum(1 for t in first_non_empty if t == 'numeric')
        second_text_count = sum(1 for t in second_non_empty if t == 'text')
        second_numeric_count = sum(1 for t in second_non_empty if t == 'numeric')
        
        # If first row is mostly text and second row has more numbers
        if first_text_count > first_numeric_count and second_numeric_count > first_numeric_count:
            return True
        
        # If first row is mostly text and second row is also text but different pattern
        if first_text_count == len(first_non_empty) and second_text_count == len(second_non_empty):
            # Both all text - check if they look different in some way
            # Headers tend to be shorter and have no special characters
            first_avg_length = np.mean([len(str(v).strip()) for v in first_row if str(v).strip()])
            second_avg_length = np.mean([len(str(v).strip()) for v in second_row if str(v).strip()])
            
            # Headers are typically shorter
            if first_avg_length < second_avg_length * 0.7:
                return True
    
    # Case 3: First row has some text, second row all text
    # This is ambiguous - could be header or data
    # Only consider as header if first row has very short values
    if first_all_text and all(t == 'text' for t in second_non_empty):
        first_avg_length = np.mean([len(str(v).strip()) for v in first_row if str(v).strip()])
        second_avg_length = np.mean([len(str(v).strip()) for v in second_row if str(v).strip()])
        
        # If first row values are much shorter, likely header
        if first_avg_length < second_avg_length * 0.5:
            return True
    
    return False


def detect_separator_and_structure(file_path):
    """
    Detect the separator and structure of the file.
    Returns: (separator, has_header, num_columns)
    """
    lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i >= 20:  # Read first 20 lines for detection
                    break
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin1') as f:
            for i, line in enumerate(f):
                if i >= 20:
                    break
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    
    if not lines:
        return '\t', False, 1
    
    # Try different separators
    separators = ['\t', ',', ';', ' ']
    best_separator = '\t'
    best_consistency = 0
    best_num_cols = 1
    
    for sep in separators:
        col_counts = []
        for line in lines:
            if sep == ' ':
                # For space separator, split on multiple spaces
                fields = re.split(r'\s+', line.strip())
            else:
                fields = line.split(sep)
            col_counts.append(len(fields))
        
        if col_counts:
            # Check consistency
            most_common = max(set(col_counts), key=col_counts.count)
            consistency = col_counts.count(most_common) / len(col_counts)
            
            if consistency > best_consistency:
                best_consistency = consistency
                best_separator = sep
                best_num_cols = most_common
    
    # Detect header by comparing first two rows
    has_header = False
    if len(lines) > 1:
        if best_separator == ' ':
            first_fields = re.split(r'\s+', lines[0].strip())
            second_fields = re.split(r'\s+', lines[1].strip())
        else:
            first_fields = lines[0].split(best_separator)
            second_fields = lines[1].split(best_separator)
        
        has_header = looks_like_header(first_fields, second_fields)
    
    return best_separator, has_header, best_num_cols


def read_file_robust(file_path):
    """
    Robust file reader that handles various formats.
    """
    separator, has_header, num_cols = detect_separator_and_structure(file_path)
    
    lines = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    except UnicodeDecodeError:
        with open(file_path, 'r', encoding='latin1') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.strip():
                    lines.append(line)
    
    # Parse lines based on detected separator
    parsed_lines = []
    for line in lines:
        if separator == ' ':
            # Split on spaces, but keep pipe-separated values together
            fields = re.split(r'\s+', line.strip())
        else:
            fields = line.split(separator)
        
        # Pad or truncate to consistent number of columns
        if len(fields) < num_cols:
            fields.extend([''] * (num_cols - len(fields)))
        elif len(fields) > num_cols:
            fields = fields[:num_cols]
        
        parsed_lines.append(fields)
    
    # Create DataFrame
    df = pd.DataFrame(parsed_lines)
    
    # If has header, set first row as column names
    if has_header and len(df) > 1:
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
    
    return df, has_header


def detect_column_type(series):
    """
    Detect the appropriate data type for a column.
    Returns: 'numeric', 'integer', 'boolean', or 'text'
    """
    # Drop null values for type detection
    non_null = series.dropna()
    if len(non_null) == 0:
        return 'text'
    
    # Convert all values to strings for analysis
    str_values = non_null.astype(str).str.strip()
    
    # Check for boolean values
    bool_values = {'true', 'false', 'yes', 'no', 'y', 'n', '1', '0', 
                   'g', 'l', 'enriched', 'depleted', 'up', 'down'}
    unique_lower = set(str_values.str.lower().unique())
    if len(unique_lower) <= 2 and unique_lower.issubset(bool_values):
        return 'boolean'
    
    # Check for numeric values
    numeric_pattern = re.compile(r'^[-+]?[\d]*\.?[\d]+(?:[eE][-+]?\d+)?$')
    numeric_count = sum(1 for v in str_values if numeric_pattern.match(v))
    
    if numeric_count == len(str_values):
        # All values are numeric
        try:
            numeric_series = pd.to_numeric(non_null, errors='raise')
            # Check if all values are integers
            if (numeric_series == numeric_series.astype(int)).all():
                return 'integer'
            else:
                return 'numeric'
        except (ValueError, TypeError):
            return 'numeric'
    elif numeric_count >= len(str_values) * 0.8:
        # At least 80% numeric - likely a numeric column with some text
        return 'numeric'
    
    return 'text'


def convert_column_type(series, col_type):
    """
    Convert a column to the detected type.
    """
    if col_type == 'numeric':
        return pd.to_numeric(series, errors='coerce')
    elif col_type == 'integer':
        # Use Int64 to handle NaN values
        numeric = pd.to_numeric(series, errors='coerce')
        return numeric.astype('Int64')
    elif col_type == 'boolean':
        bool_map = {
            'true': True, 'false': False,
            'yes': True, 'no': False,
            'y': True, 'n': False,
            '1': True, '0': False,
            'g': True, 'l': False,
            'enriched': True, 'depleted': False,
            'up': True, 'down': False
        }
        return series.astype(str).str.strip().str.lower().map(bool_map)
    else:
        # Keep as text
        return series.astype(str)


def auto_detect_and_convert(df, verbose=False):
    """
    Automatically detect and convert column types in a DataFrame.
    """
    converted_df = df.copy()
    
    for col_idx, col in enumerate(converted_df.columns):
        col_type = detect_column_type(converted_df[col])
        
        if col_type != 'text':
            converted_df[col] = convert_column_type(converted_df[col], col_type)
            
            if verbose:
                col_name = col if isinstance(col, str) else f"Column_{col_idx+1}"
                print(f"    {col_name}: {col_type}")
    
    return converted_df


def tsv_to_excel(output_file, input_specs, header=None, separator=None,
                 quotechar='"', skip_first_column=False, auto_convert_types=True):
    """
    Convert multiple delimited files to an Excel file with separate sheets.
    """
    if not input_specs:
        print("Error: No input files provided.")
        return False

    try:
        used_sheet_names = set()
        processed_count = 0

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for input_spec in input_specs:
                file_path_str, custom_sheet_name, has_header_override = parse_input_spec(input_spec)
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

                # Read file using robust reader
                df, has_header = read_file_robust(file_path)
                
                # Override header if specified
                if has_header_override is not None:
                    has_header = has_header_override
                elif header is not None:
                    has_header = header
                
                print(f"  File: '{file_path.name}'")
                print(f"    Header detected: {'Yes' if has_header else 'No'}")
                print(f"    Columns: {len(df.columns)}")
                print(f"    Rows: {len(df)}")

                if df.empty:
                    print(f"Warning: '{file_path_str}' contains no data. Skipping...")
                    continue

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
  # Auto-detect everything (recommended)
  python tsv2excel.py output.xlsx data1.tsv data2.tsv

  # Force tab separator
  python tsv2excel.py output.xlsx --separator '\\t' data1.tsv data2.tsv

  # Force header for all files
  python tsv2excel.py output.xlsx --header data1.tsv data2.tsv

  # Force no header for all files
  python tsv2excel.py output.xlsx --no-header data1.tsv data2.tsv

  # Custom sheet names with header override
  python tsv2excel.py output.xlsx data1.tsv:Table1:header data2.tsv:Table2:noheader

  # Disable automatic type detection
  python tsv2excel.py output.xlsx --no-type-detection data1.tsv data2.tsv
        """
    )

    parser.add_argument("output_file", help="Path to the output Excel file (.xlsx)")
    parser.add_argument("input_files", nargs='+',
                        help="Input files (format: file.tsv or file.tsv:sheet_name or file.tsv:sheet_name:header)")
    parser.add_argument("--header", action="store_true",
                        help="Treat first row as header for all files")
    parser.add_argument("--no-header", action="store_true",
                        help="Do not treat first row as header for any file")
    parser.add_argument("--separator", "--sep", dest="separator", default=None,
                        help="Field delimiter (default: auto-detect)")
    parser.add_argument("--quotechar", dest="quotechar", default='"',
                        help="Character used to quote fields (default: \")")
    parser.add_argument("--skip_first_column", "-s", action="store_true",
                        help="Remove the first column from each input file")
    parser.add_argument("--no-type-detection", dest="no_type_detection", action="store_true",
                        help="Disable automatic type detection (all values stored as text)")

    args = parser.parse_args()

    # Determine header mode
    if args.header and args.no_header:
        print("Error: Cannot specify both --header and --no-header")
        sys.exit(1)
    
    header_mode = None  # Auto-detect
    if args.header:
        header_mode = True
    elif args.no_header:
        header_mode = False

    output_file = args.output_file
    if not output_file.endswith('.xlsx'):
        output_file += '.xlsx'
        print(f"Note: Adding .xlsx extension to output file: {output_file}")

    success = tsv_to_excel(
        output_file,
        args.input_files,
        header=header_mode,
        separator=args.separator,
        quotechar=args.quotechar,
        skip_first_column=args.skip_first_column,
        auto_convert_types=not args.no_type_detection
    )

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
