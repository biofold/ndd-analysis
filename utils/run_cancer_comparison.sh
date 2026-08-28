#!/bin/bash
# utils/run_cancer_comparison.sh
# Complete cancer comparison pipeline: aggregate, Fisher test, and process results
# Modified to work with aggregate_pvals.py that uses stdout

set -e  # Exit on error
set -u  # Exit on undefined variable

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Complete cancer comparison pipeline: aggregate p-values, run Fisher test with BH correction,
and process results (select columns and sort).

Required arguments:
    -c, --cancer-file FILE       Cancer gene set enrichment file
    -n, --noncancer-file FILE    Non-cancer gene set enrichment file
    -a, --set-file FILE          Set genes enrichment file
    -n1, --n-cancer NUM          Number of genes in cancer set
    -n2, --n-noncancer NUM       Number of genes in non-cancer set
    -n3, --n-set NUM             Number of genes in set
    -o, --output-file FILE       Output processed file

Optional arguments:
    --pos1 COLS                  Columns for cancer contingency table (default: 3,4)
    --pos2 COLS                  Columns for non-cancer contingency table (default: 5,6)
    --pvalue-type TYPE           P-value type for BH correction (default: right)
    --select-cols COLS           Columns to select for final output (default: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
    --sort-col COL               Column to sort by (default: last column)
    --sort-order ORDER           Sort order: asc or desc (default: asc)
    --header-file FILE           File containing custom header (optional)
    --no-header                  Do not include header in output
    --keep-intermediate          Keep intermediate files (aggregate and fisher)
    --intermediate-dir DIR       Directory for intermediate files (default: same as output)
    -h, --help                   Show this help message

Example:
    $0 -c cancer_gs1.tsv -n noncancer_gs1.tsv -a gene_set1.tsv \\
       -n1 100 -n2 200 -n3 500 -o results/set1_processed.txt
EOF
    exit 0
}

# Default values
CANCER_FILE=""
NONCANCER_FILE=""
SET_FILE=""
N_CANCER=""
N_NONCANCER=""
N_SET=""
OUTPUT_FILE=""
POS1="3,4"
POS2="5,6"
PVALUE_TYPE="right"
SELECT_COLS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19"
SORT_COL="last"  # Default to sort by last column
SORT_ORDER="asc"
KEEP_INTERMEDIATE=false
INTERMEDIATE_DIR=""
HEADER_FILE=""
NO_HEADER=false

# Default header for selected columns
# This matches the fisher_bh.py output columns
DEFAULT_HEADER="GO_term\tGO_id\tcancer_overlap\tcancer_not_overlap\tnoncancer_overlap\tnoncancer_not_overlap\tset_overlap\tset_not_overlap\tcancer_pval\tnoncancer_pval\tset_pval\todds_ratio\tp_value_left\tp_value_right\tp_value_directional\tside\tp_value_adjusted"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--cancer-file) CANCER_FILE="$2"; shift 2 ;;
        -n|--noncancer-file) NONCANCER_FILE="$2"; shift 2 ;;
        -a|--set-file) SET_FILE="$2"; shift 2 ;;
        -n1|--n-cancer) N_CANCER="$2"; shift 2 ;;
        -n2|--n-noncancer) N_NONCANCER="$2"; shift 2 ;;
        -n3|--n-set) N_SET="$2"; shift 2 ;;
        -o|--output-file) OUTPUT_FILE="$2"; shift 2 ;;
        --pos1) POS1="$2"; shift 2 ;;
        --pos2) POS2="$2"; shift 2 ;;
        --pvalue-type) PVALUE_TYPE="$2"; shift 2 ;;
        --select-cols) SELECT_COLS="$2"; shift 2 ;;
        --sort-col) SORT_COL="$2"; shift 2 ;;
        --sort-order) SORT_ORDER="$2"; shift 2 ;;
        --header-file) HEADER_FILE="$2"; shift 2 ;;
        --no-header) NO_HEADER=true; shift ;;
        --keep-intermediate) KEEP_INTERMEDIATE=true; shift ;;
        --intermediate-dir) INTERMEDIATE_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: Unknown option $1" >&2; usage ;;
    esac
done

# Check required arguments
if [[ -z "$CANCER_FILE" || -z "$NONCANCER_FILE" || -z "$SET_FILE" || \
      -z "$N_CANCER" || -z "$N_NONCANCER" || -z "$N_SET" || -z "$OUTPUT_FILE" ]]; then
    echo "Error: All input files, counts, and output file are required" >&2
    usage
fi

# Check if input files exist
for file in "$CANCER_FILE" "$NONCANCER_FILE" "$SET_FILE"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: Input file not found: $file" >&2
        exit 1
    fi
done

# Set intermediate directory if not specified
if [[ -z "$INTERMEDIATE_DIR" ]]; then
    INTERMEDIATE_DIR=$(dirname "$OUTPUT_FILE")
fi

# Create directories
mkdir -p "$(dirname "$OUTPUT_FILE")"
mkdir -p "$INTERMEDIATE_DIR"

# Define intermediate files
OUTPUT_BASENAME=$(basename "$OUTPUT_FILE" .txt)
AGGREGATE_FILE="$INTERMEDIATE_DIR/${OUTPUT_BASENAME}_aggregate.txt"
FISHER_FILE="$INTERMEDIATE_DIR/${OUTPUT_BASENAME}_fisher.txt"

echo "=== Cancer Comparison Pipeline ==="
echo "Output file: $OUTPUT_FILE"
echo "Intermediate directory: $INTERMEDIATE_DIR"
echo ""

# Step 1: Aggregate p-values (using stdout redirection)
echo "Step 1: Aggregating p-values..."
python "$SCRIPT_DIR/aggregate_pvals.py" \
    "$CANCER_FILE" \
    "$NONCANCER_FILE" \
    "$SET_FILE" \
    "$N_CANCER" "$N_NONCANCER" "$N_SET" > "$AGGREGATE_FILE"

if [[ ! -s "$AGGREGATE_FILE" ]]; then
    echo "Error: Aggregation failed - output file is empty" >&2
    exit 1
fi

echo "  ✓ Aggregate file created: $AGGREGATE_FILE ($(wc -l < "$AGGREGATE_FILE") lines)"
echo ""

# Step 2: Fisher test with BH correction
echo "Step 2: Running Fisher test with BH correction..."
python "$SCRIPT_DIR/fisher_bh.py" \
    "$AGGREGATE_FILE" \
    "$POS1" "$POS2" \
    --pvalue_type "$PVALUE_TYPE" \
    --output "$FISHER_FILE"

if [[ ! -s "$FISHER_FILE" ]]; then
    echo "Error: Fisher test failed - output file is empty" >&2
    exit 1
fi

echo "  ✓ Fisher test results: $FISHER_FILE ($(wc -l < "$FISHER_FILE") lines)"
echo ""

# Step 3: Process results (select columns, merge first two, add header, and sort)
echo "Step 3: Processing results..."

# Create temporary files
TEMP_FILE="${OUTPUT_FILE}.tmp"
DATA_FILE="${OUTPUT_FILE}.data.tmp"
HEADER_FILE_TMP="${OUTPUT_FILE}.header.tmp"

# Determine the sort column number
if [[ "$SORT_COL" == "last" ]]; then
    # Calculate the number of selected columns
    IFS=',' read -ra COLS_ARRAY <<< "$SELECT_COLS"
    # After merging columns 1 and 2, we have one less column
    SORT_COL_NUM=$((${#COLS_ARRAY[@]}))
    echo "  Sorting by last column $SORT_COL_NUM (adjusted p-value)"
else
    SORT_COL_NUM="$SORT_COL"
    echo "  Sorting by column: $SORT_COL_NUM"
fi

# Determine sort flags
if [[ "$SORT_ORDER" == "desc" ]]; then
    SORT_FLAGS="-rgk"
    echo "  Sort order: descending"
else
    SORT_FLAGS="-gk"
    echo "  Sort order: ascending"
fi

# Check if Fisher file has header
FIRST_LINE=$(head -n 1 "$FISHER_FILE")
if [[ "$FIRST_LINE" == *"GO_term"* || "$FIRST_LINE" == *"odds_ratio"* || "$FIRST_LINE" == *"p_value"* ]]; then
    # File has header - skip it for data processing
    echo "  Input file has header, skipping it..."
    awk -F "\t" -v OFS='\t' \
        'NR>1 {
            # Merge first two columns with space
            go_term = $1 " " $2
            # Print merged term and selected columns
            print go_term, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $18, $19
        }' \
        "$FISHER_FILE" | \
    sort $SORT_FLAGS "$SORT_COL_NUM" > "$DATA_FILE"
else
    # No header - process all lines
    awk -F "\t" -v OFS='\t' \
        '{
            # Merge first two columns with space
            go_term = $1 " " $2
            # Print merged term and selected columns
            print go_term, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $18, $19
        }' \
        "$FISHER_FILE" | \
    sort $SORT_FLAGS "$SORT_COL_NUM" > "$DATA_FILE"
fi

# Determine header
if [[ "$NO_HEADER" == false ]]; then
    if [[ -n "$HEADER_FILE" && -f "$HEADER_FILE" ]]; then
        # Use custom header file
        HEADER=$(cat "$HEADER_FILE")
        echo "  Using custom header from: $HEADER_FILE"
    else
        # Generate header for merged output
        HEADER="GO_term\tcancer_overlap\tcancer_not_overlap\tnoncancer_overlap\tnoncancer_not_overlap\tset_overlap\tset_not_overlap\tcancer_pval\tnoncancer_pval\tset_pval\todds_ratio\tp_value_left\tp_value_right\tp_value_directional\tside\tp_value_adjusted"
        echo "  Generated header for merged columns"
    fi
    
    # Write header
    echo -e "$HEADER" > "$TEMP_FILE"
    # Append data
    cat "$DATA_FILE" >> "$TEMP_FILE"
else
    # No header - just data
    mv "$DATA_FILE" "$TEMP_FILE"
fi

# Move temp file to final output
if [[ -s "$TEMP_FILE" ]]; then
    mv "$TEMP_FILE" "$OUTPUT_FILE"
    OUTPUT_LINES=$(wc -l < "$OUTPUT_FILE")
    echo "  ✓ Processed results: $OUTPUT_FILE ($OUTPUT_LINES lines)"
    
    # Show preview
    if [[ $OUTPUT_LINES -gt 0 ]]; then
        echo "  Preview (first 5 lines):"
        head -n 5 "$OUTPUT_FILE" | while read -r line; do
            echo "    $line"
        done
    fi
else
    echo "Error: Processing failed - output file is empty" >&2
    rm -f "$TEMP_FILE" "$DATA_FILE" "$HEADER_FILE_TMP"
    exit 1
fi

# Clean up temporary files
rm -f "$DATA_FILE" "$HEADER_FILE_TMP"

echo ""

# Clean up intermediate files if not keeping them
if [[ "$KEEP_INTERMEDIATE" == false ]]; then
    echo "Cleaning up intermediate files..."
    rm -f "$AGGREGATE_FILE" "$FISHER_FILE"
    echo "  Removed intermediate files"
else
    echo "Keeping intermediate files:"
    echo "  - $AGGREGATE_FILE"
    echo "  - $FISHER_FILE"
fi

echo ""
echo "=== Pipeline Complete ==="
echo "Final output: $OUTPUT_FILE"
