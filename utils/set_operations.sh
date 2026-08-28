#!/bin/bash

# Complete set operations script with full separator support
# Usage: ./set_operations.sh <file1> <file2> <key_columns_file1> <key_columns_file2> [options]

# Default values
OPERATION="intersect"
WALL="TRUE"
DELIMITER=""

# Arrays for positional arguments
positional=()

# Parse command line
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--operation)
            if [[ -z "$2" ]]; then
                echo "Error: --operation requires a value" >&2
                exit 1
            fi
            OPERATION="$2"
            shift 2
            ;;
        -w|--wall)
            if [[ -z "$2" ]]; then
                echo "Error: --wall requires a value" >&2
                exit 1
            fi
            WALL="$2"
            shift 2
            ;;
        -d|--delimiter)
            # Delimiter may be empty, so we need to check that a second argument exists
            if [[ $# -lt 2 ]]; then
                echo "Error: --delimiter requires a value (even if empty)" >&2
                exit 1
            fi
            DELIMITER="$2"
            shift 2
            ;;
        -h|--help)
            cat << 'EOF'
Usage: $0 <file1> <file2> <key_columns_file1> <key_columns_file2> [options]

Positional arguments:
  file1            - First input file
  file2            - Second input file
  key_columns_file1 - Comma-separated column numbers for keys in file1 (1-based)
  key_columns_file2 - Comma-separated column numbers for keys in file2 (1-based)

Options:
  -o, --operation   Operation: intersect (default), diff12, diff21
                     Aliases:
                       intersect: intersect, intersection, both, INT, INTER
                       diff12: diff12, s1-s2, difference, diff, DIFF
                       diff21: diff21, s2-s1, reverse_diff, rdiff
  -w, --wall        Wall output: TRUE (default) or FALSE
                     Aliases: TRUE, True, true, T, t, 1, YES, Yes, yes, Y, y
                              FALSE, False, false, F, f, 0, NO, No, no, N, n
  -d, --delimiter   Field delimiter:
                       "" or "whitespace" = space+tab together (default)
                       "tab" = tab only
                       " " = space only
                       any other = custom delimiter
  -h, --help        Show this help

Examples:
  $0 file1.txt file2.txt 1 1
  $0 file1.txt file2.txt 1 1 -o intersect -w TRUE -d ""
  $0 file1.txt file2.txt 1 1 --operation DIFF --wall false
  $0 file1.txt file2.txt 1 1 -o diff21 -w FALSE -d tab
EOF
            exit 0
            ;;
        --)
            shift
            # Remaining arguments are positional
            while [[ $# -gt 0 ]]; do
                positional+=("$1")
                shift
            done
            break
            ;;
        -*)
            echo "Error: Unknown option: $1" >&2
            exit 1
            ;;
        *)
            positional+=("$1")
            shift
            ;;
    esac
done

# Check we have exactly 4 positional arguments
if [[ ${#positional[@]} -ne 4 ]]; then
    echo "Error: Exactly 4 positional arguments required (file1, file2, key_columns_file1, key_columns_file2)" >&2
    exit 1
fi

FILE1="${positional[0]}"
FILE2="${positional[1]}"
POS1="${positional[2]}"
POS2="${positional[3]}"

# Check files
if [ ! -f "$FILE1" ]; then echo "Error: $FILE1 not found" >&2; exit 1; fi
if [ ! -f "$FILE2" ]; then echo "Error: $FILE2 not found" >&2; exit 1; fi

# Normalize operation (convert to lowercase first)
OPERATION=$(echo "$OPERATION" | tr 'A-Z' 'a-z')

# Map operation aliases
case "$OPERATION" in
    intersect|intersection|both|int|inter)
        OPERATION="intersect"
        ;;
    diff12|s1-s2|difference|diff|s1_s2)
        OPERATION="diff12"
        ;;
    diff21|s2-s1|reverse_diff|rdiff|s2_s1|reverse)
        OPERATION="diff21"
        ;;
    *)
        echo "Error: Unknown operation: $OPERATION" >&2
        echo "Valid operations: intersect, diff12, diff21" >&2
        exit 1
        ;;
esac

# Normalize wall parameter
WALL_LOWER=$(echo "$WALL" | tr 'A-Z' 'a-z')

case "$WALL_LOWER" in
    true|t|1|yes|y)
        WALL="TRUE"
        ;;
    false|f|0|no|n)
        WALL="FALSE"
        ;;
    *)
        echo "Error: Unknown wall value: $WALL" >&2
        echo "Valid values: TRUE, FALSE (or their variations)" >&2
        exit 1
        ;;
esac

# Set field separator based on delimiter
if [ -z "$DELIMITER" ] || [ "$DELIMITER" = "whitespace" ] || [ "$DELIMITER" = "WHITESPACE" ]; then
    AWK_FS='[ \t]+'
elif [ "$DELIMITER" = "tab" ] || [ "$DELIMITER" = "TAB" ]; then
    AWK_FS='\t'
elif [ "$DELIMITER" = " " ]; then
    AWK_FS=' '
else
    AWK_FS="$DELIMITER"
fi

# Process with awk
awk -F "$AWK_FS" \
    -v op="$OPERATION" \
    -v wall="$WALL" \
    -v pos1="$POS1" \
    -v pos2="$POS2" \
    -v f1="$FILE1" \
    -v f2="$FILE2" '
BEGIN {
    split(pos1, c1, ",")
    split(pos2, c2, ",")
    
    # Read file1
    while ((getline line < f1) > 0) {
        if (line ~ /^[[:space:]]*$/) continue
        if (line ~ /^#/) continue
        
        n = split(line, fields, FS)
        
        key = ""
        for (i = 1; i <= length(c1); i++) {
            col = c1[i]
            if (col >= 1 && col <= n) {
                if (key != "") key = key "|"
                key = key fields[col]
            }
        }
        
        if (key != "") {
            data1[key] = data1[key] line "\n"
        }
    }
    close(f1)
    
    # Read file2
    while ((getline line < f2) > 0) {
        if (line ~ /^[[:space:]]*$/) continue
        if (line ~ /^#/) continue
        
        n = split(line, fields, FS)
        
        key = ""
        for (i = 1; i <= length(c2); i++) {
            col = c2[i]
            if (col >= 1 && col <= n) {
                if (key != "") key = key "|"
                key = key fields[col]
            }
        }
        
        if (key != "") {
            data2[key] = data2[key] line "\n"
        }
    }
    close(f2)
    
    # Perform operation
    if (op == "intersect") {
        for (key in data1) {
            if (key in data2) {
                n1 = split(data1[key], lines1, "\n")
                n2 = split(data2[key], lines2, "\n")
                
                for (i = 1; i <= n1; i++) {
                    if (lines1[i] == "") continue
                    if (wall == "TRUE") {
                        for (j = 1; j <= n2; j++) {
                            if (lines2[j] == "") continue
                            print lines1[i] "\t" lines2[j]
                        }
                    } else {
                        print lines1[i]
                    }
                }
            }
        }
    } else if (op == "diff12") {
        for (key in data1) {
            if (!(key in data2)) {
                n1 = split(data1[key], lines1, "\n")
                for (i = 1; i <= n1; i++) {
                    if (lines1[i] != "") print lines1[i]
                }
            }
        }
    } else if (op == "diff21") {
        for (key in data2) {
            if (!(key in data1)) {
                n2 = split(data2[key], lines2, "\n")
                for (i = 1; i <= n2; i++) {
                    if (lines2[i] != "") print lines2[i]
                }
            }
        }
    }
}
'
