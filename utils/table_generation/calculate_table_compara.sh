#!/bin/bash
# Significance is adjusted p < 0.01 (strict), on unrounded values; the
# complement is >= 0.01. Same rule as the MOE score (scripts/2_supercandidate.py).

# Default column numbers
COL9=9
COL10=10
COL11=11
COLNF="NF"  # Last column

# Parse command line options
while getopts "a:b:c:d:h" opt; do
    case $opt in
        a) COL9=$OPTARG ;;
        b) COL10=$OPTARG ;;
        c) COL11=$OPTARG ;;
        d) COLNF=$OPTARG ;;
        h) 
            echo "Usage: $0 [options] <input_file>"
            echo "Options:"
            echo "  -a <col>  Column for cancer (default: 9)"
            echo "  -b <col>  Column for non_cancer (default: 10)"
            echo "  -c <col>  Column for gene_set (default: 11)"
            echo "  -d <col>  Column for cancer_vs_non_cancer (default: NF/last column)"
            echo "  -h        Show this help"
            exit 0
            ;;
        \?) 
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

# Shift to get the input file
shift $((OPTIND-1))

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 [options] <input_file>"
    echo "Use -h for help"
    exit 1
fi

INPUT_FILE=$1

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found"
    exit 1
fi

# Function to get the value of a column (handles NF specially)
get_col_value() {
    local col=$1
    if [ "$col" == "NF" ]; then
        echo "\$NF"
    else
        echo "\$$col"
    fi
}

# Build the awk column references
AWK_COL9=$(get_col_value $COL9)
AWK_COL10=$(get_col_value $COL10)
AWK_COL11=$(get_col_value $COL11)
AWK_COLNF=$(get_col_value $COLNF)

echo -e "Matrix\tCondition"
echo -e "\t"
#echo "======================================================================"
echo -e "Matrix 1:\tGeneSet: All "
#echo "======================================================================"
#echo ""

# Matrix 1: $11 condition not considered
# Row 1: GeneSet header
# Row 2: cancer<0.01
# Row 3: cancer>=0.01
# Col 1: empty/GeneSet
# Col 2: non-cancer>=0.01
# Col 3: non-cancer<0.01

# Cell (cancer<0.01, non-cancer>=0.01): count all, (count with $NF<0.01)
m1_r1c1_all=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10>=0.01 {count++} END {print count+0}" "$INPUT_FILE")
m1_r1c1_nf=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10>=0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer<0.01, non-cancer<0.01)
m1_r1c2_all=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m1_r1c2_nf=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer>=0.01, non-cancer>=0.01)
m1_r2c1_all=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10>=0.01 {count++} END {print count+0}" "$INPUT_FILE")
m1_r2c1_nf=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10>=0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer>=0.01, non-cancer<0.01)
m1_r2c2_all=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m1_r2c2_nf=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Print Matrix 1
echo -e "GeneSet\tnon-cancer>=0.01\tnon-cancer<0.01"
echo -e "cancer<0.01\t${m1_r1c1_all} (${m1_r1c1_nf})\t${m1_r1c2_all} (${m1_r1c2_nf})"
echo -e "cancer>=0.01\t${m1_r2c1_all} (${m1_r2c1_nf})\t${m1_r2c2_all} (${m1_r2c2_nf})"
echo -e "\t"

#echo "======================================================================"
echo -e "Matrix 2:\tGeneSet: <0.01"
#echo "======================================================================"

# Matrix 2: $11<0.01 condition
# Cell (cancer<0.01, non-cancer>=0.01)
m2_r1c1_all=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10>=0.01 && $AWK_COL11<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m2_r1c1_nf=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10>=0.01 && $AWK_COL11<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer<0.01, non-cancer<0.01)
m2_r1c2_all=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10<0.01 && $AWK_COL11<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m2_r1c2_nf=$(awk "NR>0 && $AWK_COL9<0.01 && $AWK_COL10<0.01 && $AWK_COL11<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer>=0.01, non-cancer>=0.01)
m2_r2c1_all=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10>=0.01 && $AWK_COL11<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m2_r2c1_nf=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10>=0.01 && $AWK_COL11<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Cell (cancer>=0.01, non-cancer<0.01)
m2_r2c2_all=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10<0.01 && $AWK_COL11<0.01 {count++} END {print count+0}" "$INPUT_FILE")
m2_r2c2_nf=$(awk "NR>0 && $AWK_COL9>=0.01 && $AWK_COL10<0.01 && $AWK_COL11<0.01 && $AWK_COLNF<0.01 {count++} END {print count+0}" "$INPUT_FILE")

# Print Matrix 2
echo -e "GeneSet<0.01\tnon-cancer>=0.01\tnon-cancer<0.01"
echo -e "cancer<0.01\t${m2_r1c1_all} (${m2_r1c1_nf})\t${m2_r1c2_all} (${m2_r1c2_nf})"
echo -e "cancer>=0.01\t${m2_r2c1_all} (${m2_r2c1_nf})\t${m2_r2c2_all} (${m2_r2c2_nf})"
echo -e "\t"

#echo "======================================================================"
line_count=`wc -l < "$INPUT_FILE" |awk '{print $1}' `
echo -e "Total lines:\t${line_count}"
echo ""
#echo "Note: Column mapping: cancer=col$COL9, non_cancer=col$COL10, gene_set=col$COL11, cancer_vs_non_cancer=col$COLNF"
#echo "Note: Numbers in parentheses indicate counts where cancer_vs_non_cancer<0.01"
