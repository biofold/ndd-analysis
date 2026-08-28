#!/bin/bash

# Usage: ./script.sh supercandidate.tsv cancer_genes.txt mondo_genes.txt

if [ $# -ne 3 ]; then
    echo "Usage: $0 <supercandidate.tsv> <cancer_genes.txt> <mondo_genes.txt>"
    echo ""
    echo "Note: Lines starting with # are treated as comments and skipped"
    echo "      Cancer and MONDO gene files can have multiple columns (first column = gene ID)"
    exit 1
fi

SUPERCANDIDATE_FILE="$1"
CANCER_GENES="$2"
MONDO_GENES="$3"

# Check if files exist
for file in "$SUPERCANDIDATE_FILE" "$CANCER_GENES" "$MONDO_GENES"; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file"
        exit 1
    fi
done

awk -F '\t' -v cancer_file="$CANCER_GENES" -v mondo_file="$MONDO_GENES" '
function format_number(n) {
    # Add commas to numbers for thousands
    if (n >= 1000000) {
        return sprintf("%d,%03d,%03d", int(n/1000000), int((n%1000000)/1000), n%1000)
    } else if (n >= 1000) {
        return sprintf("%d,%03d", int(n/1000), n%1000)
    }
    return sprintf("%d", n)
}

function is_comment(line) {
    # Check if line starts with # (optionally after whitespace)
    return (line ~ /^[[:space:]]*#/)
}

BEGIN {
    # Read cancer genes (first column only)
    while ((getline line < cancer_file) > 0) {
        # Skip comments and empty lines
        if (!is_comment(line) && line != "") {
            # Split by tab to get first column
            split(line, fields, "\t")
            gene = fields[1]
            
            # Remove any leading/trailing whitespace
            gsub(/^[[:space:]]+/, "", gene)
            gsub(/[[:space:]]+$/, "", gene)
            
            if (gene != "") {
                cancer_genes[gene] = 1
            }
        }
    }
    close(cancer_file)

    # Read MONDO genes (first column only)
    while ((getline line < mondo_file) > 0) {
        # Skip comments and empty lines
        if (!is_comment(line) && line != "") {
            # Split by tab to get first column
            split(line, fields, "\t")
            gene = fields[1]
            
            # Remove any leading/trailing whitespace
            gsub(/^[[:space:]]+/, "", gene)
            gsub(/[[:space:]]+$/, "", gene)
            
            if (gene != "") {
                mondo_genes[gene] = 1
            }
        }
    }
    close(mondo_file)

    # Initialize counters
    for (i = 0; i <= 5; i++) {
        mondo_count[i] = 0
        cancer_count[i] = 0
        total_count[i] = 0
    }

    # Initialize totals
    total_mondo = 0
    total_cancer = 0
    total_all = 0
    
    # Flag to track if header has been seen
    header_seen = 0
}

{
    # Skip comments and empty lines
    if (is_comment($0) || $0 == "") {
        next
    }

    # Skip header (first non-comment line)
    if (!header_seen) {
        header_seen = 1
        next
    }

    gene = $1
    score = $2

    if (score >= 0 && score <= 5) {
        total_count[score]++
        total_all++

        if (gene in cancer_genes) {
            cancer_count[score]++
            total_cancer++
        }

        if (gene in mondo_genes) {
            mondo_count[score]++
            total_mondo++
        }
    }
}

END {
    # Print header
    printf "%-10s\t%-25s\t%-25s\t%-25s\n", "Score", "MONDO Genes", "Cancer Genes", "Total Genes"
    #printf "%-10s\t%-25s\t%-25s\t%-25s\n", "-----", "-----------", "------------", "-----------"

    # Print rows
    for (score = 5; score >= 0; score--) {
        mondo_pct = (total_mondo > 0) ? (mondo_count[score] / total_mondo) * 100 : 0
        cancer_pct = (total_cancer > 0) ? (cancer_count[score] / total_cancer) * 100 : 0
        total_pct = (total_all > 0) ? (total_count[score] / total_all) * 100 : 0

        printf "%-10s\t%-25s\t%-25s\t%-25s\n",
            score,
            sprintf("%s (%.1f%%)", format_number(mondo_count[score]), mondo_pct),
            sprintf("%s (%.1f%%)", format_number(cancer_count[score]), cancer_pct),
            sprintf("%s (%.1f%%)", format_number(total_count[score]), total_pct)
    }

    # Print total row (sum of all rows)
    printf "%-10s\t%-25s\t%-25s\t%-25s\n",
        "Total",
        sprintf("%s (100.0%%)", format_number(total_mondo)),
        sprintf("%s (100.0%%)", format_number(total_cancer)),
        sprintf("%s (100.0%%)", format_number(total_all))
}
' "$SUPERCANDIDATE_FILE"
