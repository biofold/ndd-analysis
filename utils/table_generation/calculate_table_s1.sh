#!/bin/bash

# Usage: ./generate_comparison.sh <mondo_genes_codes.txt> <cancer_genes.txt>
# mondo_genes_codes.txt: File with gene and associated MONDO codes (pipe-separated)
# cancer_genes.txt: List of cancer genes (one per line)

if [ $# -lt 2 ]; then
    echo "Usage: $0 <cancer_genes.txt> <mondo_genes_codes.txt>" >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  cancer_genes.txt       - List of cancer genes (one per line)" >&2
    echo "  mondo_genes_codes.txt  - File with columns: gene, code1|code2|code3" >&2
    echo "" >&2
    echo "Note: Lines starting with # are treated as comments and skipped" >&2
    echo "      Output: Table sorted by MONDO gene count (descending) to stdout" >&2
    exit 1
fi

CANCER_GENES="$1"
MONDO_GENES_CODES="$2"

# Check if files exist and are not empty
for file in "$MONDO_GENES_CODES" "$CANCER_GENES"; do
    if [ ! -f "$file" ]; then
        echo "Error: File not found: $file" >&2
        exit 1
    fi
    if [ ! -s "$file" ]; then
        echo "Error: File is empty: $file" >&2
        exit 1
    fi
done

# Generate comparison table to stdout
awk -F '\t' -v cancer_file="$CANCER_GENES" -v mondo_file="$MONDO_GENES_CODES" '
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
    # Read cancer genes, skipping comments
    if (cancer_file == "") {
        print "Error: Cancer file not specified" > "/dev/stderr"
        exit 1
    }
    
    while ((getline line < cancer_file) > 0) {
        if (!is_comment(line) && line != "") {
            # Remove any leading/trailing whitespace
            gsub(/^[[:space:]]+/, "", line)
            gsub(/[[:space:]]+$/, "", line)
            # Replace spaces with underscores to match gene names
            gsub(/ /, "_", line)
            cancer_genes[line] = 1
        }
    }
    close(cancer_file)
    
    # Initialize arrays
    num_classes = 0
    total_mondo = 0
    
    # Read MONDO genes and codes, skipping comments
    if (mondo_file == "") {
        print "Error: MONDO file not specified" > "/dev/stderr"
        exit 1
    }
    
    while ((getline line < mondo_file) > 0) {
        if (is_comment(line) || line == "") {
            continue
        }
        
        split(line, fields, "\t")
        gene = fields[1]
        codes = fields[2]
        
        if (gene == "" || codes == "") {
            continue
        }
        
        # Store unique MONDO genes
        if (!(gene in mondo_genes)) {
            mondo_genes[gene] = 1
            total_mondo++
        }
        
        # Split codes by pipe
        n = split(codes, code_array, "|")
        for (i = 1; i <= n; i++) {
            code = code_array[i]
            
            # Store gene-code relationships
            gene_code[gene, code] = 1
            
            # Store unique codes
            if (!(code in all_codes)) {
                all_codes[code] = 1
                code_list[++num_classes] = code
            }
        }
    }
    close(mondo_file)
    
    # Initialize counters
    for (i = 1; i <= num_classes; i++) {
        code = code_list[i]
        mondo_count[code] = 0
        cancer_count[code] = 0
        total_count[code] = 0
    }
    
    # Count total cancer genes
    total_cancer = 0
    for (gene in cancer_genes) total_cancer++
    
    # Count MONDO genes per code
    for (gene in mondo_genes) {
        for (i = 1; i <= num_classes; i++) {
            code = code_list[i]
            if ((gene, code) in gene_code) {
                mondo_count[code]++
                total_count[code]++
            }
        }
    }
    
    # Count cancer genes per code (only if they are also MONDO genes)
    for (gene in cancer_genes) {
        if (gene in mondo_genes) {
            for (i = 1; i <= num_classes; i++) {
                code = code_list[i]
                if ((gene, code) in gene_code) {
                    cancer_count[code]++
                }
            }
        }
    }
    
    # Sort codes by MONDO count (descending)
    for (i = 1; i <= num_classes; i++) {
        for (j = i + 1; j <= num_classes; j++) {
            if (mondo_count[code_list[j]] > mondo_count[code_list[i]]) {
                temp = code_list[i]
                code_list[i] = code_list[j]
                code_list[j] = temp
            }
        }
    }
}

END {
    # Print header to stdout
    printf "%-30s\t%-25s\t%-25s\t%-25s\n", "MONDO Code", "Cancer Genes", "MONDO Genes", "Total Genes"
    #printf "%-30s\t%-25s\t%-25s\t%-25s\n", "----------", "-----------", "------------", "-----------"
    
    # Print rows for each code (sorted by MONDO count) to stdout
    for (i = 1; i <= num_classes; i++) {
        code = code_list[i]
        
        mondo_pct = (total_mondo > 0) ? (mondo_count[code] / total_mondo) * 100 : 0
        cancer_pct = (total_cancer > 0) ? (cancer_count[code] / total_cancer) * 100 : 0
        total_pct = (total_mondo > 0) ? (total_count[code] / total_mondo) * 100 : 0
        
        printf "%-30s\t%-25s\t%-25s\n", 
            code,
            sprintf("%s (%.1f%%)", format_number(cancer_count[code]), cancer_pct),
            sprintf("%s (%.1f%%)", format_number(mondo_count[code]), mondo_pct)
    }
    
    # Print total row to stdout
    printf "%-30s\t%-25s\t%-25s\n", 
        "Total",
        sprintf("%s (100.0%%)", format_number(total_cancer)),
        sprintf("%s (100.0%%)", format_number(total_mondo))
}' "$MONDO_GENES_CODES"
