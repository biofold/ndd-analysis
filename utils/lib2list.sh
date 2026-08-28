#!/bin/bash

# Usage: ./extract_mondo_genes.sh <mondo_library.gmt>
# mondo_library.gmt: GMT format file with columns: class_code, class_name, gene1, gene2, ...

if [ $# -lt 1 ]; then
    echo "Usage: $0 <mondo_library.gmt>" >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  mondo_library.gmt  - GMT format file with columns:" >&2
    echo "                       Column 1: Class code (e.g., MONDO:0002051)" >&2
    echo "                       Column 2: Class name" >&2
    echo "                       Column 3+: Genes in the class" >&2
    echo "" >&2
    echo "Note: Lines starting with # are treated as comments and skipped" >&2
    echo "      Spaces will be replaced with underscores" >&2
    echo "      Output: gene<TAB>code1|code2|code3" >&2
    exit 1
fi

MONDO_LIBRARY="$1"

# Check if file exists
if [ ! -f "$MONDO_LIBRARY" ]; then
    echo "Error: File not found: $MONDO_LIBRARY" >&2
    exit 1
fi

# Extract unique genes with associated codes to stdout
awk -F '\t' '
function is_comment(line) {
    # Check if line starts with # (optionally after whitespace)
    return (line ~ /^[[:space:]]*#/)
}

{
    # Skip comments and empty lines
    if (is_comment($0) || $0 == "") {
        next
    }
    
    class_code = $1
    
    # Skip if missing class code
    if (class_code == "") {
        next
    }
    
    # Replace spaces with underscores in class code
    gsub(/ /, "_", class_code)
    
    # Process genes starting from column 3
    for (i = 3; i <= NF; i++) {
        gene = $i
        
        # Skip empty gene fields
        if (gene != "") {
            # Remove any spaces from gene names
            gsub(/ /, "_", gene)
            print gene "\t" class_code
        }
    }
}' "$MONDO_LIBRARY" | \
    sort -u | \
    awk -F '\t' '
    {
        gene = $1
        code = $2
        
        # Accumulate codes for each gene
        if (gene != prev_gene) {
            if (prev_gene != "") {
                print prev_gene "\t" codes
            }
            prev_gene = gene
            codes = code
        } else {
            codes = codes "|" code
        }
    }
    END {
        if (prev_gene != "") {
            print prev_gene "\t" codes
        }
    }'
