#!/bin/bash

# Usage: ./generate_summary.sh <input.tsv> [output.txt]
# input.tsv: File with 4 columns: gene, class, score, pipe-separated subgroups
# output.txt: Optional output file (default: stdout)

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.tsv> [output.txt]" >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  input.tsv   - File with 4 columns: gene, class, score, pipe-separated subgroups" >&2
    echo "  output.txt  - Optional output file (default: stdout)" >&2
    echo "" >&2
    echo "Note: Lines starting with # are treated as comments and skipped" >&2
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found: $INPUT_FILE" >&2
    exit 1
fi

# Process with awk
awk -F '\t' '
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
    # Initialize counters
    valid_classes["curated"] = 1
    valid_classes["candidate"] = 1
    valid_classes["no_evidence"] = 1
    
    # Initialize evidence categories
    categories[1] = "hc"
    categories[2] = "candidate"
    categories[3] = "sanchis"
    categories[4] = "develop"
    categories[5] = "neuro"
    categories[6] = "sfari+"
    categories[7] = "sfari"
    num_categories = 7
}

{
    # Skip comments and empty lines
    if (is_comment($0) || $0 == "") {
        next
    }
    
    # Check if we have at least 4 columns
    if (NF < 4) {
        next
    }
    
    gene = $1
    gene_class = $2
    score = $3
    evidence = $4
    
    # Skip lines that do not have a valid class
    if (!(gene_class in valid_classes)) {
        next
    }
    
    # Count genes by class
    class_counts[gene_class]++
    
    # Process evidence for curated and candidate classes
    if (gene_class == "curated" || gene_class == "candidate") {
        if (evidence != "none" && evidence != "") {
            # Split evidence by pipe
            n = split(evidence, evidence_array, "|")
            for (i = 1; i <= n; i++) {
                category = tolower(evidence_array[i])
                # Remove leading/trailing whitespace
                gsub(/^[[:space:]]+/, "", category)
                gsub(/[[:space:]]+$/, "", category)
                
                if (category != "") {
                    evidence_counts[gene_class, category]++
                }
            }
        }
    }
}

END {
    # Get class counts
    curated = (class_counts["curated"] != "") ? class_counts["curated"] : 0
    candidate = (class_counts["candidate"] != "") ? class_counts["candidate"] : 0
    no_evidence = (class_counts["no_evidence"] != "") ? class_counts["no_evidence"] : 0
    total = curated + candidate + no_evidence
    
    # Helper function to get counts
    # (awk does not have functions that can access arrays globally easily,
    #  so we will calculate directly)
    
    # Curated row
    hc_curated = evidence_counts["curated", "hc"]
    lc_curated = evidence_counts["curated", "candidate"]
    sanchis_curated = evidence_counts["curated", "sanchis"]
    develop_curated = evidence_counts["curated", "develop"]
    neuro_curated = evidence_counts["curated", "neuro"]
    sfari_plus_curated = evidence_counts["curated", "sfari+"]
    sfari_curated = evidence_counts["curated", "sfari"]
    
    # Candidate row
    hc_candidate = evidence_counts["candidate", "hc"]
    lc_candidate = evidence_counts["candidate", "candidate"]
    sanchis_candidate = evidence_counts["candidate", "sanchis"]
    develop_candidate = evidence_counts["candidate", "develop"]
    neuro_candidate = evidence_counts["candidate", "neuro"]
    sfari_plus_candidate = evidence_counts["candidate", "sfari+"]
    sfari_candidate = evidence_counts["candidate", "sfari"]
    
    # Totals
    total_hc = hc_curated + hc_candidate
    total_lc = lc_curated + lc_candidate
    total_sanchis = sanchis_curated + sanchis_candidate
    total_develop = develop_curated + develop_candidate
    total_neuro = neuro_curated + neuro_candidate
    total_sfari_plus = sfari_plus_curated + sfari_plus_candidate
    total_sfari = sfari_curated + sfari_candidate
    
    # Build output
    output = ""
    output = output "Gene Set\t# genes\tHC\tLC\tSanchis\tdevelop\tneuro\tsfari+\tsfari\n"
    output = output sprintf("curated\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        format_number(curated), format_number(hc_curated), format_number(lc_curated),
        format_number(sanchis_curated), format_number(develop_curated), 
        format_number(neuro_curated), format_number(sfari_plus_curated), 
        format_number(sfari_curated))
    
    output = output sprintf("candidate\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        format_number(candidate), format_number(hc_candidate), format_number(lc_candidate),
        format_number(sanchis_candidate), format_number(develop_candidate), 
        format_number(neuro_candidate), format_number(sfari_plus_candidate), 
        format_number(sfari_candidate))
    
    output = output sprintf("no_evidence\t%s\t-\t-\t-\t-\t-\t-\t-\n", 
        format_number(no_evidence))
    
    output = output sprintf("Total\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        format_number(total), format_number(total_hc), format_number(total_lc),
        format_number(total_sanchis), format_number(total_develop), 
        format_number(total_neuro), format_number(total_sfari_plus), 
        format_number(total_sfari))
    
    # Print output
    printf "%s", output
}' "$INPUT_FILE" > "${OUTPUT_FILE:-/dev/stdout}"

# If output file is specified, also display to stdout
if [ -n "$OUTPUT_FILE" ]; then
    cat "$OUTPUT_FILE"
fi
