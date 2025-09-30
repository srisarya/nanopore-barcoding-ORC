#!/bin/bash

# Usage: ./adapter_search_dir.sh <adapter_file> <input_directory> <output_file>
# adapter_file: FASTA file with adapter sequences
# input_directory: directory containing .gz files to search
# output_file: results will be written here

if [ $# -ne 3 ]; then
    echo "Usage: $0 <adapter_file> <input_directory> <output_file>"
    echo "Example: $0 adapters.fa /path/to/fastq_files/ results.txt"
    exit 1
fi

ADAPTER_FILE="$1"
INPUT_DIR="$2"
OUTPUT_FILE="$3"

# Check if adapter file exists
if [ ! -f "$ADAPTER_FILE" ]; then
    echo "Error: Adapter file '$ADAPTER_FILE' not found"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' not found"
    exit 1
fi

# Remove trailing slash from directory path if present
INPUT_DIR="${INPUT_DIR%/}"

echo "Adapter file: $ADAPTER_FILE"
echo "Input directory: $INPUT_DIR"
echo "Output file: $OUTPUT_FILE"
echo ""

# Initialize output file with header
echo "# Adapter search results - $(date)" > "$OUTPUT_FILE"
echo "# Format: ADAPTER_SEQUENCE | FILENAME | COUNT" >> "$OUTPUT_FILE"
echo "# Input directory: $INPUT_DIR" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Read adapters line by line and track current adapter label
current_label=""
while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    if [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]]; then
        continue
    fi
    
    # If it's a FASTA header, store the label and continue
    if [[ "$line" =~ ^\> ]]; then
        current_label=$(echo "$line" | sed 's/^>//')
        continue
    fi
    
    # Remove any whitespace from sequence
    adapter=$(echo "$line" | tr -d '[:space:]')
    
    # Skip empty sequences
    if [[ -z "$adapter" ]]; then
        continue
    fi
    
    echo "Searching for adapter: $current_label ($adapter)"
    echo "## Adapter: $current_label" >> "$OUTPUT_FILE"
    
    # Search in all .gz files in the specified directory that don't contain 'unknown' and are not empty
    file_count=0
    for file in "$INPUT_DIR"/*.gz; do
        # Check if glob matched any files
        if [ ! -e "$file" ]; then
            continue
        fi
        
        # Skip files with 'unknown' in name and empty files
        if [[ "$(basename "$file")" =~ unknown ]] || [[ ! -s "$file" ]]; then
            continue
        fi
        
        # Search for adapter in file using grep
        count=$(zless "$file" | grep -c "$adapter" 2>/dev/null || echo "0")
        
        # Get just the filename for output
        filename=$(basename "$file")
        
        # Output all results (including 0 counts)
        echo "$current_label | $filename | $count" >> "$OUTPUT_FILE"
        if [ "$count" -gt 0 ]; then
            echo "  Found $count matches in $filename"
        else
            echo "  No matches in $filename"
        fi
        
        ((file_count++))
    done
    
    if [ $file_count -eq 0 ]; then
        echo "  No .gz files found in directory"
        echo "$current_label | NO_FILES_FOUND | 0" >> "$OUTPUT_FILE"
    fi
    
    echo "" >> "$OUTPUT_FILE"
    
done < "$ADAPTER_FILE"

echo ""
echo "Search complete. Results written to: $OUTPUT_FILE"
