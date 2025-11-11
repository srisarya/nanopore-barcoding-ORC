#!/bin/bash
#SBATCH --job-name=primer_search
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-169

# Array job wrapper for primer search
# Usage: sbatch --array=1-N primer_search_wrapper.sh <input_dir>
# Where N is the number of FASTA/FASTQ files to process

# Get arguments
INPUT_DIR="${1}"
PRIMER_FILE="nanopore-barcoding-ORC/adapters_hmms/amplicon_primer_sequences.fa"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' not found!"
    exit 1
fi

# Create output directories
OUTPUT_DIR="${INPUT_DIR}/primer_search_results"
mkdir -p "$OUTPUT_DIR"

# Find all FASTA/FASTQ files
mapfile -t FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( \
    -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o \
    -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort)

# Check if files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "Error: No FASTA files found in $INPUT_DIR"
    exit 1
fi

# Get the file for this array task
FILE_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
INPUT_FILE="${FILES[$FILE_INDEX]}"

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File index $FILE_INDEX not found in file list"
    exit 1
fi

# Get basename without extension
BASENAME=$(basename "$INPUT_FILE")
# Remove all extensions (.fastq.gz, .fasta.gz, etc.)
BASENAME_NO_EXT="${BASENAME%%.*}"
OUTPUT_CSV="${OUTPUT_DIR}/${BASENAME_NO_EXT}_primersearch.csv"

echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing file: $INPUT_FILE"
echo "Output: $OUTPUT_CSV"
echo "Primer file: $PRIMER_FILE"

# Activate seqkit environment
source activate seqkit

# Run the primer search script
python3 nanopore-barcoding-ORC/scripts/checks_balances/primer_amplicon_search.py \
    --primers "$PRIMER_FILE" \
    --file "$INPUT_FILE" \
    --output "$OUTPUT_CSV"

# Check if output was created
if [ -f "$OUTPUT_CSV" ]; then
    echo "Success! Results written to: $OUTPUT_CSV"
else
    echo "Error: Output file was not created!"
    exit 1
fi
