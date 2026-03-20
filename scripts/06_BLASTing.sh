#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --mem=24G
#SBATCH --cpus-per-task=24
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch $0 <dataset> <gene>"
    echo "Example: sbatch $0 Lakes_day1 COI_gene"
    exit 1
fi

# Activate the environment
source activate blast # version 2.0.1
BLASTDB=/mnt/shared/datasets/databases/ncbi/nt

# Define directories from command line arguments
dataset="$1"
gene="$2"
output_dir="${dataset}/BLAST"

echo "Dataset: $dataset"
echo "Gene: $gene"
echo "Output directory: $output_dir"

# Make output directory
mkdir -p "${output_dir}"

# Concatenate only top-level FASTA files (avoid subdirectories)
concatenated_file="${output_dir}/${dataset}_${gene}_concatenated.fa"
echo "Concatenating top-level .fa/.fasta files into: $concatenated_file"

echo "Files being included:"
find "${dataset}/${gene}" \
    -mindepth 2 -maxdepth 2 \
    -type f \( -name "*.fa" -o -name "*.fasta" \) -print

# Actual concatenation
find "${dataset}/${gene}" \
    -mindepth 2 -maxdepth 2 \
    -type f \( -name "*.fa" -o -name "*.fasta" \) \
    -exec cat {} + > "$concatenated_file"

# Check if concatenation was successful
if [ ! -s "$concatenated_file" ]; then
    echo "Error: Concatenated file is empty or does not exist"
    echo "No top-level .fa or .fasta files found in: ${dataset}/${gene}/*/"
    exit 1
fi

echo "Running BLAST on concatenated file"

# Run blast
blastn \
 -max_target_seqs 500 \
 -out "${output_dir}/${dataset}_${gene}_blast.tsv" \
 -outfmt "6 qseqid qlen sseqid evalue bitscore pident staxids" \
 -db $BLASTDB \
 -num_threads 12 \
 -query "$concatenated_file"

echo "BLAST completed successfully"

# filter for top 5 e-value hits
# Claude Sonnet v4.5 was used to help debug ask lines in this section
sort -k1,1 -k4,4g "${output_dir}/${dataset}_${gene}_blast.tsv" \
| awk '{ key=$1; if (key!=prev) {count=0; prev=key} count++; if (count<=5) print }' \
> "${output_dir}/${dataset}_${gene}_blast_top5.tsv"

echo "BLAST filtered successfully"
