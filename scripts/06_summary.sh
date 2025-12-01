#!/bin/bash
#SBATCH --job-name=summary
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <base_dir> <amplicon> [gene]"
    echo "  base_dir: Base directory to process"
    echo "  amplicon: Amplicon type (e.g., rRNA or COI)"
    echo "  gene: Gene name (optional, e.g., 18S if rRNA amplicons are being processed)"
    exit 1
fi

# Assign command line arguments
base_dir=$1
amplicon=$2
gene=${3}  # Use $3 if provided, otherwise default to "18S"

source activate R

# Check if base_dir exists and is a directory
if [ ! -d "$base_dir" ]; then
    echo "Error: Directory '$base_dir' does not exist"
    exit 1
fi

echo "Processing directory: $base_dir"
echo "Amplicon: $amplicon"
echo "Gene: $gene"

Rscript nanopore-barcoding-ORC/scripts/auxiliary_code/amplicon_summary.R "$base_dir" "$amplicon" "$gene"