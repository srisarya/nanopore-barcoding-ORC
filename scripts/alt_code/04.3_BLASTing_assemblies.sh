#!/bin/bash
#SBATCH --job-name=blast_array
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

# note, your number of files may vary, so adjust the array job line accordingly

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input directory provided"
    echo "Usage: sbatch $0 <in_directory>"
    echo "Example: sbatch $0 /path/to/indir"
    exit 1
fi

# Activate the environment
source activate blast # version 2.0.1
BLASTDB=/mnt/shared/datasets/databases/ncbi/nt

# Define input and output directories from command line argument
input_dir="$1"  # where the sequences are stored
basedir=$(dirname "$input_dir")  # Go up one level from indir to get base directory
output_dir="${basedir}/blast"

echo "Input directory: $input_dir"
echo "Base directory: $basedir"
echo "Output directory: $output_dir"

# make output dir
mkdir -p "${output_dir}"

files=( $(ls "${input_dir}"/*.fa | sort) )
total_files=${#files[@]}

# Get the filename for the current task (array is 1-indexed, bash arrays are 0-indexed)
file="${files[$((SLURM_ARRAY_TASK_ID-1))]}"
filename=$(basename "$file")  # Extract just the filename without path
echo "file being processed: $file"

# Run blast
blastn \
 -max_target_seqs 1 \
 -out "${output_dir}/${filename}_blast.tsv" \
 -outfmt "6 qseqid sseqid length evalue bitscore pident" \
 -db $BLASTDB \
 -query "$file"