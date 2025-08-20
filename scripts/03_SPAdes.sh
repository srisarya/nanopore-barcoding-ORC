#!/bin/bash
#SBATCH --job-name=SPAdes_array
#SBATCH --mem=24G
#SBATCH --cpus-per-task=24
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input directory provided"
    echo "Usage: sbatch $0 <demuxed_SP27_directory>"
    echo "Example: sbatch $0 /path/to/demuxed/SP27"
    exit 1
fi

# Define input and output directories from command line argument
input_dir="$1"  # SP27 directory from command line
basedir=$(dirname $(dirname "$input_dir"))  # Go up two levels from SP27 to get base directory
output_base_dir="${basedir}/assembled"

echo "Input directory: $input_dir"
echo "Base directory: $basedir"
echo "Output base directory: $output_base_dir"

# make output basedir : megahit will complain if this does not exist
mkdir -p "${output_base_dir}"

# Get the filename for the current task
# file=$(ls "${input_dir}"/*.fastq.gz | perl -ne 'print if /SP27_(\d+)_SP5_\1_/' | sed -n "${SLURM_ARRAY_TASK_ID}p")
file=$(ls "${input_dir}"/*.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "file being processed:" ${file}

# Check if file exists
if [ ! -f "$file" ]; then
    echo "Error: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract identifier from filename
identifier=$(basename "$file" _pass.fastq.gz)

# Extract dataset name from identifier by removing the SP27_XXX_SP5_XXX_ pattern
dataset=$(echo "$identifier" | sed 's/.*SP27_[0-9]*_SP5_[0-9]*_//')

echo "Identifier: $identifier"
echo "Dataset: $dataset"

output_dir="${output_base_dir}/${identifier}"

mkdir -p "${output_dir}"

# Run SPAdes: this needs installing prior to running this code and adding to $PATH
spades.py \
 -s ${file} \
 -t 24 \
 -o "${output_dir}"

# Define the contigs file path (SPAdes output)
contigs_file="${output_dir}/contigs.fasta"

# Check if contigs file exists
if [ ! -f "$contigs_file" ]; then
    echo "Error: contigs.fasta not found at $contigs_file"
    exit 1
fi

# Rename contigs.fasta to include identifier
renamed_contigs="${output_dir}/${identifier}_contigs.fa"
mv "$contigs_file" "$renamed_contigs"

# Use seqkit replace to modify headers
# First replace: capture NODE info and reformat, then prepend identifier
seqkit replace -p '^>NODE_(\d+)_length_(\d+)_cov_([0-9.]+)' -r ">${identifier}_NODE\${1}_length\${2}_cov\${3}" "$renamed_contigs" | \
seqkit replace -p ' ' -r '_' > "${output_dir}/${identifier}.contigs.renamed.fa"
