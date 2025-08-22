#!/bin/bash
#SBATCH --job-name=barrnap
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input directory provided"
    echo "Usage: sbatch $0 <assembled_directory>"
    echo "Example: sbatch $0 /path/to/assembled"
    exit 1
fi

# Define input and output directories from command line argument
input_base_dir="$1"  # assembled directory from command line
basedir=$(dirname "$input_base_dir")  # Go up one level to get base directory
output_base_dir="${basedir}/rRNAs"
split_files="${basedir}/split_files"
combined_rRNA="${basedir}/combined_rRNA"

echo "Input base dir: ${input_base_dir}"
echo "Base directory: ${basedir}"
echo "Output base dir: ${output_base_dir}"
echo "Split files dir: ${split_files}"
echo "Combined rRNA dir: ${combined_rRNA}"

# Create output directories
mkdir -p "${output_base_dir}" "${split_files}" "${combined_rRNA}"

# Get the filename for the current task
input_file=$(ls "${input_base_dir}"/*.contigs.renamed.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "File being processed: ${input_file}"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract identifier from filename
identifier=$(basename "$input_file" .contigs.renamed.fa)

echo "Identifier: ${identifier}"

# Activate conda environment
source activate barrnap

# Run barrnap
barrnap -k euk -t 1 \
    -o "${output_base_dir}/${identifier}_euk.fa" \
    "${input_file}" > "${output_base_dir}/${identifier}_euk.gff3"

echo "Barrnap completed for ${identifier}"
