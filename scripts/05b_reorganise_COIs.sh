#!/bin/bash
#SBATCH --job-name=copy_COI
#SBATCH --mem=100M
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-96

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No working directory provided"
    echo "Usage: sbatch $0 /path/to/dataset/primerless"
    echo "Example: sbatch $0 Lakes_day1/primerless"
    exit 1
fi

# Define directories
workdir="$1"
parent_dir=$(dirname "$workdir")
input_pattern="${workdir}/*/COIs/cleaned*.fasta"
output_base_dir="${parent_dir}/COI_gene"

echo "Working directory: ${workdir}"
echo "Parent directory: ${parent_dir}"
echo "Input pattern: ${input_pattern}"
echo "Output base dir: ${output_base_dir}"

# Get the filename for the current task
input_file=$(ls ${input_pattern} 2>/dev/null | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "File being processed: ${input_file}"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract sample identifier from path
# From: <workdir>/<sample>/COIs/cleaned*.fasta
sample_path=$(dirname $(dirname "$input_file"))
identifier=$(basename "$sample_path")

echo "Identifier: ${identifier}"

# Create output directory for this sample
sample_output_dir="${output_base_dir}/${identifier}"
mkdir -p "${sample_output_dir}"

# Copy the file with new name
output_file="${sample_output_dir}/${identifier}_COI.fasta"
cp "${input_file}" "${output_file}"

echo "Copied ${input_file} to ${output_file}"
echo "Copy completed for ${identifier}"