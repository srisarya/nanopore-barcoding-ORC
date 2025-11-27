#!/bin/bash
#SBATCH --job-name=barrnap
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-96

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No working directory provided"
    echo "Usage: sbatch $0 /path/to/workdir/primerless"
    echo "Example: sbatch $0 Lakes_day1/primerless"
    exit 1
fi

# Define directories
workdir="$1"
# Get the parent directory for output (Lakes_day1 instead of Lakes_day1/primerless)
parent_dir=$(dirname "$workdir")
input_pattern="${workdir}/*/rRNAs/cleaned_amplicon*.fasta"
output_base_dir="${parent_dir}/rRNA_genes"

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
# From: <workdir>/<sample><dataset>/rRNAs/cleaned_amplicon*.fasta
sample_path=$(dirname $(dirname "$input_file"))
identifier=$(basename "$sample_path")

echo "Identifier: ${identifier}"

# Create output directory for this sample
sample_output_dir="${output_base_dir}/${identifier}"
mkdir -p "${sample_output_dir}"

# Activate conda environment
source activate barrnap

# Create temporary directory for barrnap output
temp_dir="${sample_output_dir}/temp"
mkdir -p "${temp_dir}"

# Run barrnap
barrnap -k euk -t 1 \
    -o "${temp_dir}/${identifier}_euk.fa" \
    "${input_file}" > "${temp_dir}/${identifier}_euk.gff3"

echo "Barrnap completed for ${identifier}"

# Split the output by rRNA type (keep only 18S and 28S)
if [ -f "${temp_dir}/${identifier}_euk.fa" ]; then
    # Extract 18S sequences
    grep -A 1 "18S_rRNA" "${temp_dir}/${identifier}_euk.fa" | grep -v "^--$" > "${sample_output_dir}/${identifier}_18S.fa"
    
    # Extract 28S sequences
    grep -A 1 "28S_rRNA" "${temp_dir}/${identifier}_euk.fa" | grep -v "^--$" > "${sample_output_dir}/${identifier}_28S.fa"
    
    echo "Split rRNA sequences for ${identifier}"
    echo "18S sequences: $(grep -c "^>" "${sample_output_dir}/${identifier}_18S.fa" 2>/dev/null || echo 0)"
    echo "28S sequences: $(grep -c "^>" "${sample_output_dir}/${identifier}_28S.fa" 2>/dev/null || echo 0)"
else
    echo "Warning: No barrnap output found for ${identifier}"
fi

# Clean up temporary files
rm -rf "${temp_dir}"

echo "Processing completed for ${identifier}"