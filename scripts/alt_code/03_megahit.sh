#!/bin/bash
#SBATCH --job-name=megahit_array
#SBATCH --mem=2G
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-6

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input directory provided"
    echo "Usage: sbatch $0 <demuxed_SP27_directory>"
    echo "Example: sbatch $0 /path/to/demuxed/SP27"
    exit 1
fi

# Activate the environment
source activate megahit # version 2.0.1

# Define input and output directories from command line argument
input_dir="$1"  # SP27 directory from command line
basedir=$(dirname $(dirname "$input_dir"))  # Go up two levels from SP27 to get base directory
output_base_dir="${basedir}/assembled_megahit"

echo "Input directory: $input_dir"
echo "Base directory: $basedir"
echo "Output base directory: $output_base_dir"

# make output basedir : megahit will complain if this does not exist
mkdir -p "${output_base_dir}"

# Get the filename for the current task
file=$(ls "${input_dir}"/*.fastq.gz | perl -ne 'print if /SP27_(\d+)_SP5_\1_/' | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "file being processed:" ${file}

# Check if file exists
if [ ! -f "$file" ]; then
    echo "Error: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract identifier from filename
identifier=$(basename "$file" .fastq.gz)

# Extract dataset name from identifier by removing the SP27_XXX_SP5_XXX_ pattern
dataset=$(echo "$identifier" | sed 's/.*SP27_[0-9]*_SP5_[0-9]*_//')

echo "Identifier: $identifier"
echo "Dataset: $dataset"

# the variable for the output dir. DO NOT MAKE THIS DIR, MEGAHIT WILL DO THIS, AND IF YOU MAKE IT BEFOREHAND IT WILL COMPLAIN
output_dir="${output_base_dir}/${identifier}_${SLURM_JOB_ID}"

# Run megahit
megahit --read "${file}" -t 24 --out-dir "${output_dir}" --out-prefix "${identifier}"

conda deactivate && source activate seqkit

# Extract the SP5_XXX-SP27_XXX pattern from the identifier
sample_id=$(echo "$identifier" | grep -o 'SP5_[0-9]*-SP27_[0-9]*')

# Define the contigs file path
contigs_file="${output_dir}/${identifier}.contigs.fa"
echo "Contigs file: ${contigs_file}"

# Extract basename by removing the dataset suffix
basename=$(echo "$identifier" | sed "s/_${dataset}$//")
echo "Basename: ${basename}"

# Use seqkit replace to modify headers
# Replace everything before 'flag' in the header with the basename
seqkit replace -p '^.*?(flag.*)' -r "${basename}_\${1}" "$contigs_file" | \
seqkit replace -p ' ' -r '_' > "${output_dir}/${basename}.contigs.renamed.fa"

echo "Assembly and renaming complete for ${identifier}!"
