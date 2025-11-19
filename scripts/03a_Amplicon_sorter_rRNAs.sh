#!/bin/bash
#SBATCH --job-name=AS_rRNA_array
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

# Record start time
start_time=$(date +%s)
start_datetime=$(date '+%Y-%m-%d %H:%M:%S')
echo "=== Script started at: $start_datetime ==="
echo ""

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
output_base_dir="${basedir}/amplicon_sorted"

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

# Extract identifier_nopass from filename
identifier_nosuff=$(basename "$file" fastq.gz)
identifier_nopass=$(basename "$identifier_nosuff" _pass.)

# Extract dataset name from identifier_nopass by removing the SP27_XXX_SP5_XXX_ pattern
dataset=$(echo "$identifier_nopass" | sed 's/.*SP27_[0-9]*_SP5_[0-9]*_//')

echo "identifier_nopass: $identifier_nopass"
echo "Dataset: $dataset"

outfolder_rRNAs="${output_base_dir}/${identifier_nopass}/rRNAs"

mkdir -p "${outfolder_rRNAs}"

# Run Amplicon_sorter: this is both a conda environment and a script
echo ""
echo "--- Starting Amplicon_sorter at $(date '+%H:%M:%S') ---"
source activate amplicon_sorter
python3 nanopore-barcoding-ORC/scripts/for_AS/amplicon_sorter.py \
 -i ${file} \
 -o ${outfolder_rRNAs} \
 -ar \
 -np 4 \
 -min 3000
echo "--- Amplicon_sorter completed at $(date '+%H:%M:%S') ---"

# Verify the expected output consensus file exists; exit with non-zero code if missing
consensus_file="${outfolder_rRNAs}/consensusfile.fasta"
if [ ! -f "${consensus_file}" ]; then
    echo "Error: output file '${consensus_file}' not made!! See $outfolder_rRNAs/$identifier_nosuff/results.txt for details'"
    exit 2
fi

# Use seqkit replace to modify headers; want to replace parantheses with underscores and explanations
echo ""
echo "--- Starting header renaming at $(date '+%H:%M:%S') ---"
conda deactivate && source activate seqkit

# Step 1: Replace parentheses with _readcount_
seqkit replace -p '\((\d+)\)$' -r '_readcount_$1' "${consensus_file}" > "${consensus_file}.tmp1"

# Step 2: Dynamically replace _X_Y patterns with sequential _groupN
# This uses a counter that increments for each header found
awk '
BEGIN {counter = 1}
/^>/ {
    # Match the pattern _X_Y_readcount and replace with _groupN_readcount
    if (match($0, /_[0-9]+_[0-9]+_readcount/)) {
        sub(/_[0-9]+_[0-9]+_readcount/, "_group" counter "_readcount")
        counter++
    }
}
{print}
' "${consensus_file}.tmp1" > "${identifier_nopass}_consensus_rRNAs.fasta"

# Clean up temporary files and keep renamed files
rm "${consensus_file}.tmp1"
rm "${consensus_file}"

echo "--- Header renaming completed at $(date '+%H:%M:%S') ---"

# Calculate and display end time and duration
end_time=$(date +%s)
end_datetime=$(date '+%Y-%m-%d %H:%M:%S')
duration=$((end_time - start_time))

# Convert duration to hours, minutes, seconds
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo "=== Script completed at: $end_datetime ==="
echo "=== Total runtime: ${hours}h ${minutes}m ${seconds}s (${duration} seconds) ==="