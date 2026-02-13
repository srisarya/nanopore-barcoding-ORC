#!/bin/bash
#SBATCH --job-name=cleaning_28S
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-96

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No plate directory provided"
    echo "Usage: sbatch $0 <plate_name>"
    echo "Example: sbatch $0 Lakes_day1"
    exit 1
fi

# Define directories
plate="$1"
input_pattern="${plate}/rRNA_genes/*/*.fa"

echo "Plate directory: ${plate}"
echo "Input pattern: ${input_pattern}"

# Get the filename for the current task
input_file=$(ls ${input_pattern} 2>/dev/null | grep "_28S.fa$" | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "File being processed: ${input_file}"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: No file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract sample identifier from path
# From: <plate>/rRNA_genes/<sample>/<sample>_28S.fa
sample_dir=$(dirname "$input_file")
sample=$(basename "$sample_dir")

echo "Sample: ${sample}"

# Create output directory for this sample
output_dir="${sample_dir}/28S_pybarrnap_checks"
mkdir -p "${output_dir}"

#~~~~~~~~~~~~#
# Run pybarrnap
# Activate conda environment
source activate pybarrnap

# Run pybarrnap
pybarrnap \
 --accurate \
 -k euk \
 -t 4 \
 -o "${output_dir}/${sample}_pybarrnap_raw.fa" \
 "${input_file}" > "${output_dir}/${sample}_pybarrnap_raw.gff3"

echo "pybarrnap completed for ${sample}"

#~~~~~~~~~~~~#
# Filter GFF3 to retain only 28S hits
echo "Filtering GFF3 for 28S hits only..."

awk '$3 == "rRNA" && $9 ~ /28S_rRNA/ {print}' \
    "${output_dir}/${sample}_pybarrnap_raw.gff3" \
    > "${output_dir}/${sample}_28S_only.gff3"

# Count hits
count_28S=$(grep -c "28S_rRNA" "${output_dir}/${sample}_pybarrnap_raw.gff3" || echo 0)
count_5_8S=$(grep -c "5_8S_rRNA" "${output_dir}/${sample}_pybarrnap_raw.gff3" || echo 0)
count_18S=$(grep -c "18S_rRNA" "${output_dir}/${sample}_pybarrnap_raw.gff3" || echo 0)

echo "Found ${count_28S} 28S hits, ${count_5_8S} 5.8S hits, and ${count_18S} 18S hits (5.8S and 18S will be filtered out)"

#~~~~~~~~~~~~#
# Filter FASTA to retain only 28S sequences (exclude 5.8S and 18S)
conda deactivate && source activate seqkit

if [ -f "${output_dir}/${sample}_pybarrnap_raw.fa" ]; then
    # Extract only 28S sequences (exclude 5.8S and 18S)
    seqkit grep -r -p "28S_rRNA" \
        "${output_dir}/${sample}_pybarrnap_raw.fa" | \
    seqkit grep -r -v -p "5.8S_rRNA" \
        > "${output_dir}/${sample}_28S_only.fa"
    
    # Verify output
    if [ -f "${output_dir}/${sample}_28S_only.fa" ]; then
        seq_count=$(grep -c "^>" "${output_dir}/${sample}_28S_only.fa" || echo 0)
        echo "✓ 28S-only FASTA created with ${seq_count} sequences: ${output_dir}/${sample}_28S_only.fa"
    else
        echo "WARNING: 28S-only FASTA file not created!"
    fi
    
    # Count what was found in FINAL output
    check_28S=$(seqkit grep -r -p "28S_rRNA" "${output_dir}/${sample}_28S_only.fa" 2>/dev/null | seqkit grep -r -v -p "5.8S_rRNA" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
    check_5_8S=$(seqkit grep -r -p "5.8S_rRNA" "${output_dir}/${sample}_28S_only.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
    check_18S=$(seqkit grep -r -p "18S_rRNA" "${output_dir}/${sample}_28S_only.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
    
    # Remove any whitespace/newlines from the counts
    check_5_8S=$(echo "$check_5_8S" | tr -d '\n\r' | awk '{print $1}')
    check_18S=$(echo "$check_18S" | tr -d '\n\r' | awk '{print $1}')
    
    echo "Post filtering, found: ${check_28S} 28S hits, ${check_5_8S} 5.8S hits (filtered), ${check_18S} 18S hits (filtered)"

    # Check for contamination in final output
    if [ "$check_5_8S" -gt 0 ] 2>/dev/null || [ "$check_18S" -gt 0 ] 2>/dev/null; then
        echo "!!! ERROR !!! Final FASTA contains non-28S sequences!"
        echo "  - 5.8S sequences found: ${check_5_8S}"
        echo "  - 18S sequences found: ${check_18S}"
        echo "Filtering failed - please check ${output_dir}/${sample}_28S_only.fa"
        exit 1
    else
        echo "✓ Verified: Final FASTA contains only 28S sequences"
    fi    
    # Optional: Remove raw files to save space
    # rm -f "${output_dir}/${sample}_pybarrnap_raw.fa" "${output_dir}/${sample}_pybarrnap_raw.gff3"
else
    echo "!!! ERROR !!! No pybarrnap output found for ${sample}"
    exit 1
fi


echo "Processing completed for ${sample}"
echo "Output files:"
echo "  - GFF3: ${output_dir}/${sample}_28S_only.gff3"
echo "  - FASTA: ${output_dir}/${sample}_28S_only.fa"
