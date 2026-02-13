#!/bin/bash
#SBATCH --job-name=split_28S_amplicons
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
# Activate seqkit for filtering
conda deactivate && source activate seqkit

if [ ! -f "${output_dir}/${sample}_pybarrnap_raw.fa" ]; then
    echo "!!! ERROR !!! No pybarrnap output found for ${sample}"
    exit 1
fi

#~~~~~~~~~~~~#
# Extract contig IDs (group[number]_readcount_[number]) from pybarrnap output

# Get contigs that have 18S or 5.8S hits (contaminated contigs)
grep "18S_rRNA\|5.8S_rRNA\|5_8S_rRNA" "${output_dir}/${sample}_pybarrnap_raw.fa" | \
    grep "^>" | \
    sed 's/^>//' | \
    grep -o "group[0-9]*_readcount_[0-9]*" | \
    sort | uniq > "${output_dir}/contaminated_contig_ids.txt"

# Get all contig IDs from the original input file
grep "^>" "${input_file}" | \
    sed 's/^>//' | \
    grep -o "group[0-9]*_readcount_[0-9]*" | \
    sort | uniq > "${output_dir}/all_contig_ids.txt"

# Contigs from 18S+ amplicon: have contamination
cp "${output_dir}/contaminated_contig_ids.txt" "${output_dir}/from_18S_plus_amplicon.txt"

# Contigs from 28S-only amplicon: NO contamination
comm -23 "${output_dir}/all_contig_ids.txt" "${output_dir}/contaminated_contig_ids.txt" > "${output_dir}/from_28S_only_amplicon.txt"

#~~~~~~~~~~~~#
# Extract sequences from ORIGINAL input file based on contig classification

echo "Splitting original contigs..."

# 18S+ amplicon contigs (contaminated)
if [ -s "${output_dir}/from_18S_plus_amplicon.txt" ]; then
    # Extract full contigs from original input
    seqkit grep -r -f "${output_dir}/from_18S_plus_amplicon.txt" \
        "${input_file}" \
        > "${output_dir}/${sample}_contigs_from_18S_plus_amplicon.fa"
    
    count_18S_contigs=$(grep -c "^>" "${output_dir}/${sample}_contigs_from_18S_plus_amplicon.fa" 2>/dev/null || echo 0)
    echo "✓ Found ${count_18S_contigs} contigs from 18S+ amplicon (had 18S/5.8S contamination)"
else
    touch "${output_dir}/${sample}_contigs_from_18S_plus_amplicon.fa"
    echo "No contaminated contigs found"
fi

# 28S-only amplicon contigs (clean)
if [ -s "${output_dir}/from_28S_only_amplicon.txt" ]; then
    # Extract full contigs from original input
    seqkit grep -r -f "${output_dir}/from_28S_only_amplicon.txt" \
        "${input_file}" \
        > "${output_dir}/${sample}_contigs_from_28S_only_amplicon.fa"
    
    count_28S_contigs=$(grep -c "^>" "${output_dir}/${sample}_contigs_from_28S_only_amplicon.fa" 2>/dev/null || echo 0)
    echo "✓ Found ${count_28S_contigs} contigs from 28S-only amplicon (no contamination)"
else
    touch "${output_dir}/${sample}_contigs_from_28S_only_amplicon.fa"
    echo "No clean 28S-only contigs found"
fi

#~~~~~~~~~~~~#
# Extract corresponding 28S regions from pybarrnap output for each amplicon type

# For 18S+ amplicon: get only the 28S portions
if [ -s "${output_dir}/from_18S_plus_amplicon.txt" ]; then
    seqkit grep -r -f "${output_dir}/from_18S_plus_amplicon.txt" \
        "${output_dir}/${sample}_pybarrnap_raw.fa" | \
    seqkit grep -r -p "28S_rRNA" \
        > "${output_dir}/${sample}_28S_from_18S_plus_amplicon.fa"
    
    # Extract GFF3 for these contigs - 28S only
    grep -f "${output_dir}/from_18S_plus_amplicon.txt" \
        "${output_dir}/${sample}_pybarrnap_raw.gff3" | \
    awk '$9 ~ /28S_rRNA/ {print}' \
        > "${output_dir}/${sample}_28S_from_18S_plus_amplicon.gff3"
    
    count_28S_18S=$(grep -c "^>" "${output_dir}/${sample}_28S_from_18S_plus_amplicon.fa" 2>/dev/null || echo 0)
    echo "✓ Extracted ${count_28S_18S} 28S regions from 18S+ amplicon contigs"
fi

# For 28S-only amplicon: get the 28S portions
if [ -s "${output_dir}/from_28S_only_amplicon.txt" ]; then
    seqkit grep -r -f "${output_dir}/from_28S_only_amplicon.txt" \
        "${output_dir}/${sample}_pybarrnap_raw.fa" | \
    seqkit grep -r -p "28S_rRNA" \
        > "${output_dir}/${sample}_28S_from_28S_only_amplicon.fa"
    
    # Extract GFF3 for these contigs - 28S only
    grep -f "${output_dir}/from_28S_only_amplicon.txt" \
        "${output_dir}/${sample}_pybarrnap_raw.gff3" | \
    awk '$9 ~ /28S_rRNA/ {print}' \
        > "${output_dir}/${sample}_28S_from_28S_only_amplicon.gff3"
    
    count_28S_only=$(grep -c "^>" "${output_dir}/${sample}_28S_from_28S_only_amplicon.fa" 2>/dev/null || echo 0)
    echo "✓ Extracted ${count_28S_only} 28S regions from 28S-only amplicon contigs"
fi

#~~~~~~~~~~~~#
# Verify no contamination in the 28S-extracted outputs
check_5_8S_in_18S=$(seqkit grep -r -p "5.8S_rRNA\|5_8S_rRNA" "${output_dir}/${sample}_28S_from_18S_plus_amplicon.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
check_18S_in_18S=$(seqkit grep -r -p "18S_rRNA" "${output_dir}/${sample}_28S_from_18S_plus_amplicon.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
check_5_8S_in_28S=$(seqkit grep -r -p "5.8S_rRNA\|5_8S_rRNA" "${output_dir}/${sample}_28S_from_28S_only_amplicon.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)
check_18S_in_28S=$(seqkit grep -r -p "18S_rRNA" "${output_dir}/${sample}_28S_from_28S_only_amplicon.fa" 2>/dev/null | grep -c "^>" 2>/dev/null || echo 0)

# Clean up counts
check_5_8S_in_18S=$(echo "$check_5_8S_in_18S" | tr -d '\n\r' | awk '{print $1}')
check_18S_in_18S=$(echo "$check_18S_in_18S" | tr -d '\n\r' | awk '{print $1}')
check_5_8S_in_28S=$(echo "$check_5_8S_in_28S" | tr -d '\n\r' | awk '{print $1}')
check_18S_in_28S=$(echo "$check_18S_in_28S" | tr -d '\n\r' | awk '{print $1}')

echo "Contamination check on 28S-extracted outputs:"
echo "  18S+ amplicon 28S file: ${check_5_8S_in_18S} 5.8S hits, ${check_18S_in_18S} 18S hits"
echo "  28S-only amplicon 28S file: ${check_5_8S_in_28S} 5.8S hits, ${check_18S_in_28S} 18S hits"

# Error if ANY 28S output file has contamination
error=0
if [ "$check_5_8S_in_18S" -gt 0 ] 2>/dev/null || [ "$check_18S_in_18S" -gt 0 ] 2>/dev/null; then
    echo "!!! ERROR !!! 18S+ amplicon 28S output contains non-28S sequences!"
    error=1
fi

if [ "$check_5_8S_in_28S" -gt 0 ] 2>/dev/null || [ "$check_18S_in_28S" -gt 0 ] 2>/dev/null; then
    echo "!!! ERROR !!! 28S-only amplicon 28S output contains non-28S sequences!"
    error=1
fi

if [ $error -eq 1 ]; then
    exit 1
else
    echo "✓ Verified: Both 28S output files contain only 28S sequences"
fi

# Cleanup intermediate files
#rm -f "${output_dir}/contaminated_contig_ids.txt" \
#      "${output_dir}/all_contig_ids.txt" \
#      "${output_dir}/from_18S_plus_amplicon.txt" \
#      "${output_dir}/from_28S_only_amplicon.txt"

echo "Processing completed for ${sample}"
echo "Output files:"
echo "  Original contigs from 18S+ amplicon: ${output_dir}/${sample}_contigs_from_18S_plus_amplicon.fa"
echo "  Original contigs from 28S-only amplicon: ${output_dir}/${sample}_contigs_from_28S_only_amplicon.fa"
echo "  28S regions from 18S+ amplicon: ${output_dir}/${sample}_28S_from_18S_plus_amplicon.fa + .gff3"
echo "  28S regions from 28S-only amplicon: ${output_dir}/${sample}_28S_from_28S_only_amplicon.fa + .gff3"
