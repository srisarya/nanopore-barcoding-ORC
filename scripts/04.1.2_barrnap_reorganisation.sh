#!/bin/bash
#SBATCH --job-name=barrnap
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log

basedir=/mnt/shared/projects/nhm/sarya/amplicons/07_08_25_DK_BW_CEL_6MOD_AMP_CO1_RRNA
input_base_dir="${basedir}/assembled"
output_base_dir="${basedir}/rRNAs"
split_files=${basedir}/split_files    # Directory for split fasta files
combined_rRNA=${basedir}/combined_rRNA # Directory for combined rRNA files

# Create output directory if it doesn't exist
mkdir -p  "${split_files}" "${combined_rRNA}"

source activate UCSC_utils

# Process each multifasta file for eukaryotic rRNAs
for file in ${output_base_dir}/*.fa; do
    # faSplit byname splits the file into separate files named based on their headers
    faSplit byname "$file" ${split_files}/
done

# Concatenate the longest rRNA per type per sample
rrna_types=("18S" "28S" "5_8S")

# Echo out the rRNA types being processed
echo "rRNA types to process: ${rrna_types[@]}"

# Loop through each rRNA type
for rRNA in "${rrna_types[@]}"; do
    echo "Processing ${rRNA}_rRNA"

    # Concatenate only eukaryotic rRNA files
    find ${split_files} -type f -name "${rRNA}*" -print0 | xargs -0 cat > "${combined_rRNA}/${rRNA}.fasta"
    echo "Concatenated $(find ${split_files} -type f -name "${rRNA}*" | wc -l) files for ${rRNA}"

done
