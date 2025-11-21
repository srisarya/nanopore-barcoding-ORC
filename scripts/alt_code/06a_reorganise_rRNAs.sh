#!/bin/bash
#SBATCH --job-name=rRNA_reorg
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log

# Check if input directory argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input directory provided"
    echo "Usage: sbatch $0 <input directory>"
    echo "Example: sbatch $0 /path/to/rRNAs"
    exit 1
fi

# Define input and output directories from command line argument
input_dir="$1"  # assembled directory from command line
split_files="${input_dir}/split_files"
combined_rRNA="${input_dir}/combined_rRNA"

# Create output directory if it doesn't exist
mkdir -p "${split_files}" "${combined_rRNA}"

source activate UCSC_utils

# Process each multifasta file for eukaryotic rRNAs
for file in "${input_dir}"/*.fa; do
    # Extract sample name from filename (everything before the first period)
    filename=$(basename "$file")
    sample_name="${filename%%.*}"
    echo "Processing sample: $sample_name"

    # Create sample-specific directory
    sample_dir="${split_files}/${sample_name}"
    mkdir -p "$sample_dir"
    # faSplit byname splits the file into separate files named based on their headers
    faSplit byname "$file" "${sample_dir}/"


    #change conda env
    conda deactivate && source activate seqkit

    # Rename headers in split files to include sample name using seqkit
    for split_file in "${sample_dir}"/*.fa; do
        if [[ -f "$split_file" ]]; then
            # Get the base name without extension for the rRNA type
            rRNA_file=$(basename "$split_file" .fa)

            # Use seqkit to rename headers: extract rRNA type, length, cov, and strand info
            # Pattern matches: >18S_rRNA::NODE_1_length_3896_cov_104.111435:35-1825(+)
            # Captures: rRNA type, length, coverage, strand
            seqkit replace -p '^(\w+)_rRNA::NODE_\d+_length_(\d+)_cov_([\d.]+):\d+-\d+\(([+-])\)' \
                          -r "${sample_name}_\${1}_length_\${2}_cov_\${3}_\${4}" "$split_file" > "${split_file}.tmp"
            mv "${split_file}.tmp" "$split_file"

            echo "Renamed headers in $split_file"
        fi
    done
    # swiching back to UCSC_utils for other operations
    conda deactivate && source activate UCSC_utils
done

# # Group rRNAs per sample (instead of one massive file per rRNA type)
# rrna_types=("18S" "28S" "5_8S")

# # Echo out the rRNA types being processed
# echo "rRNA types to process: ${rrna_types[@]}"

# # Process each sample directory
# for sample_dir in "${split_files}"/*/; do
#     if [[ -d "$sample_dir" ]]; then
#         sample_name=$(basename "$sample_dir")
#         echo "Processing rRNAs for sample: $sample_name"
        
#         # Loop through each rRNA type for this sample
#         for rRNA in "${rrna_types[@]}"; do
#             # Find all rRNA files of this type for this sample
#             rRNA_files=(${sample_dir}${rRNA}*.fa)
            
#             if [[ -e "${rRNA_files[0]}" ]]; then
#                 # Concatenate all rRNA sequences of this type for this sample
#                 cat "${rRNA_files[@]}" > "${combined_rRNA}/${sample_name}_${rRNA}.fasta"
#                 echo "Created ${combined_rRNA}/${sample_name}_${rRNA}.fasta with ${#rRNA_files[@]} sequences"
#             else
#                 echo "No ${rRNA} sequences found for sample ${sample_name}"
#             fi
#         done
#     fi
# done

# echo "Processing complete!"
