#!/bin/bash
#SBATCH --job-name=cutadapt_variables
#SBATCH --mem=4G
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log


# Check if input file argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <pychopped_input_file> <output_directory> <error rate>"
    exit 1
fi

# Define variables
infile="$1"  # Input file from command line argument
outdir="$2" # Output directory from command line argument
e_rate="$3" # variable error rate from command line

# Extract directory and filename components
indir=$(dirname "$infile")
filename=$(basename "$infile")

# Extract dataset name by removing "pychopped_" prefix and ".gz" suffix
dataset_intermediate="${filename#pychopped_}"  # Remove "pychopped_" prefix
dataset="${dataset_intermediate%.gz}"          # Remove ".gz" suffix

# Create output directory in the parent directory of the pychopped folder
parent_dir=$(dirname "$indir")

# relative paths for adapters
adapters_SP5=nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_forward.fa
adapters_SP27=nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_reverse_rc.fa

echo "Processing: $infile"
echo "Dataset name: $dataset"
echo "Output directory: $outdir"

# Create output directories: SP5 and SP27 are our primer names, so we have directories named after them
mkdir -p "$outdir"/SP5 "$outdir"/SP27

source activate cutadapt # v4.9

# First round of demultiplexing with SP5 adapters
echo "Starting first round of demultiplexing..."
cutadapt \
 --action=none \
 -e "$e_rate" -j 24 \
 -g file:"$adapters_SP5" \
 -o "$outdir"/SP5/{name}_"$dataset".fastq.gz \
 "$infile" \
 --json="$outdir"/SP5/cutadapt_SP5_"$dataset".json

# Making the file with identifiers from SP5 output
identifier_file="$outdir"/SP27/cutadapt_SP27_identifier_"$dataset".txt
if [ ! -f "$identifier_file" ]; then
    for file1 in "$outdir"/SP5/*_"$dataset".fastq.gz; do
        if [ -f "$file1" ]; then
            identifier=$(basename "$file1" _"$dataset".fastq.gz)
            echo "$identifier" >> "$identifier_file"
        fi
    done
fi

# Check if identifier file exists and has content
if [ ! -s "$identifier_file" ]; then
    echo "Error: No identifiers found or identifier file is empty"
    exit 1
fi

# Loop through each identifier and run second round of cutadapt, looping through each demuxed SP5 and splitting it by SP27
echo "Starting second round of demultiplexing..."
while IFS= read -r identifier; do
    echo "Processing: $identifier"
    # Second round of demultiplexing with SP27 adapters
    cutadapt \
     --action=none \
     -e "$e_rate" -j 24 \
     -a file:"$adapters_SP27" \
     -o "$outdir"/SP27/{name}_"$identifier"_"$dataset".fastq.gz \
     "$outdir"/SP5/"$identifier"_"$dataset".fastq.gz \
     --json="$outdir"/SP27/"$identifier"_"$dataset".json
done < "$identifier_file"

echo "Demultiplexing complete!"
