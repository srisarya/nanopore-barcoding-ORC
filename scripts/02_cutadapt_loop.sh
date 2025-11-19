#!/bin/bash
#SBATCH --job-name=demux_and_clean
#SBATCH --mem=4G
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log

# Check if input file argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <pychopped_input_file>"
    echo "Example: sbatch $0 /path/to/pychopped/pychopped_sample1.gz"
    exit 1
fi

# Define variables
infile="$1"  # Input file from command line argument

# Extract directory and filename components
indir=$(dirname "$infile")
filename=$(basename "$infile")

# Extract dataset name by removing "pychopped_" prefix and common file extensions
dataset="${filename#pychopped_}"      # Remove "pychopped_" prefix
# strip common compound extension first
dataset="${dataset%.fastq.gz}"
dataset="${dataset%.fastq}"
dataset="${dataset%.fq}"
dataset="${dataset%.gz}"
# remove trailing _pass if present
dataset="${dataset%_pass}"

# Create output directory in the parent directory of the pychopped folder
parent_dir=$(dirname "$indir")
outdir="$parent_dir/demuxed"

# relative paths for adapters
adapters_SP5=nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_forward.fa
adapters_SP27=nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_reverse_rc.fa
e_rate=0.1 # this can be changed if you like, but 0.2 is a good value for our adapter seqs where each is >=3bp different from the any other

# relative paths for primers
primers_fwd=nanopore-barcoding-ORC/adapters_hmms/amplicon_primers_forward.fa
primers_rvs=nanopore-barcoding-ORC/adapters_hmms/amplicon_primers_reverse.fa

echo "Processing: $infile"
echo "Dataset name: $dataset"
echo "Output directory: $outdir"

# Verify required files exist (input + adapters + primers)
for f in "$infile" "$adapters_SP5" "$adapters_SP27" "$primers_fwd" "$primers_rvs"; do
    if [ ! -f "$f" ]; then
        echo "Error: required file not found: $f"
        exit 1
    fi
done

# Create output directories: SP5 and SP27 are our primer names, so we have directories named after them
mkdir -p "$outdir"/SP5 "$outdir"/SP27 "$outdir"/primerless

source activate cutadapt # v4.9

# First round of demultiplexing with SP5 adapters
echo "Starting first round of demultiplexing..."
cutadapt \
 --action=trim \
 -e "$e_rate" -j 24 --rc \
 -g file:"$adapters_SP5" \
 -o "$outdir"/SP5/{name}_"$dataset".fastq.gz \
 "$infile" \
 --json="$outdir"/SP5/cutadapt_SP5_"$dataset".json

# Making the file with identifiers from SP5 output
identifier_file="$outdir"/SP27/cutadapt_SP27_identifier_"$dataset".txt
if [ ! -f "$identifier_file" ]; then
    suffix="_${dataset}.fastq.gz"
    for file1 in "$outdir"/SP5/*"$suffix"; do
        if [ -f "$file1" ]; then
            base1=$(basename "$file1")
            # remove the suffix to get the adapter name
            identifier="${base1%$suffix}"
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
     --action=trim \
     -e "$e_rate" -j 24 --rc \
     -a file:"$adapters_SP27" \
     -o "$outdir"/SP27/{name}_"$identifier"_"$dataset".fastq.gz \
     "$outdir"/SP5/"$identifier"_"$dataset".fastq.gz \
     --json="$outdir"/SP27/"$identifier"_"$dataset".json
     
done < "$identifier_file" && echo "Demultiplexing complete!"

