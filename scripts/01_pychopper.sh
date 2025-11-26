#!/bin/bash
#SBATCH --job-name=pychopper
#SBATCH --mem=2G
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log

# Check if input file argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No input file provided"
    echo "Usage: sbatch $0 <input_fastq_file>"
    exit 1
fi

# variables
Qscore=10  # Set the quality score threshold
infile="$1"  # Input FASTQ file from command line argument

# Extract directory and filename components
indir=$(dirname "$infile")
filename=$(basename "$infile")
# Remove .gz extension if present, then remove .fastq/.fq extension
basename_no_ext="${filename%.gz}"
basename_no_ext="${basename_no_ext%.fastq}"
basename_no_ext="${basename_no_ext%.fq}"

# Create output directory and filenames based on input
outdir="$indir/pychopped"
outfile="$outdir/pychopped_${basename_no_ext}.gz"

# Static paths for primer sequences and config
primer_seqs=nanopore-barcoding-ORC/adapters_primers/M13_seqs_for_pychopper.fa # primer sequences plus N wildcards for variable section
config=nanopore-barcoding-ORC/adapters_primers/M13_config_for_pychopper.txt # configuration file for sequence orientation

# Make directories if not existing already
mkdir -p "$outdir"

# for my specific cluster, I need source activate my environments. Most clusters use conda activate (or mamba activate)
source activate pychopper_v2 # v2.7.10

echo "Processing: $infile"
echo "Output directory: $outdir"
echo "Output file: $outfile"

pychopper \
 -b "$primer_seqs" \
 -c "$config" \
 -k LSK114 \
 -Q "$Qscore" \
 -w "$outdir"/${basename_no_ext}_rescued.fastq \
 -u "$outdir"/${basename_no_ext}_unclass.fastq \
 -l "$outdir"/${basename_no_ext}_short.fastq \
 -S "$outdir"/${basename_no_ext}_stats.out \
 -p \
 -t 24 \
 -m edlib \
 "$infile" > "$outdir"/${basename_no_ext}_pass.fastq