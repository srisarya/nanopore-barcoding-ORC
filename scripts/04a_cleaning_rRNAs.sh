#!/bin/bash
#SBATCH --job-name=trim_consensus_rRNA
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

start_time=$(date +%s)
start_datetime=$(date '+%Y-%m-%d %H:%M:%S')
echo "--- Primer removal started at: $start_datetime ---"
echo ""

# Check input
if [ $# -lt 1 ]; then
    echo "Error: Missing required argument"
    echo "Usage: sbatch $0 <amplicon_sorted_directory>"
    exit 1
fi

amplicon_sorted_dir="$1"

echo "Amplicon sorted directory: $amplicon_sorted_dir"

# Get *_consensus_rRNAs.fasta for this array index
consensus_file=$(find "$amplicon_sorted_dir" -type f -name "*_consensus_rRNAs.fasta" -path "*/rRNAs/*" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Processing file: $consensus_file"

if [ ! -f "$consensus_file" ]; then
    echo "Error: No consensus file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Extract identifier
identifier_dir=$(dirname "$consensus_file")  
sample_dir=$(dirname "$identifier_dir")      
identifier=$(basename "$sample_dir")         

echo "Sample identifier: $identifier"

workdir=$(dirname "$amplicon_sorted_dir") 

# Create subdirectories for round 1 and round 2
output_dir_r1="${workdir}/primerless/${identifier}/rRNAs/round1"
output_dir_r2="${workdir}/primerless/${identifier}/COIs/round2"
mkdir -p "$output_dir_r1"
mkdir -p "$output_dir_r2"

echo "Round 1 output directory: $output_dir_r1"
echo "Round 2 output directory: $output_dir_r2"
echo ""

# Output files
output_fasta_r1="${output_dir_r1}/primerless_rRNAs_${identifier}.fasta"
untrimmed_fasta_r1="${output_dir_r1}/untrimmed_rRNAs_${identifier}.fasta"
output_fasta_r2="${output_dir_r2}/recategorised_COIs_${identifier}.fasta"

# Run cutadapt
source activate cutadapt

# Round 1: trim rRNA primers and put untrimmed sequences into a separate file
echo "--- Round 1: Trimming rRNA primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -g "GCTTGTCTCAAAGATTAAGCC...ACCCGCTGAAYTTAAGCATAT" \
 -g "TTTTGGTAAGCAGAACTGGYG...CTGAACGCCTCTAAGKYRGWA" \
 --untrimmed-output="$untrimmed_fasta_r1" \
 -o "$output_fasta_r1" \
 "$consensus_file"

echo "Round 1 completed. Trimmed sequences: $output_fasta_r1"
echo "Round 1 untrimmed sequences: $untrimmed_fasta_r1"
echo ""

# Round 2: test untrimmed sequences for COI primers and recategorise if found
echo "--- Round 2: Testing untrimmed sequences for COI primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -g "TNTCNACNAAYCAYAARGAYATTGG...TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 -g "GGDRCWGGWTGAACWGTWTAYCCNCC...TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 --discard-untrimmed \
 -o "$output_fasta_r2" \
 "$untrimmed_fasta_r1"

echo "--- Primer trimming completed ---"
echo "Round 1 (rRNA) output: $output_fasta_r1"
echo "Round 2 (recategorised COI) output: $output_fasta_r2"
echo ""
echo "--- Job completed at: $(date '+%Y-%m-%d %H:%M:%S') ---"