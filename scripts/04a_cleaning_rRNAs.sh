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
output_dir_round1="${workdir}/primerless/${identifier}/rRNAs/rRNAs_3kbplus"
output_dir_round2="${workdir}/primerless/${identifier}/rRNAs/COIs_3kbplus"
mkdir -p "$output_dir_round1"
mkdir -p "$output_dir_round2"

echo "Round 1 output directory: $output_dir_round1"
echo "Round 2 output directory: $output_dir_round2"
echo ""

# Output files
primerless_fasta_round1="${output_dir_round1}/primerless_rRNAs_${identifier}.fasta" # the correct output
untrimmed_fasta_round1="${output_dir_round1}/untrimmed_rRNAs_${identifier}.fasta" # the sequences that didn't match COI primers
primerless_fasta_round2="${output_dir_round2}/COI_ish_${identifier}.fasta" # the sequences that matched rRNA primers

# Run cutadapt
source activate cutadapt

# Round 1: trim rRNA primers and put untrimmed sequences into a separate file
echo "--- Round 1: Trimming rRNA primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -g "GCTTGTCTCAAAGATTAAGCC...ACCCGCTGAAYTTAAGCATAT" \
 -g "TTTTGGTAAGCAGAACTGGYG...CTGAACGCCTCTAAGKYRGWA" \
 --untrimmed-output="$untrimmed_fasta_round1" \
 -o "$primerless_fasta_round1" \
 "$consensus_file"

echo "Round 1 completed. Trimmed sequences: $primerless_fasta_round1"
echo "Round 1 untrimmed sequences: $untrimmed_fasta_round1"
echo ""

# Round 2: test untrimmed sequences for COI primers and recategorise if found
echo "--- Round 2: Testing untrimmed sequences for COI primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -e 0.4 \
 -g "TNTCNACNAAYCAYAARGAYATTGG" \
 -a "TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 -g "GGDRCWGGWTGAACWGTWTAYCCNCC" \
 -a "TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 --discard-untrimmed \
 -o "$primerless_fasta_round2" \
 "$untrimmed_fasta_round1"

echo ""
echo "--- All processing completed ---"
echo "Final outputs:"
echo "  Round 1 primer trimming rRNAs: $primerless_fasta_round1"
echo "  Round 2 primer trimming COIs: $primerless_fasta_round2"
echo ""*
echo "--- Job completed at: $(date '+%Y-%m-%d %H:%M:%S') ---"