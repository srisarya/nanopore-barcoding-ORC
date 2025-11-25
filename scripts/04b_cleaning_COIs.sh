#!/bin/bash
#SBATCH --job-name=trim_consensus_COI
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

# Get *_consensus_COIs.fasta for this array index
consensus_file=$(find "$amplicon_sorted_dir" -type f -name "*_consensus_COIs.fasta" -path "*/COIs/*" | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")

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
output_dir_round1="${workdir}/primerless/${identifier}/COIs/COIs_300_900bp"
output_dir_round2="${workdir}/primerless/${identifier}/COIs/rRNA_300_900bp"
mkdir -p "$output_dir_round1"
mkdir -p "$output_dir_round2"

echo "Round 1 output directory: $output_dir_round1"
echo "Round 2 output directory: $output_dir_round2"
echo ""

# Output files
primerless_fasta_round1="${output_dir_round1}/primerless_COIs_${identifier}.fasta" # the correct output
untrimmed_fasta_round1="${output_dir_round1}/untrimmed_COIs_${identifier}.fasta" # the sequences that didn't match COI primers
output_fasta_round1="${output_dir_round1}/nr_COIs_${identifier}.fasta" # clustered COIs
primerless_fasta_round2="${output_dir_round2}/rRNA_ish_${identifier}.fasta" # the sequences that matched rRNA primers

# Run cutadapt
source activate cutadapt

# Round 1: trim COI primers
echo "--- Round 1: Trimming COI primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -g "TNTCNACNAAYCAYAARGAYATTGG...TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 -g "GGDRCWGGWTGAACWGTWTAYCCNCC...TGRTTYTTYGGNCAYCCNGNRGTNTA" \
 --untrimmed-output="$untrimmed_fasta_round1" \
 -o "$primerless_fasta_round1" \
 "$consensus_file"

echo "Round 1 completed. Trimmed sequences: $primerless_fasta_round1"
echo "Round 1 untrimmed sequences: $untrimmed_fasta_round1"
echo ""

# Round 2: test untrimmed sequences for rRNA primers and recategorise if found
echo "--- Round 2: Testing untrimmed sequences for rRNA primers ---"
cutadapt \
 -j "$SLURM_CPUS_PER_TASK" \
 -e 0.4 \
 -g "GCTTGTCTCAAAGATTAAGCC" \
 -a "ACCCGCTGAAYTTAAGCATAT" \
 -g "TTTTGGTAAGCAGAACTGGYG" \
 -a "CTGAACGCCTCTAAGKYRGWA" \
 --discard-untrimmed \
 -o "$primerless_fasta_round2" \
 "$untrimmed_fasta_round1"

echo "Round 2 completed. Recategorised sequences: $primerless_fasta_round2"
echo ""

# Cluster hits for first round only
conda deactivate && source activate cd-hit

# Cluster round 1 COIs
if [ -f "$primerless_fasta_round1" ] && [ -s "$primerless_fasta_round1" ]; then
    echo "--- Clustering Round 1 COIs ---"
    cd-hit-est \
     -i "$primerless_fasta_round1" \
     -o "$output_fasta_round1" \
     -T "$SLURM_CPUS_PER_TASK"
    echo "Round 1 clustering completed: $output_fasta_round1"
    echo "Round 1 cluster info: ${output_fasta_round1}.clstr"
else
    echo "No sequences in Round 1 to cluster"
fi

echo ""
echo "--- All processing completed ---"
echo "Final outputs:"
echo "  Round 1 primer trimming COIs: $primerless_fasta_round1"
echo "  Round 1 non-redundant COIs: $output_fasta_round1"
echo "  Round 2 primer trimming rRNAs: $primerless_fasta_round2"
echo ""*
echo "--- Job completed at: $(date '+%Y-%m-%d %H:%M:%S') ---"