#!/bin/bash
#SBATCH --job-name=trim_consensus_array
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

# Primer files
primers_fwd=nanopore-barcoding-ORC/adapters_hmms/amplicon_primers_forward.fa
primers_rvs=nanopore-barcoding-ORC/adapters_hmms/amplicon_primers_reverse.fa
e_rate="0.1"

echo "Amplicon sorted directory: $amplicon_sorted_dir"
echo ""

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
output_dir="${workdir}/primerless/${identifier}"
mkdir -p "$output_dir"

echo "Output directory: $output_dir"
echo ""

output_fasta="${output_dir}/primerless_rRNAs_${identifier}.fasta"
output_json="${output_dir}/primerless_rRNAs_${identifier}.json"

# Run cutadapt
source activate cutadapt

cutadapt \
 --action=trim \
 -e "$e_rate" \
 -j "$SLURM_CPUS_PER_TASK" \
 -g file:"$primers_fwd" \
 -a file:"$primers_rvs" \
 -o "$output_fasta" \
 "$consensus_file" \
 --json="$output_json"

echo "--- Primer trimming completed ---"
echo "Output saved to: $output_fasta"
