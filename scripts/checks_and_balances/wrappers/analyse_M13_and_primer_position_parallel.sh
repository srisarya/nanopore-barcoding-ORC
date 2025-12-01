#!/bin/bash
#SBATCH --job-name=positions
#SBATCH --mem=1G
#SBATCH --cpus-per-task=3
#SBATCH --time=00:15:00

# Activate environment
source activate python_stuff

# Define the Python script path
SCRIPT="nanopore-barcoding-ORC/scripts/checks_balances/adapter_position.py"

# Run variable adapter analysis in parallel
#parallel -j 3 --line-buffer \
#  "python ${SCRIPT} \
#   nanopore-barcoding-ORC/adapters_hmms/M13_segment_only.fa \
#   day{1}_fastq_pass/demuxed/SP27 \
#   day{1}_M13positions_variable.csv" \
#  ::: 1 2 3

parallel -j 3 --line-buffer \
  "python ${SCRIPT} \
   nanopore-barcoding-ORC/adapters_hmms/amplicon_primer_sequences.fa \
   day{1}_fastq_pass/demuxed/SP27 \
   day{1}_primer_positions_variable.csv" \
  ::: 1 2 3
