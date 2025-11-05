#!/bin/bash
#SBATCH --job-name=adapter_position
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G

# Activate environment
source activate python_stuff

# Define the Python script path
SCRIPT="nanopore-barcoding-ORC/scripts/checks_balances/adapter_search.py"

# Define adapters
FULL_ADAPTER="nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_all.fa"
VARIABLE_ADAPTER="nanopore-barcoding-ORC/adapters_hmms/M13_variable_indices_all.fa"

# Run all 8 adapter position analyses in parallel
parallel -j 8 --line-buffer \
  "python ${SCRIPT} \
   {1} \
   day1_fastq_pass/{3}/SP27 \
   day1_{3}_{2}_presence.csv" \
  ::: "${FULL_ADAPTER}" "${VARIABLE_ADAPTER}" \
  :::+ full variable \
  ::: demuxed_keepadapters demuxed_keepadapters_lowererate demuxed_keepadapters_norescue_lowererate demuxed_keepadapters_norescue

echo "All adapter position analyses complete!"
