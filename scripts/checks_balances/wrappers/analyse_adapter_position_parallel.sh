#!/bin/bash

# Activate environment
source activate python_stuff

# Define the Python script path
SCRIPT="nanopore-barcoding-ORC/scripts/checks_balances/adapter_position.py"

# Run full adapter analysis in parallel
parallel -j 3 --line-buffer \
  "python ${SCRIPT} \
   nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_all.fa \
   day{1}_fastq_pass/demuxed_keepadapters/SP27 \
   day{1}_adapterpositions_full.csv" \
  ::: 1 2 3

# Run variable adapter analysis in parallel
parallel -j 3 --line-buffer \
  "python ${SCRIPT} \
   nanopore-barcoding-ORC/adapters_hmms/M13_variable_indices_all.fa \
   day{1}_fastq_pass/demuxed_keepadapters/SP27 \
   day{1}_adapterpositions_variable.csv" \
  ::: 1 2 3

echo "All adapter position analyses complete!"
