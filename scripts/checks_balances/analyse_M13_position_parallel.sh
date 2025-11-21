#!/bin/bash

# Activate environment
source activate python_stuff

# Define the Python script path
SCRIPT="nanopore-barcoding-ORC/scripts/checks_balances/adapter_position.py"

# Run variable adapter analysis in parallel
parallel -j 3 --line-buffer \
  "python ${SCRIPT} \
   nanopore-barcoding-ORC/adapters_hmms/M13_segment_only.fa \
   day{1}_fastq_pass/demuxed/SP27 \
   day{1}_M13positions_variable.csv" \
  ::: 1 2 3

