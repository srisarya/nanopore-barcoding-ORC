#!/bin/bash
#SBATCH --job-name=adapter_search
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=4G

# Activate environment
source activate python_stuff

# Define the Python script path
SCRIPT="nanopore-barcoding-ORC/scripts/checks_balances/adapter_search.py"

# Define adapters
FULL_ADAPTER="nanopore-barcoding-ORC/adapters_hmms/M13_amplicon_indices_all.fa"
VARIABLE_ADAPTER="nanopore-barcoding-ORC/adapters_hmms/M13_variable_indices_all.fa"

# Run all 6 adapter searches in parallel (3 days × 2 adapters)
parallel -j 6 --line-buffer \
  "python ${SCRIPT} \
   {1} \
   day{3}_fastq_pass/demuxed_keepadapters/SP27 \
   day{3}_{2}_presence.csv" \
  ::: "${FULL_ADAPTER}" "${VARIABLE_ADAPTER}" \
  :::+ full variable \
  ::: 1 2 3

echo "All adapter search analyses complete!"
