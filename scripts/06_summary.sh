#!/bin/bash
#SBATCH --job-name=summary
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log

source activate R

for i in * ; do
    if [ -d "$i" ]; then
        base_dir=$i
        amplicon="COI"
        Rscript nanopore-barcoding-ORC/scripts/auxiliary_code/amplicon_summary.R "$base_dir" "$amplicon"
    fi
done       

for i in * ; do
    if [ -d "$i" ]; then
        base_dir=$i
        amplicon="rRNA"
        gene="18S"
        Rscript nanopore-barcoding-ORC/scripts/auxiliary_code/amplicon_summary.R "$base_dir" "$amplicon" "$gene"
    fi
done  

for i in * ; do
    if [ -d "$i" ]; then
        base_dir=$i
        amplicon="rRNA"
        gene="28S"
        Rscript nanopore-barcoding-ORC/scripts/auxiliary_code/amplicon_summary.R "$base_dir" "$amplicon" "$gene"
    fi
done  


