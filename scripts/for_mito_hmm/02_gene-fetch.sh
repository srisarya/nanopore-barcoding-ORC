#!/bin/bash
#SBATCH --job-name=gene_fetch_COI
#SBATCH --mem=2G
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-604

# directories
basedir=Lakes_meiofauna_workshop/adapters_hmms/COI_seqs

# variables
taxon=metazoa
level=order
gene=COI  # You had ${gene} which was undefined
taxid_file=${basedir}/${taxon}_orders.txt
type=nucleotide
max_sequences=5

email=bty636@qmul.ac.uk
api_key=3e5fa6db78edc8fafad2ba363cfa0fb63e0a

# Get the taxid for this array task (FIXED: use SLURM_ARRAY_TASK_ID instead of SGE_TASK_ID)
this_taxid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${taxid_file}")

echo "Processing taxid: ${this_taxid}"

# make a subdirectory for each taxid
taxid_outdir="${basedir}/${this_taxid}"
mkdir -p "${taxid_outdir}"

source activate gene-fetch
gene-fetch \
    --gene "${gene}" \
    -s "${this_taxid}" \
    --type "${type}" \
    --out "${taxid_outdir}" \
    --max-sequences "${max_sequences}" \
    --email "${email}" \
    --api-key "${api_key}"
