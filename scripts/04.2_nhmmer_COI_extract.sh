#!/bin/bash
#SBATCH --job-name=extract_COIs
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
#SBATCH --array=1-169

# Check arguments
if [ $# -eq 0 ]; then
    echo "Error: No input assembly directory provided"
    echo "Usage: sbatch $0 <assembly_directory>"
    echo "Example: sbatch $0 /mnt/shared/projects/.../assembled"
    exit 1
fi

# Setup variables
gene="cox1"
assembly_dir="$1"               # comes from command line
basedir=$(dirname "$assembly_dir")
hmmdir=Lakes_meiofauna_workshop/adapters_hmms
output_base_dir="${basedir}/COIs"

echo "Assembly directory: ${assembly_dir}"
echo "Base directory: ${basedir}"
echo "HMM dir: ${hmmdir}"
echo "Output base dir: ${output_base_dir}"

# Get input fasta file for this array task
input_fasta=$(find "${assembly_dir}" -name "*contigs.renamed.fa" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$input_fasta" ]; then
    echo "Error: No input file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

sample=$(basename "${input_fasta}" .contigs.renamed.fa)
outdir="${output_base_dir}/${sample}"
hmm_file="${hmmdir}/COI.hmm"

mkdir -p "${outdir}"

# Define output file paths
output_log=${outdir}/${sample}_${gene}.log
output_tbl=${outdir}/${sample}_${gene}.tbl
output_aliscores=${outdir}/${sample}_${gene}_aliscores.tsv
output_bed=${outdir}/${sample}_${gene}.bed
output_extracted=${outdir}/${sample}_${gene}_extracted_sequences.fna
input_fasta_copied="${outdir}/$(basename "${input_fasta}")"

# Copy and clean input fasta wiht seqtk
source activate seqtk
seqtk seq -l 0 "${input_fasta}" > "${input_fasta_copied}"
sed -i 's/\r$//' "${input_fasta_copied}"

# Run nhmmer
conda deactivate && source activate hmmer
nhmmer \
 -o "${output_log}" \
 -T 10 \
 --tblout "${output_tbl}" \
 --aliscoresout "${output_aliscores}" \
 "${hmm_file}" \
 "${input_fasta_copied}"

# Parse tblout -> BED file
echo "Parsing tblout file..."
awk '
BEGIN { OFS="\t" }
!/^#/ && NF > 0 {
    target = $1
    alifrom = $7
    alito = $8
    if (alifrom <= alito) {
        start = alifrom - 1
        end = alito
        strand = "+"
    } else {
        start = alito - 1
        end = alifrom
        strand = "-"
    }
    print target, start, end, target"_"NR, "0", strand
}' "${output_tbl}" > "${output_bed}"

# Extract sequences using bedtools
conda deactivate && source activate bedtools
bedtools getfasta \
 -fi "${input_fasta_copied}" \
 -bed "${output_bed}" \
 -s \
 -name \
 -fo "${output_extracted}"

# remove copied fasta to save space 
rm -rf "${input_fasta_copied}"
