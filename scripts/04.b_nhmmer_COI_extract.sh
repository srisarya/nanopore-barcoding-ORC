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
hmmdir=nanopore-barcoding-ORC/adapters_hmms
output_base_dir="${basedir}/COIs"

echo "Assembly directory: ${assembly_dir}"
echo "Base directory: ${basedir}"
echo "HMM dir: ${hmmdir}"
echo "Output base dir: ${output_base_dir}"

# Get input fasta file for this array task - USING SAME METHOD AS BARRNAP SCRIPT
input_fasta=$(ls "${assembly_dir}"/*.contigs.renamed.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Improved error checking
if [ -z "$input_fasta" ]; then
    echo "Error: No input file found for array task ${SLURM_ARRAY_TASK_ID}"
    echo "Available files:"
    ls "${assembly_dir}"/*.contigs.renamed.fa | nl
    exit 1
fi

if [ ! -f "$input_fasta" ]; then
    echo "Error: Input file does not exist: $input_fasta"
    exit 1
fi

echo "Processing file: $input_fasta"

sample=$(basename "${input_fasta}" .contigs.renamed.fa)
outdir="${output_base_dir}/${sample}"
hmm_file="${hmmdir}/COI.hmm"

echo "Sample identifier: $sample"
echo "Output directory: $outdir"

# Check if HMM file exists
if [ ! -f "$hmm_file" ]; then
    echo "Error: HMM file not found: $hmm_file"
    exit 1
fi

mkdir -p "${outdir}"

# Define output file paths
output_log=${outdir}/${sample}_${gene}.log
output_tbl=${outdir}/${sample}_${gene}.tbl
output_aliscores=${outdir}/${sample}_${gene}_aliscores.tsv
output_bed=${outdir}/${sample}_${gene}.bed
output_extracted=${outdir}/${sample}_${gene}_extracted_sequences.fna
input_fasta_copied="${outdir}/$(basename "${input_fasta}")"

echo "Starting COI extraction for sample: $sample"

# Copy and clean input fasta with seqtk
echo "Cleaning input fasta with seqtk..."
source activate seqtk
seqtk seq -l 0 "${input_fasta}" > "${input_fasta_copied}"
sed -i 's/\r$//' "${input_fasta_copied}"

# Verify copied file exists and has content
if [ ! -s "$input_fasta_copied" ]; then
    echo "Error: Failed to create or copied fasta is empty: $input_fasta_copied"
    exit 1
fi

echo "Running nhmmer..."
# Run nhmmer
conda deactivate && source activate hmmer
nhmmer \
 -o "${output_log}" \
 -T 10 \
 --tblout "${output_tbl}" \
 --aliscoresout "${output_aliscores}" \
 "${hmm_file}" \
 "${input_fasta_copied}"

# Check if nhmmer completed successfully
if [ $? -ne 0 ]; then
    echo "Error: nhmmer failed for sample $sample"
    rm -f "${input_fasta_copied}"
    exit 1
fi

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

# Check if any hits were found
hit_count=$(wc -l < "${output_bed}")
echo "Found $hit_count COI hits"

if [ "$hit_count" -gt 0 ]; then
    # Extract sequences using bedtools
    echo "Extracting sequences with bedtools..."
    conda deactivate && source activate bedtools
    bedtools getfasta \
     -fi "${input_fasta_copied}" \
     -bed "${output_bed}" \
     -s \
     -name \
     -fo "${output_extracted}"
    
    # Check if extraction was successful
    if [ $? -ne 0 ]; then
        echo "Error: bedtools getfasta failed for sample $sample"
        rm -f "${input_fasta_copied}"
        exit 1
    fi
    
    # Check if extracted file has content
    extracted_count=$(grep -c "^>" "${output_extracted}" 2>/dev/null || echo "0")
    echo "Extracted $extracted_count COI sequences"
else
    echo "No COI hits found - creating empty output file"
    touch "${output_extracted}"
fi

# Remove copied fasta to save space 
rm -f "${input_fasta_copied}"

echo "COI extraction completed successfully for sample: $sample"
