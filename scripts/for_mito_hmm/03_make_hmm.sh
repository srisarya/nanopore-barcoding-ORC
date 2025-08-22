#!/bin/bash
#SBATCH --job-name=makehmm
#SBATCH --mem=4G
#SBATCH --cpus-per-task=12
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log

# setup variables: directories
gene="cox1"
basedir=Lakes_meiofauna_workshop/adapters_hmms/COI_seqs

# setup variables: files
input_raw=${basedir}/combined_${gene}.fna
input_filtered=${basedir}/"${gene}_lenfiltered.fna"
out_alignment=${basedir}/"${gene}_aligned.fna"
out_hmm=${basedir}/"${gene}.hmm"
hmm_log=${basedir}/"${gene}_hmmbuild.log"

# setup variables: temporary files
seqlengths=${basedir}/${gene}_seq_lengths.tsv
keep_headers=${basedir}/${gene}_keep_headers.txt

# 1. concatenate gene-fetch_results
cat ${basedir}/*/nucleotide/*fasta > ${input_raw}

# 2. filtering the gene-fetch file
    # Activate seqkit environment
    source activate seqkit
    # 1.1. Create a TSV with headers and lengths
        seqkit fx2tab -nl ${input_raw} > ${seqlengths}
    # 1.2. Extract headers of sequences from genomic sequences (so, full CDSs, not barcodes. also, no length filter for now since each gene is a different length)
        awk -F'\t' '$1 ~ /genome/ {print $1}' ${seqlengths} > ${keep_headers}
    # 1.3. Filter the original fasta to keep only those headers
        seqkit grep --by-name --pattern-file ${keep_headers} ${input_raw} > ${input_filtered}

# 3. aligning the sequences with MAFFT
conda deactivate && source activate mafft
mafft --auto \
 --thread 12 \
 ${input_filtered} > \
 ${out_alignment}

# 4. building the HMM with HMMER
conda deactivate && source activate hmmer
hmmbuild \
 --cpu 12 \
 -o ${hmm_log} \
 --dna \
 ${out_hmm} \
 ${out_alignment}
