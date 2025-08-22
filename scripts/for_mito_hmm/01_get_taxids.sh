#!/bin/bash
#SBATCH --job-name=get_taxids
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log

# setup variables: directories
basedir=Lakes_meiofauna_workshop/adapters_hmms
taxid_dir=${basedir}/COI_seqs

mkdir  -p ${taxid_dir}

# setup variables: options
taxon="metazoa"
level="Order"
taxid=33208 # metazoan order taxid, needs to be found before the run

# setup variables: files
norank_file=${taxid_dir}/${taxon}_taxids_norank.txt
orders_file=${taxid_dir}/${taxon}_orders.txt

echo "taxid_dir: $taxid_dir"
echo "taxid_file: $taxid_file"

source activate taxonkit
# all metazoan order taxids
taxonkit list --ids "${taxid}" > ${norank_file}
cat ${norank_file} | taxonkit filter -E ${level} > ${orders_file}
