# Field workshop on the taxonomy and natural history of freshwater and limno-terrestrial meiofauna 🪱 🏞️
## *This is the code used to process COI/28S/18S amplicon data from the 2025 NERC/NHM Lake District Freshwater Meiofauna workshop*

## Installation and setup
* This workflow is written for a SLURM HPC system (I used the CropDiversity cluster)
* To run this code, git clone this whole repo, keeping the adapter directory and the scripts directory in the relative places they're in so that the codes run smoothly (or, edit the code yourself only in your cloned remote repo!
* You'll need this software:
  * pychopper v2.7.0 : `conda install bioconda::pychopper`
  * cutadapt v4.9.0: `conda install bioconda::cutadapt`
  * SPAdes v4.2.0: see [SPAdes docs](https://ablab.github.io/spades/installation.html) for guidance
  * seqkit v2.9.0: `conda install bioconda::seqkit`
  * barrnap v0.9: `conda install bioconda::barrnap`
  * seqtk v1.4-r122: `conda install bioconda::seqtk`
  * hmmer v3.1b2: `conda install bioconda::hmmer`
  * bedtools v2.31.1: `conda install bioconda::bedtools`

## Running the code
- Make sure to edit the code IN YOUR REMOTE GITHUB REPO ONLY and change the queueing system section at the top accordingly
- Run the scripts in the parent directory to the cloned github repo (so the filepaths remain relative)
- The order of job submissions is as follows:
  1. 01_pychopper.sh: This script reorients and quality-score trims cDNA reads into the same orientation
     * This script takes in 2 extra files, which are hardcoded into the script:
       * one which has the orientation of adapters sequences listed (`M13_config_for_pychopper.txt`)
       * one which has the non-variable sequences for the adapters (`M13_seqs_for_pychopper.fa`)
  2. 02_cutadapt_loop.sh: This script is a loop which first demultiplexes the reads based on the 5' SP5 primers, and then one-by-one, demultiplexes each output from the SP5 demultiplexing by the SP27 3' primers. It also trims adapters on-the-fly.
       * This script takes in 2 extra files, which are hardcoded into the script: the 5' and 3' adapter sequences, `M13_amplicon_indices_forward.fa` and the reverse complement of the 3' sequences (since we reoriented before with pychopper, `M13_amplicon_indices_reverse_rc.fa`)
       * You may notice that we do bin reads into adapter combinations which do not exist! Though this is inefficient, it's not the end of the world, since we have sequenced deeply enough per individual that enough reads are binned into the correct adapter combinations
       * Speak to me if you'd like to know more about this step :)
  4. 03_SPAdes.sh: This script takes each cleaned, demultiplexed, trimmed sample of reads and assembles the contigs into our 3 amplicons (and whatever else may be in your sample...)
  5. 04.1.1_barrnap_rRNA_extract.sh: this script uses barrnap version 0.9. It takes the assembled contigs, and uses a HMM to extract sequences matching 28S and 18S rRNA profiles.
  6. 04.1.2_barrnap_reorganisation.sh: this script reorganises barrnap outputs from having all the rRNAs from one sample in one file to per-rRNA-type files (optional)
  7. 04.2_nhmmer_COI_extract.sh: this script runs an nhmmer search against custom pre-built COI HMM profiles and extracts the sequences
       * Speak to me if you'd like to know how I did this, since it's part of a larger project :) 
