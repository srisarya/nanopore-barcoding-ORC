# Field workshop on the taxonomy and natural history of freshwater and limno-terrestrial meiofauna 🪱 🏞️
## *This is the code used to process COI/28S/18S amplicon data from the 2025 NERC/NHM Lake District Freshwater Meiofauna workshop*

## Installation and setup
* You'll need this software:
  * pychopper v2.7.0 : `conda install bioconda::pychopper`
  * cutadapt v4.9.0: `conda install bioconda::cutadapt`
  * Amplicon Sorter (no releases, separate conda packages to install)
    `conda create -n amplicon_sorter && conda install bioconda::python-edlib biopython matplotlib` 
  * seqkit v2.9.0: `conda install bioconda::seqkit`
  * cd-hit v4.8.1: `conda install bioconda::cd-hit`
  * barrnap v0.9: `conda install bioconda::barrnap`
  * seqtk v1.4-r122: `conda install bioconda::seqtk`
  * hmmer v3.1b2: `conda install bioconda::hmmer`
  * bedtools v2.31.1: `conda install bioconda::bedtools`

## How to structure
- Run the scripts in the parent directory to the cloned github repo (so the filepaths remain relative)
- This workflow is written for a SLURM HPC system (I used the CropDiversity cluster)
- To run this code, git clone this whole repo, keeping the adapter directory and the scripts directory in the relative places they're in so that the codes run smoothly (or, edit the code yourself only in your cloned remote repo!

## The workflow
![bioinformatic workflow for processing DNA barcoding, from raw data to reorienting, demuxing, consensus building, cleaning, getting gene sequences, and selecting the best hits](<new strategy DNA Barcoding.png>)

  1. `01_pychopper.sh`: This script reorients and quality-score trims cDNA reads into the same orientation
     * This script takes in 2 extra files, which are hardcoded into the script:
       * one which has the orientation of adapters sequences listed (`M13_config_for_pychopper.txt`)
       * one which has the non-variable sequences for the adapters (`M13_seqs_for_pychopper.fa`)
  
  2. `02_cutadapt_loop.sh`: This script is a loop which first demultiplexes the reads based on the 5' SP5 primers, and then one-by-one, demultiplexes each output from the SP5 demultiplexing by the SP27 3' primers. It also trims adapters on-the-fly.
       * This script takes in 2 extra files, which are hardcoded into the script: the 5' and 3' adapter sequences, `M13_amplicon_indices_forward.fa` and the reverse complement of the 3' sequences (since we reoriented before with pychopper, `M13_amplicon_indices_reverse_rc.fa`)
       * You may notice that we do bin reads into adapter combinations which do not exist! Though this is inefficient, it's not the end of the world, since we have sequenced deeply enough per individual that enough reads are binned into the correct adapter combinations
  
  3. `03a_Amplicon_sorter_rRNAs.sh` and/or `03b_Amplicon_sorter_COIs.sh`
     * The workflow branches here, to either cluster all the rRNA reads together, before using a HMM to extract exact genes, OR to cluster all COI fragments together. In the wet-lab protocol:
      * For the rRNAs, we use 2 sets of primers to amplify overlapping segments of the nuclear rRNA cistron 
      * For the COIs, we use 2 alternate options for a forward primer, and one option for a reverse primer, to amplify 'redundant' sequences of a COI barcode segment (the 2 alternate options for the forward primer allow for matching more taxa than just one)
     * The amplicon sorter step splits samples' reads by size and clusters them by sequence
  
  4. `04a_cleaning_rRNAs.sh` and/or `04b_cleaning_COIs.sh`
       * This script takes each clustered (amplicon-sorted) file and removes the primer sequences from the amplicon sequences, since these are synthetic
       * for the 04b_cleaning_COIs.sh, we add an additional step to collapse sequences to make them non-redundant
  
  5. `05a_barrnap_rRNA_extract.sh` and/or `05b_nhmmer_COI_extract`: 
       * 05a_barrnap_rRNA_extract.sh: this script uses barrnap version 0.9. It takes the assembled contigs, and uses a HMM to extract sequences matching 28S and 18S rRNA profiles.
       * 05b_nhmmer_COI_extract.sh: this script runs an nhmmer search against custom pre-built COI HMM profiles and extracts the sequences
        * Speak to me if you'd like to know how I did this, since it's part of a larger project :) 

## Primer schematics
### rRNA amplification
![schematic for amplification of 18S/5.8S/partial 28S rRNA cistron, showing wiggly purple line representing DNA, yelllow boxes where genes are (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify genes, not to scale](18S-58S_28Spartial.png)
![schematic for amplification of full 28S rRNA gene, showing wiggly purple line representing DNA, yelllow boxes where genes are (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify genes, not to scale](28S_full.png)

### COI amplification
![schematic for amplification of partial COI gene, showing wiggly purple line representing DNA, yelllow boxes where COI gene is (overlaid on the purple DNA), and red arrows labelled with primer sequence name at the positions where the primers sit to amplify COI, not to scale](COI.png)