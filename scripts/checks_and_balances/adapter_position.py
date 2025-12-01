#!/usr/bin/env python3

"""
Script to find adapter positions in FASTQ files and categorize them as start/middle/end.
This version checks BOTH the forward sequence and reverse complement of each adapter.

Usage: ./script.py adapters.fasta fastq_dir/ output.csv

Requirements: pip install biopython
"""

import sys
import gzip
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq  # For reverse complement functionality

def reverse_complement(seq):
    """
    Generate the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string (e.g., 'ACGT')
    
    Returns:
        Reverse complement string (e.g., 'ACGT' -> 'ACGT')
    
    Example:
        'ATCG' -> complement is 'TAGC' -> reverse is 'CGAT'
    """
    # Bio.Seq handles the reverse complement for us
    return str(Seq(seq).reverse_complement())

def parse_fasta(fasta_file):
    """
    Parse a multifasta file and extract adapter names and sequences.
    Also generates reverse complements for each adapter.
    
    Args:
        fasta_file: Path to the fasta file containing adapter sequences
    
    Returns:
        Dictionary mapping adapter names to their sequences (both strands)
        e.g., {'SP5_001_fwd': 'ACGTACGT', 'SP5_001_rev': 'ACGTACGT'}
    """
    adapters = {}
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            # Store forward sequence
            adapters[f"{record.id}_fwd"] = str(record.seq)
            # Store reverse complement
            adapters[f"{record.id}_rev"] = reverse_complement(str(record.seq))
    return adapters

def process_fastq(fastq_file, adapter_name, adapter_seq, out_file, threshold=50):
    """
    Process one FASTQ file to find all occurrences of one adapter.
    
    Args:
        fastq_file: Path to FASTQ file (can be .gz compressed)
        adapter_name: Name of the adapter being searched (includes _fwd or _rev suffix)
        adapter_seq: Sequence of the adapter to search for
        out_file: Open file handle to write results to
        threshold: Number of bases from start/end to categorize as "start"/"end" (default: 50)
    
    Returns:
        Number of hits found in this file
    """
    filename = os.path.basename(fastq_file)
    hit_count = 0
    
    # Determine whether to use gzip or regular file opening based on extension
    opener = gzip.open if fastq_file.endswith('.gz') else open
    
    with opener(fastq_file, 'rt') as f:
        # SeqIO.parse reads the FASTQ file and returns an iterator of SeqRecord objects
        for record in SeqIO.parse(f, "fastq"):
            read_id = record.id
            read_seq = str(record.seq).upper()
            read_length = len(read_seq)
            
            # Find all occurrences of the adapter in this read
            adapter_upper = adapter_seq.upper()
            pos = 0
            
            while pos < len(read_seq):
                # Find the next occurrence of the adapter
                pos = read_seq.find(adapter_upper, pos)
                
                if pos == -1:
                    break
                
                # Convert to 1-indexed position
                position = pos + 1
                
                # Categorize the position
                if position <= threshold:
                    category = "start"
                elif position >= read_length - threshold:
                    category = "end"
                else:
                    category = "middle"
                
                # Write this match to the output CSV
                # adapter_name will include _fwd or _rev suffix
                out_file.write(f"{adapter_name},{filename},{position},{read_length},{category}\n")
                hit_count += 1
                
                pos += 1
    
    return hit_count

def main():
    """
    Main function to orchestrate the adapter position analysis.
    Searches for both forward and reverse complement of each adapter.
    """
    if len(sys.argv) != 4:
        print("Usage: ./script.py adapters.fasta fastq_dir/ output.csv")
        sys.exit(1)
    
    adapter_fasta = sys.argv[1]
    fastq_dir = sys.argv[2]
    output_csv = sys.argv[3]
    
    # Step 1: Parse adapters (this now includes both forward and reverse complement)
    print("Parsing adapters and generating reverse complements...")
    adapters = parse_fasta(adapter_fasta)
    # adapters dict now has twice as many entries (fwd and rev for each original)
    original_count = len(adapters) // 2
    print(f"Found {original_count} adapters ({len(adapters)} including reverse complements)\n")
    
    # Step 2: Get list of all FASTQ files
    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq.gz"))
    fastq_files.extend(glob.glob(os.path.join(fastq_dir, "*.fastq")))
    
    print(f"Found {len(fastq_files)} FASTQ files\n")
    
    # Step 3: Process each adapter (both fwd and rev) against each FASTQ file
    with open(output_csv, 'w') as out:
        # Write CSV header
        # Note: adapter column will now include _fwd or _rev suffix
        out.write("adapter,file,position,read_length,category\n")
        
        for adapter_name, adapter_seq in adapters.items():
            print(f"Processing adapter: {adapter_name}")
            
            for fastq_file in fastq_files:
                filename = os.path.basename(fastq_file)
                print(f"  {filename}...", end='', flush=True)
                
                hits = process_fastq(fastq_file, adapter_name, adapter_seq, out, threshold=50)
                
                print(f" {hits} hits")
    
    print(f"\nComplete! Results in {output_csv}")
    print("Note: Adapter names have _fwd or _rev suffix to indicate strand")

if __name__ == "__main__":
    main()