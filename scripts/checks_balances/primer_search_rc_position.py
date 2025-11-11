#!/usr/bin/env python3
"""
Primer Search Script for FASTA/FASTQ Files

Searches for primer sequences in FASTA or FASTQ files (compressed or uncompressed)
and outputs which primers are found in each read, their positions, and read lengths.
Searches both forward and reverse complement orientations.

Usage:
    # Single file
    python3 primer_search.py -p primers.fa -f sample.fastq.gz -o results.csv
    
    # Directory of files
    python3 primer_search.py -p primers.fa -d /path/to/files/ -o results.csv

Requirements:
    - seqkit (command-line tool)
"""

import subprocess
import csv
import argparse
from pathlib import Path
from collections import defaultdict
import re


def reverse_complement(seq):
    """
    Generate reverse complement of a DNA sequence.
    Handles standard IUPAC ambiguity codes.
    """
    complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'N': 'N',
        # Lowercase versions
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
        'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
        'd': 'h', 'h': 'd', 'n': 'n'
    }
    return ''.join(complement.get(base, base) for base in reversed(seq))


def parse_fasta(fasta_file):
    """
    Read primer sequences from a FASTA file.
    
    Returns a dictionary: {primer_name: sequence}
    Removes asterisks (*) from sequences and takes first part of header before |
    """
    primers = {}
    current_name = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            if line.startswith('>'):
                # Header line: extract primer name (text before first |)
                current_name = line[1:].split('|')[0]
            elif current_name:
                # Sequence line: store sequence without asterisks
                primers[current_name] = line.replace('*', '')
    
    return primers


def detect_file_type(file_path):
    """
    Detect if file is FASTA or FASTQ, and if it's gzipped.
    
    Returns: 'fasta', 'fasta.gz', 'fastq', 'fastq.gz', or None
    """
    name_lower = Path(file_path).name.lower()
    
    # Check gzipped files
    if name_lower.endswith('.gz'):
        if any(ext in name_lower for ext in ['.fasta.gz', '.fa.gz', '.fna.gz']):
            return 'fasta.gz'
        elif any(ext in name_lower for ext in ['.fastq.gz', '.fq.gz']):
            return 'fastq.gz'
    
    # Check uncompressed files
    if any(name_lower.endswith(ext) for ext in ['.fasta', '.fa', '.fna']):
        return 'fasta'
    elif any(name_lower.endswith(ext) for ext in ['.fastq', '.fq']):
        return 'fastq'
    
    return None


def get_sequences_from_file(file_path, file_type):
    """
    Extract all sequences from a FASTA/FASTQ file with their headers.
    Returns dict: {header: sequence}
    """
    sequences = {}
    
    try:
        # Use seqkit to extract sequences in tab-delimited format
        if file_type in ['fasta.gz', 'fastq.gz']:
            cmd = f"zcat {file_path} | seqkit fx2tab"
        else:
            cmd = f"seqkit fx2tab {file_path}"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"  Warning: seqkit fx2tab returned error code {result.returncode}")
            return sequences
        
        # Parse tab-delimited output: header \t sequence
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    # Extract just the first word of the header (before any spaces)
                    header = parts[0].split()[0]
                    sequence = parts[1]
                    sequences[header] = sequence
        
        return sequences
        
    except Exception as e:
        print(f"Warning: Error reading sequences from {file_path}: {e}")
        return sequences


def find_primer_in_sequence(sequence, primer_seq):
    """
    Find positions where primer matches in the sequence using simple string search.
    Returns list of (position, strand) tuples where position is 0-based.
    strand is 'fwd' or 'rc' (reverse complement)
    
    Note: This does exact matching. For IUPAC ambiguity codes, we rely on
    seqkit grep to identify which reads contain the primer, then do simple
    string matching here for position finding.
    """
    positions = []
    seq_upper = sequence.upper()
    primer_upper = primer_seq.upper()
    
    # Search forward strand
    start = 0
    while True:
        pos = seq_upper.find(primer_upper, start)
        if pos == -1:
            break
        positions.append((pos, 'fwd'))
        start = pos + 1
    
    # Search reverse complement
    primer_rc = reverse_complement(primer_seq).upper()
    start = 0
    while True:
        pos = seq_upper.find(primer_rc, start)
        if pos == -1:
            break
        positions.append((pos, 'rc'))
        start = pos + 1
    
    return positions


def find_reads_with_primer(file_path, primer_seq, file_type):
    """
    Search for primer sequence in a file using seqkit grep.
    
    Returns a dict: {read_header: sequence} for reads containing the primer.
    seqkit grep -s -d searches both strands and handles IUPAC codes.
    """
    matches = {}
    
    try:
        # Build seqkit command - use fx2tab to get both header and sequence
        if file_type in ['fasta.gz', 'fastq.gz']:
            cmd = f"zcat {file_path} | seqkit grep -s -d -p {primer_seq} | seqkit fx2tab"
        else:
            cmd = f"seqkit grep -s -d -p {primer_seq} {file_path} | seqkit fx2tab"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"  Warning: seqkit returned error code {result.returncode}")
            if result.stderr:
                print(f"  Error message: {result.stderr}")
            return matches
        
        # Parse tab-delimited output: header \t sequence
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    header = parts[0].split()[0]
                    sequence = parts[1]
                    matches[header] = sequence
        
        return matches
        
    except Exception as e:
        print(f"Warning: Error searching {file_path}: {e}")
        return matches


def process_single_file(primers_file, input_file, output_file):
    """Process one FASTA/FASTQ file and write results to CSV."""
    
    # Step 1: Load primers
    primers = parse_fasta(primers_file)
    print(f"Loaded {len(primers)} primers")
    
    # Step 2: Detect file type
    file_path = Path(input_file)
    filename = file_path.name
    file_type = detect_file_type(input_file)
    
    if file_type is None:
        print(f"Error: Unsupported file type for {filename}")
        print("Supported: .fasta, .fa, .fastq, .fq (and .gz versions)")
        return
    
    print(f"File type: {file_type}")
    print(f"Processing: {filename}\n")
    
    # Step 3: Search for each primer and record positions
    # Store results as: read_header -> [(primer_name, position, strand, read_length)]
    read_primers = defaultdict(list)
    
    for i, (primer_name, primer_seq) in enumerate(primers.items(), 1):
        primer_rc = reverse_complement(primer_seq)
        
        print(f"Searching primer {i}/{len(primers)}: {primer_name}")
        print(f"  Forward:  {primer_seq}")
        print(f"  Rev Comp: {primer_rc}")
        
        # Get reads that contain this primer (with sequences)
        matches = find_reads_with_primer(str(file_path), primer_seq, file_type)
        
        print(f"  Found in {len(matches)} reads")
        
        # For each matching read, find primer position(s)
        for header, sequence in matches.items():
            positions = find_primer_in_sequence(sequence, primer_seq)
            
            # If positions not found by simple search, just record position as 0
            # (this can happen with IUPAC codes where exact match won't work)
            if not positions:
                # Record as unknown position
                read_primers[header].append((primer_name, 0, 'unknown', len(sequence)))
            else:
                for pos, strand in positions:
                    read_primers[header].append((primer_name, pos, strand, len(sequence)))
    
    # Step 4: Write results to CSV
    write_results_to_csv(output_file, {filename: read_primers})
    
    print(f"\nDone! Results written to: {output_file}")
    print(f"Total reads with primers: {len(read_primers)}")


def process_directory(primers_file, input_dir, output_file):
    """Process all FASTA/FASTQ files in a directory and write combined results."""
    
    # Step 1: Load primers
    primers = parse_fasta(primers_file)
    print(f"Loaded {len(primers)} primers\n")
    
    # Step 2: Find all supported files in directory
    input_path = Path(input_dir)
    supported_extensions = [
        '.fasta', '.fa', '.fna', '.fastq', '.fq',
        '.fasta.gz', '.fa.gz', '.fna.gz', '.fastq.gz', '.fq.gz'
    ]
    
    files_to_process = [
        f for f in input_path.iterdir()
        if f.is_file() and any(f.name.lower().endswith(ext) for ext in supported_extensions)
    ]
    
    if not files_to_process:
        print(f"No FASTA/FASTQ files found in {input_dir}")
        return
    
    print(f"Found {len(files_to_process)} files to process\n")
    
    # Step 3: Process each file
    all_results = defaultdict(lambda: defaultdict(list))
    
    for file_idx, file_path in enumerate(files_to_process, 1):
        filename = file_path.name
        file_type = detect_file_type(file_path)
        
        print(f"[{file_idx}/{len(files_to_process)}] Processing: {filename} ({file_type})")
        
        # Search for each primer
        for primer_name, primer_seq in primers.items():
            matches = find_reads_with_primer(str(file_path), primer_seq, file_type)
            
            for header, sequence in matches.items():
                positions = find_primer_in_sequence(sequence, primer_seq)
                
                if not positions:
                    all_results[filename][header].append((primer_name, 0, 'unknown', len(sequence)))
                else:
                    for pos, strand in positions:
                        all_results[filename][header].append((primer_name, pos, strand, len(sequence)))
        
        print(f"  Found primers in {len(all_results[filename])} reads\n")
    
    # Step 4: Write all results to CSV
    write_results_to_csv(output_file, all_results)
    
    total_reads = sum(len(reads) for reads in all_results.values())
    print(f"Done! Results written to: {output_file}")
    print(f"Total reads with primers: {total_reads}")


def write_results_to_csv(output_file, results_dict):
    """
    Write search results to CSV file.
    
    results_dict format: {filename: {read_header: [(primer_name, position, strand, read_length)]}}
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header row
        header_row = ['file', 'read_header', 'primer', 'position', 'strand', 'read_length']
        writer.writerow(header_row)
        
        # Write data rows (one per primer match)
        for filename in sorted(results_dict.keys()):
            for read_header, primer_info_list in sorted(results_dict[filename].items()):
                for primer_name, position, strand, read_length in primer_info_list:
                    # Position is 0-based, convert to 1-based for output
                    # If position is 0 (unknown), keep it as 0 in output
                    output_pos = position + 1 if position > 0 else 0
                    row = [filename, read_header, primer_name, output_pos, strand, read_length]
                    writer.writerow(row)


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Find primers in FASTA/FASTQ files with position information',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single file
  python3 primer_search.py -p primers.fa -f sample.fastq.gz -o results.csv
  
  # Process all files in a directory
  python3 primer_search.py -p primers.fa -d /path/to/files/ -o results.csv

Note: This script uses seqkit grep -s -d which searches both forward and reverse
      complement strands and handles IUPAC ambiguity codes.
      Positions are reported as 1-based (first base = position 1).
      Position 0 means the exact position couldn't be determined (IUPAC ambiguity).
        """
    )
    
    parser.add_argument('-p', '--primers', required=True, 
                        help='FASTA file containing primer sequences')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', 
                       help='Single FASTA/FASTQ file to search')
    group.add_argument('-d', '--directory', 
                       help='Directory containing multiple files to search')
    
    parser.add_argument('-o', '--output', required=True, 
                        help='Output CSV file for results')
    
    args = parser.parse_args()
    
    if args.file:
        process_single_file(args.primers, args.file, args.output)
    else:
        process_directory(args.primers, args.directory, args.output)


if __name__ == "__main__":
    main()