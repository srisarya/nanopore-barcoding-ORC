#!/usr/bin/env python3
"""
Primer Search Script for FASTA/FASTQ Files

Searches for primer sequences in FASTA or FASTQ files (compressed or uncompressed)
and outputs which primers are found in each read. Searches both forward and reverse
complement orientations automatically via seqkit.

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


def find_reads_with_primer(file_path, primer_seq, file_type):
    """
    Search for primer sequence in a file using seqkit.
    
    Returns a list of read headers that contain the primer sequence.
    
    The seqkit flags do all the work:
    -s: search in sequence (not header)
    -d: treat pattern as degenerate sequence (handles IUPAC codes like N, R, Y, etc.)
    -p: the pattern to search for
    
    seqkit automatically searches BOTH forward and reverse complement strands,
    so we don't need to manually compute or search reverse complements.
    """
    try:
        # Build seqkit command based on file type
        if file_type in ['fasta.gz', 'fastq.gz']:
            # For gzipped files: decompress first, then search
            cmd = f"zcat {file_path} | seqkit grep -s -d -p {primer_seq}"
        else:
            # For regular files: search directly
            cmd = f"seqkit grep -s -d -p {primer_seq} {file_path}"
        
        # Run the command and capture output
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # Check for errors (but don't stop - some primers might just not be found)
        if result.returncode != 0:
            print(f"  Warning: seqkit returned error code {result.returncode}")
            if result.stderr:
                print(f"  Error message: {result.stderr}")
        
        # Extract read headers from seqkit output
        # seqkit grep returns the matching reads in FASTA/FASTQ format
        headers = []
        for line in result.stdout.split('\n'):
            # FASTA headers start with >, FASTQ headers start with @
            if line.startswith('>') or line.startswith('@'):
                # Remove the > or @ and take first word (header ID)
                header = line[1:].split()[0]
                headers.append(header)
        
        return headers
        
    except Exception as e:
        print(f"Warning: Error searching {file_path}: {e}")
        return []


def process_single_file(primers_file, input_file, output_file):
    """
    Process one FASTA/FASTQ file and write results to CSV.
    
    For each primer, we use seqkit to find all reads containing that primer,
    then store the results in a dictionary mapping read headers to the list
    of primers found in that read.
    """
    
    # Step 1: Load primers from the FASTA file
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
    
    # Step 3: Search for each primer
    # Store results as: read_header -> [list of primers found in that read]
    read_primers = defaultdict(list)
    
    for i, (primer_name, primer_seq) in enumerate(primers.items(), 1):
        print(f"Searching primer {i}/{len(primers)}: {primer_name}")
        print(f"  Sequence: {primer_seq}")
        
        # Find all reads containing this primer
        # seqkit handles IUPAC codes and searches both strands automatically
        headers = find_reads_with_primer(str(file_path), primer_seq, file_type)
        
        # Add this primer to each read's list
        for header in headers:
            read_primers[header].append(primer_name)
        
        print(f"  Found in {len(headers)} reads")
    
    # Step 4: Write results to CSV
    write_results_to_csv(output_file, {filename: read_primers})
    
    print(f"\nDone! Results written to: {output_file}")
    print(f"Total reads with primers: {len(read_primers)}")


def process_directory(primers_file, input_dir, output_file):
    """
    Process all FASTA/FASTQ files in a directory and write combined results.
    
    This is useful when you have multiple samples and want to search them all
    in one go. The output CSV will include a 'file' column so you can tell
    which reads came from which file.
    """
    
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
    # Store results as: filename -> {read_header -> [primers]}
    all_results = defaultdict(lambda: defaultdict(list))
    
    for file_idx, file_path in enumerate(files_to_process, 1):
        filename = file_path.name
        file_type = detect_file_type(file_path)
        
        print(f"[{file_idx}/{len(files_to_process)}] Processing: {filename} ({file_type})")
        
        # Search for each primer in this file
        for primer_name, primer_seq in primers.items():
            headers = find_reads_with_primer(str(file_path), primer_seq, file_type)
            
            for header in headers:
                all_results[filename][header].append(primer_name)
        
        print(f"  Found primers in {len(all_results[filename])} reads\n")
    
    # Step 4: Write all results to CSV
    write_results_to_csv(output_file, all_results)
    
    total_reads = sum(len(reads) for reads in all_results.values())
    print(f"Done! Results written to: {output_file}")
    print(f"Total reads with primers: {total_reads}")


def write_results_to_csv(output_file, results_dict):
    """
    Write search results to CSV file.
    
    The output format is designed to be flexible for different numbers of primers
    per read. Each read gets one row, with columns for each primer found.
    
    results_dict format: {filename: {read_header: [list of primers]}}
    """
    # Figure out the maximum number of primers found in any single read
    # This tells us how many primer columns we need
    max_primers = 0
    for file_results in results_dict.values():
        for primer_list in file_results.values():
            max_primers = max(max_primers, len(primer_list))
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header row
        header_row = ['file', 'read_header'] + [f'primer_{i+1}' for i in range(max_primers)]
        writer.writerow(header_row)
        
        # Write data rows (one per read that has at least one primer)
        for filename in sorted(results_dict.keys()):
            for read_header, primer_list in sorted(results_dict[filename].items()):
                # Build row: filename, read_header, primer1, primer2, ...
                row = [filename, read_header] + primer_list
                # Pad with empty strings if this read has fewer primers than max
                row += [''] * (len(header_row) - len(row))
                writer.writerow(row)


def main():
    """
    Main entry point for the script.
    
    Sets up command-line arguments and calls the appropriate processing function
    based on whether the user provided a single file or a directory.
    """
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description='Find primers in FASTA/FASTQ files (searches both strands, handles IUPAC codes)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single file
  python3 primer_search.py -p primers.fa -f sample.fastq.gz -o results.csv
  
  # Process all files in a directory
  python3 primer_search.py -p primers.fa -d /path/to/files/ -o results.csv

Note: This script uses seqkit grep which automatically searches both forward 
      and reverse complement strands, and properly handles IUPAC ambiguity codes
      (N, R, Y, W, S, K, M, B, D, H, V) in primer sequences.
        """
    )
    
    parser.add_argument('-p', '--primers', required=True, 
                        help='FASTA file containing primer sequences')
    
    # User must specify either -f OR -d, but not both
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--file', 
                       help='Single FASTA/FASTQ file to search')
    group.add_argument('-d', '--directory', 
                       help='Directory containing multiple files to search')
    
    parser.add_argument('-o', '--output', required=True, 
                        help='Output CSV file for results')
    
    args = parser.parse_args()
    
    # Run the appropriate function based on input type
    if args.file:
        process_single_file(args.primers, args.file, args.output)
    else:
        process_directory(args.primers, args.directory, args.output)


if __name__ == "__main__":
    main()
