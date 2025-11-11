#!/usr/bin/env python3
"""
Amplicon Search Script for FASTA/FASTQ Files

Searches for amplicon regions using primer pairs in FASTA or FASTQ files 
(compressed or uncompressed) and outputs which amplicons are found in each read.

Usage:
    # Single file
    python3 primer_amplicon_search.py -p primers.fa -f sample.fastq.gz -o results.csv
    
    # Directory of files
    python3 primer_amplicon_search.py -p primers.fa -d /path/to/files/ -o results.csv

Requirements:
    - seqkit (command-line tool with amplicon subcommand)
"""

import subprocess
import csv
import argparse
import tempfile
from pathlib import Path
from collections import defaultdict


def parse_primer_fasta(fasta_file):
    """
    Read primer sequences from a FASTA file and organize into primer pairs.
    
    Returns a list of tuples: [(pair_name, fwd_name, fwd_seq, rvs_name, rvs_seq), ...]
    
    Expected format:
    - Headers with format: >primer_name|description
    - Primer pairs are defined by shared descriptions or manual grouping
    """
    primers = {}
    current_name = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line:
                continue
            
            if line.startswith('>'):
                # Extract full header (name|description)
                header = line[1:]
                parts = header.split('|')
                current_name = parts[0]
                # Store both name and full description
                primers[current_name] = {'header': header, 'seq': ''}
            elif current_name:
                primers[current_name]['seq'] = line.replace('*', '')
    
    # Define primer pairs based on the provided primers
    # Format: (pair_name, [fwd_primer_names], rvs_primer_name)
    primer_pairs = [
        # Moorea uses jgLCO1490 forward
        ("Moorea", ["jgLCO1490"], "jgHCO2198"),
        # Sauron uses COI-Sauron-S878 forward
        ("Sauron", ["COI-Sauron-S878"], "jgHCO2198"),
        # 18S_5.8S_28S_part
        ("18S_5.8S_28S_part", ["SSU_F04"], "28S_3RC"),
        # 28S
        ("28S", ["F63.2"], "R3264.2"),
    ]
    
    # Create list of (pair_name, fwd_name, fwd_seq, rvs_name, rvs_seq)
    pairs = []
    for pair_name, fwd_names, rvs_name in primer_pairs:
        for fwd_name in fwd_names:
            if fwd_name in primers and rvs_name in primers:
                pairs.append((
                    pair_name,
                    fwd_name,
                    primers[fwd_name]['seq'],
                    rvs_name,
                    primers[rvs_name]['seq']
                ))
            else:
                print(f"Warning: Missing primers for pair {pair_name} ({fwd_name}/{rvs_name})")
    
    return pairs


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


def get_all_read_headers(file_path, file_type):
    """
    Get all read headers from a FASTA/FASTQ file.
    
    Returns a set of read headers (without > or @).
    """
    try:
        # Use seqkit seq to extract just the IDs
        if file_type in ['fasta.gz', 'fastq.gz']:
            cmd = f"zcat {file_path} | seqkit seq -n -i"
        else:
            cmd = f"seqkit seq -n -i {file_path}"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Warning: Error getting read headers: {result.stderr}")
            return set()
        
        # Each line is a read ID
        headers = set(line.strip() for line in result.stdout.split('\n') if line.strip())
        return headers
        
    except Exception as e:
        print(f"Warning: Error getting read headers: {e}")
        return set()


def find_amplicons(file_path, fwd_seq, rvs_seq, file_type):
    """
    Search for amplicon regions using seqkit amplicon.
    
    Returns a set of read headers that contain the amplicon.
    """
    try:
        # Build seqkit amplicon command
        # --only-positive-strand: only search forward strand
        # --immediate-output: output results immediately
        # --bed: output in BED format (includes read ID)
        if file_type in ['fasta.gz', 'fastq.gz']:
            cmd = f"zcat {file_path} | seqkit amplicon --only-positive-strand --immediate-output --bed -F '{fwd_seq}' -R '{rvs_seq}'"
        else:
            cmd = f"seqkit amplicon --only-positive-strand --immediate-output --bed -F '{fwd_seq}' -R '{rvs_seq}' {file_path}"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # seqkit amplicon returns exit code 0 even with no matches
        # BED format: chrom (read ID), chromStart, chromEnd, name, score, strand
        headers = set()
        for line in result.stdout.split('\n'):
            line = line.strip()
            if line and not line.startswith('#'):
                # First column is the read ID
                parts = line.split('\t')
                if parts:
                    headers.add(parts[0])
        
        return headers
        
    except Exception as e:
        print(f"Warning: Error searching amplicons: {e}")
        return set()


def process_single_file(primers_file, input_file, output_file):
    """Process one FASTA/FASTQ file and write results to CSV."""
    
    # Step 1: Load primer pairs
    primer_pairs = parse_primer_fasta(primers_file)
    print(f"Loaded {len(primer_pairs)} primer pairs")
    
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
    
    # Step 3: Get all read headers to track which have no amplicons
    print("Getting all read headers...")
    all_headers = get_all_read_headers(str(file_path), file_type)
    print(f"Total reads in file: {len(all_headers)}\n")
    
    # Step 4: Search for each primer pair
    # Store results as: read_header -> [list of amplicon pairs found]
    read_amplicons = defaultdict(list)
    reads_with_amplicons = set()
    
    for i, (pair_name, fwd_name, fwd_seq, rvs_name, rvs_seq) in enumerate(primer_pairs, 1):
        print(f"Searching amplicon {i}/{len(primer_pairs)}: {pair_name}")
        print(f"  Forward ({fwd_name}): {fwd_seq}")
        print(f"  Reverse ({rvs_name}): {rvs_seq}")
        
        # Find all reads containing this amplicon
        headers = find_amplicons(str(file_path), fwd_seq, rvs_seq, file_type)
        
        # Add this amplicon to each read's list
        for header in headers:
            read_amplicons[header].append(pair_name)
            reads_with_amplicons.add(header)
        
        print(f"  Found in {len(headers)} reads\n")
    
    # Step 5: Identify reads with no amplicons
    reads_no_amplicons = all_headers - reads_with_amplicons
    print(f"Reads with amplicons: {len(reads_with_amplicons)}")
    print(f"Reads with NO amplicons: {len(reads_no_amplicons)}\n")
    
    # Step 6: Write results to CSV
    write_results_to_csv(output_file, {filename: read_amplicons}, {filename: reads_no_amplicons})
    
    print(f"Done! Results written to: {output_file}")


def process_directory(primers_file, input_dir, output_file):
    """Process all FASTA/FASTQ files in a directory and write combined results."""
    
    # Step 1: Load primer pairs
    primer_pairs = parse_primer_fasta(primers_file)
    print(f"Loaded {len(primer_pairs)} primer pairs\n")
    
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
    all_no_amplicons = defaultdict(set)
    
    for file_idx, file_path in enumerate(files_to_process, 1):
        filename = file_path.name
        file_type = detect_file_type(file_path)
        
        print(f"[{file_idx}/{len(files_to_process)}] Processing: {filename} ({file_type})")
        
        # Get all read headers
        all_headers = get_all_read_headers(str(file_path), file_type)
        reads_with_amplicons = set()
        
        # Search for each primer pair
        for pair_name, fwd_name, fwd_seq, rvs_name, rvs_seq in primer_pairs:
            headers = find_amplicons(str(file_path), fwd_seq, rvs_seq, file_type)
            for header in headers:
                all_results[filename][header].append(pair_name)
                reads_with_amplicons.add(header)
        
        # Track reads with no amplicons
        all_no_amplicons[filename] = all_headers - reads_with_amplicons
        
        print(f"  Total reads: {len(all_headers)}")
        print(f"  With amplicons: {len(all_results[filename])}")
        print(f"  No amplicons: {len(all_no_amplicons[filename])}\n")
    
    # Step 4: Write all results to CSV
    write_results_to_csv(output_file, all_results, all_no_amplicons)
    
    total_reads = sum(len(reads) for reads in all_results.values())
    total_no_amplicons = sum(len(reads) for reads in all_no_amplicons.values())
    print(f"Done! Results written to: {output_file}")
    print(f"Total reads with amplicons: {total_reads}")
    print(f"Total reads with NO amplicons: {total_no_amplicons}")


def write_results_to_csv(output_file, results_dict, no_amplicons_dict):
    """
    Write search results to CSV file.
    
    results_dict format: {filename: {read_header: [list of amplicons]}}
    no_amplicons_dict format: {filename: set(read_headers)}
    """
    # Calculate how many amplicon columns we need
    max_amplicons = 0
    for file_results in results_dict.values():
        for amplicon_list in file_results.values():
            max_amplicons = max(max_amplicons, len(amplicon_list))
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header row
        header_row = ['file', 'read_header', 'status'] + [f'amplicon_{i+1}' for i in range(max_amplicons)]
        writer.writerow(header_row)
        
        # Write reads with amplicons
        for filename in sorted(results_dict.keys()):
            for read_header, amplicon_list in sorted(results_dict[filename].items()):
                row = [filename, read_header, 'amplicon_found'] + amplicon_list
                row += [''] * (len(header_row) - len(row))
                writer.writerow(row)
        
        # Write reads with no amplicons
        for filename in sorted(no_amplicons_dict.keys()):
            for read_header in sorted(no_amplicons_dict[filename]):
                row = [filename, read_header, 'no_amplicon'] + [''] * max_amplicons
                writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description='Find amplicons in FASTA/FASTQ files using primer pairs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single file
  python3 primer_amplicon_search.py -p primers.fa -f sample.fastq.gz -o results.csv
  
  # Process all files in a directory
  python3 primer_amplicon_search.py -p primers.fa -d /path/to/files/ -o results.csv
  
Primer File Format:
  The script expects primers in FASTA format and will automatically pair them
  based on the primer names in the script. For custom pairing, edit the
  parse_primer_fasta() function.
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