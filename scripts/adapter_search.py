#!/usr/bin/env python3

"""
Adapter Search Script with Reverse Complement Support

Usage: python adapter_search_dir.py <adapter_file> <input_directory> <output_file>
- adapter_file: FASTA file with adapter sequences
- input_directory: directory containing .gz files to search
- output_file: results will be written here

Requirements: biopython
Install with: pip install biopython
"""

import sys
import os
import gzip
import re
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def validate_requirements():
    """Check if required modules are available"""
    try:
        import Bio
    except ImportError:
        print("Error: Biopython is required but not installed")
        print("Install with: pip install biopython")
        sys.exit(1)

def parse_arguments():
    """Parse and validate command line arguments"""
    if len(sys.argv) != 4:
        print("Usage: python3 adapter_search_dir.py <adapter_file> <input_directory> <output_file>")
        print("Example: python3 adapter_search_dir.py adapters.fa /path/to/fastq_files/ results.txt")
        sys.exit(1)
    
    adapter_file = sys.argv[1]
    input_dir = sys.argv[2]
    output_file = sys.argv[3]
    
    # Validate inputs
    if not os.path.isfile(adapter_file):
        print(f"Error: Adapter file '{adapter_file}' not found")
        sys.exit(1)
    
    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' not found")
        sys.exit(1)
    
    return adapter_file, input_dir, output_file

def load_adapters(adapter_file):
    """Load adapter sequences from FASTA file"""
    adapters = {}
    
    try:
        with open(adapter_file, 'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq_str = str(record.seq).upper()
                
                # Validate sequence contains only valid DNA bases
                valid_bases = set('ATGCNRYSWKMBDHV')
                seq_bases = set(seq_str)
                
                if not seq_bases.issubset(valid_bases):
                    invalid = seq_bases - valid_bases
                    print(f"Warning: Invalid DNA bases found in {record.id}: {', '.join(invalid)}")
                    print(f"Skipping adapter: {record.id}")
                    continue
                
                adapters[record.id] = seq_str
                
    except Exception as e:
        print(f"Error reading adapter file: {e}")
        sys.exit(1)
    
    if not adapters:
        print("Error: No valid adapters found in the input file")
        sys.exit(1)
    
    return adapters

def get_reverse_complement(sequence):
    """Generate reverse complement using Biopython"""
    try:
        seq = Seq(sequence)
        return str(seq.reverse_complement())
    except Exception as e:
        print(f"Error generating reverse complement for {sequence}: {e}")
        return None

def search_in_file(filepath, sequence):
    """Search for a sequence in a gzipped file"""
    count = 0
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') as f:
            for line in f:
                # Count occurrences in the line (case-insensitive)
                count += len(re.findall(re.escape(sequence), line, re.IGNORECASE))
    except Exception as e:
        print(f"Warning: Error reading {filepath}: {e}")
        return 0
    
    return count

def get_gz_files(input_dir):
    """Get list of .gz files, excluding those with 'unknown' in name and empty files"""
    gz_files = []
    input_path = Path(input_dir)
    
    for filepath in input_path.glob("*.gz"):
        # Skip files with 'unknown' in name
        if 'unknown' in filepath.name.lower():
            continue
        
        # Skip empty files
        if filepath.stat().st_size == 0:
            continue
            
        gz_files.append(filepath)
    
    return sorted(gz_files)

def write_header(output_file, input_dir):
    """Write header to output file"""
    with open(output_file, 'w') as f:
        f.write(f"# Adapter search results - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# Format: ADAPTER_SEQUENCE | FILENAME | FORWARD_COUNT | REVERSE_COMPLEMENT_COUNT | TOTAL_COUNT\n")
        f.write(f"# Input directory: {input_dir}\n")
        f.write("# Generated using Python with Biopython for reverse complement calculation\n")
        f.write("\n")

def main():
    # Validate requirements and parse arguments
    validate_requirements()
    adapter_file, input_dir, output_file = parse_arguments()
    
    print(f"Adapter file: {adapter_file}")
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print()
    
    # Load adapters
    adapters = load_adapters(adapter_file)
    print(f"Loaded {len(adapters)} adapter sequences")
    
    # Get list of files to search
    gz_files = get_gz_files(input_dir)
    if not gz_files:
        print("Warning: No .gz files found in the input directory")
    else:
        print(f"Found {len(gz_files)} .gz files to search")
    print()
    
    # Initialize output file
    write_header(output_file, input_dir)
    
    # Process each adapter
    with open(output_file, 'a') as f:
        for adapter_name, adapter_seq in adapters.items():
            print(f"Searching for adapter: {adapter_name}")
            print(f"  Forward sequence: {adapter_seq}")
            
            # Generate reverse complement
            adapter_rc = get_reverse_complement(adapter_seq)
            if adapter_rc is None:
                print(f"  Skipping {adapter_name} due to reverse complement error")
                continue
                
            print(f"  Reverse complement: {adapter_rc}")
            
            # Write adapter info to output
            f.write(f"## Adapter: {adapter_name}\n")
            f.write(f"## Forward: {adapter_seq}\n")
            f.write(f"## Reverse complement: {adapter_rc}\n")
            
            if not gz_files:
                print("  No .gz files found in directory")
                f.write(f"{adapter_name} | NO_FILES_FOUND | 0 | 0 | 0\n")
                f.write("\n")
                continue
            
            # Search in each file
            total_files_with_matches = 0
            for filepath in gz_files:
                filename = filepath.name
                
                # Search for both forward and reverse complement
                forward_count = search_in_file(filepath, adapter_seq)
                reverse_count = search_in_file(filepath, adapter_rc)
                total_count = forward_count + reverse_count
                
                # Write results
                f.write(f"{adapter_name} | {filename} | {forward_count} | {reverse_count} | {total_count}\n")
                
                if total_count > 0:
                    print(f"  Found {total_count} total matches in {filename} (forward: {forward_count}, reverse: {reverse_count})")
                    total_files_with_matches += 1
                else:
                    print(f"  No matches in {filename}")
            
            print(f"  Summary: Found matches in {total_files_with_matches}/{len(gz_files)} files")
            f.write("\n")
            print()
    
    print(f"Search complete. Results written to: {output_file}")

if __name__ == "__main__":
    main()
