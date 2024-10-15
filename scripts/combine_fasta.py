# scripts/combine_fasta.py
'''
* Author: Evan Komp
* Created: 10/15/2024
* Company: National Renewable Energy Lab, Bioeneergy Science and Technology
* License: MIT
'''
"""Combine multiple FASTA files into a single file."""

import sys
from pathlib import Path
from Bio import SeqIO

def combine_fasta_files(input_files, output_file):
    """
    Combine multiple FASTA files into a single file.

    Args:
        input_files (list): List of input FASTA file paths.
        output_file (str): Path to the output combined FASTA file.
    """
    with open(output_file, 'w') as outfile:
        for infile in input_files:
            for record in SeqIO.parse(infile, 'fasta'):
                SeqIO.write(record, outfile, 'fasta')

def main():
    """Parse arguments and combine FASTA files."""
    if len(sys.argv) < 4:
        print("Usage: python combine_fasta.py <output_file> <input_file1> <input_file2> ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = [Path(f) for f in sys.argv[2:] if Path(f).exists()]

    combine_fasta_files(input_files, output_file)
    print(f"Combined {len(input_files)} files into {output_file}")

if __name__ == "__main__":
    main()