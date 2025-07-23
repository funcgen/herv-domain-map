#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA file to retain only sequences with 'ENV' in the header.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with only ENV sequences")

    args = parser.parse_args()

    with open(args.output, "w") as out_f:
        for record in SeqIO.parse(args.input, "fasta"):
            if "ENV" in record.id or "ENV" in record.description:
                SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    main()

