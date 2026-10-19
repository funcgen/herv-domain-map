#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA file to retain records whose domain gene field is env.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with only ENV sequences")

    args = parser.parse_args()

    with open(args.output, "w") as out_f:
        for record in SeqIO.parse(args.input, "fasta"):
            if len(record.id.rsplit('|',3)) == 4 and record.id.rsplit('|',3)[1].lower() == 'env':
                SeqIO.write(record, out_f, "fasta")

if __name__ == "__main__":
    main()

