#!/usr/bin/env python3
"""
Extracts strand-corrected sequences from a merged BED file
of ERV internal regions and writes them to a multi-FASTA.
"""

import argparse
from pyfaidx import Fasta
from Bio.Seq import Seq
from tqdm import tqdm

def rc(seq: str) -> str:
    """Reverse-complement a DNA string."""
    return str(Seq(seq).reverse_complement())

def extract_sequences_from_bed(bed_file, genome_fasta, output_fasta):
    genome = Fasta(genome_fasta)

    # First pass to count total entries for progress bar
    with open(bed_file) as f:
        total = sum(1 for _ in f)

    with open(bed_file) as f, open(output_fasta, "w") as out:
        pbar = tqdm(total=total, desc="Extracting sequences", unit="seq")

        for line in f:
            if line.strip() == "":
                continue

            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue

            chrom, start, end, name, _, strand = fields
            start = int(start)
            end = int(end)

            if chrom not in genome:
                print(f"Warning: {chrom} not found in FASTA – skipping.")
                pbar.update(1)
                continue

            seq = genome[chrom][start:end].seq
            if strand == '-':
                seq = rc(seq)

            out.write(f">{name}\n{seq}\n")
            pbar.update(1)

        pbar.close()

def main():
    parser = argparse.ArgumentParser(description="Extract strand-corrected sequences from BED file.")
    parser.add_argument("--bed", "-b", required=True, help="Input BED file (6-column, merged ERV internal regions)")
    parser.add_argument("--genome", "-g", required=True, help="Reference genome FASTA file (indexed with pyfaidx)")
    parser.add_argument("--output", "-o", required=True, help="Output FASTA file")
    args = parser.parse_args()

    extract_sequences_from_bed(args.bed, args.genome, args.output)

if __name__ == "__main__":
    main()
