#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Extract domain sequences from ORF protein sequences based on HMMER hit coordinates.")
    parser.add_argument("-t", "--tsv", required=True, help="Filtered HMMER hits TSV file (e.g., domains_gydb_filtered.tsv)")
    parser.add_argument("-f", "--fasta", required=True, help="ORF protein FASTA file (e.g., orfs_fixed_headers.fa)")
    parser.add_argument("-o", "--output", required=True, help="Output domain sequences FASTA file")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load ORF protein sequences into a dictionary
    orf_dict = {record.id: record.seq for record in SeqIO.parse(args.fasta, "fasta")}
    print(f"Loaded {len(orf_dict)} ORF protein sequences from {args.fasta}")

    count_written = 0
    skipped = 0

    with open(args.tsv) as tsv, open(args.output, "w") as out_f:
        for line in tsv:
            if line.startswith("#") or line.startswith("Query"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            query = parts[0]
            domain_name = parts[2]
            gene_name = parts[1]
            try:
                line_fields = parts[12].split()
                ali_from = int(line_fields[17])
                ali_to = int(line_fields[18])
            except (IndexError, ValueError):
                print(f"Skipping {query}: Could not parse coordinates from Line field")
                skipped += 1
                continue

            if query not in orf_dict:
                print(f"Skipping {query}: not found in ORF FASTA")
                skipped += 1
                continue

            orf_seq = orf_dict[query]
            start = min(ali_from, ali_to)
            end = max(ali_from, ali_to)
            domain_seq = orf_seq[start - 1 : end]
            header = f">{query}|{gene_name}|{domain_name}|{ali_from}-{ali_to}"
            out_f.write(f"{header}\n{domain_seq}\n")
            count_written += 1

    print(f"Domain sequences written to {args.output}")
    print(f"{count_written} domain sequences extracted")
    if skipped > 0:
        print(f"{skipped} entries skipped due to missing data")

if __name__ == "__main__":
    main()
