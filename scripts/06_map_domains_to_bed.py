#!/usr/bin/env python3
"""
Convert filtered HMMscan TSV hits on ORFs back to genome BED, with gene type, conservedness, and coverage.
"""

import argparse
import re
import csv

# ---------- helper functions ----------
def parse_header(name):
    m = re.match(
        r'.*?_pos_(?P<chr>.+?)_(?P<int_start>\d+)_(?P<int_end>\d+)_strand_(?P<strand>[+-]).*?\|(?P<orf_start>\d+)-(?P<orf_end>\d+)$',
        name
    )
    if not m:
        raise ValueError(f"Header not understood: {name}")
    d = m.groupdict()
    return {
        "chr": d["chr"],
        "int_start": int(d["int_start"]),
        "int_end": int(d["int_end"]),
        "strand": d["strand"],
        "orf_start_in_int": int(d["orf_start"]),
        "orf_end_in_int": int(d["orf_end"]),
    }

def extract_ali_coords(line):
    fields = re.split(r'\s+', line.strip(), maxsplit=22)
    try:
        ali_from_aa = int(fields[17])
        ali_to_aa = int(fields[18])
        # fix reversed coordinates
        ali_start = min(ali_from_aa, ali_to_aa)
        ali_end   = max(ali_from_aa, ali_to_aa)
        return ali_start, ali_end
    except (IndexError, ValueError):
        raise ValueError(f"Could not extract alignment coordinates from line: {line}")

# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(description="Map filtered HMMscan hits on ORFs back to genome BED format, including gene type.")
    parser.add_argument("-i", "--input", required=True, help="Input filtered TSV file (from filter_hmmscan_hits_per_protein.py)")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()
       
    with open(args.input) as fh, open(args.output, "w", newline="") as out:
        bed = csv.writer(out, delimiter="\t")

        header = fh.readline()
        if header.lower().startswith("query"):
            pass  # skip header line

        for line in fh:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            query_name = fields[0]
            gene_class = fields[1]
            target_hmm = fields[2]
            score = int(float(fields[4]))
            coverage = fields[7]
            conservedness = fields[8]
            line_raw = fields[-1]

            # Extract alignment coordinates
            ali_from_aa, ali_to_aa = extract_ali_coords(line_raw)
            # Parse header info
            info = parse_header(query_name)

            # aa → nt within ORF (0-based)
            nt_from_orf = (ali_from_aa - 1) * 3
            nt_to_orf = ali_to_aa * 3



            # ORF offset to INT (0-based)
            nt_from_int = nt_from_orf + (info["orf_start_in_int"] - 1)
            nt_to_int   = nt_to_orf   + (info["orf_start_in_int"] - 1)


            # INT → genome
            if info["strand"] == "+":
                genome_from = info["int_start"] + nt_from_int
                genome_to   = info["int_start"] + nt_to_int
            else:
                genome_from = info["int_end"] - nt_to_int
                genome_to   = info["int_end"] - nt_from_int

            # BED output: 0-based start, 1-based end
            bed_start = min(genome_from, genome_to)
            bed_end   = max(genome_from, genome_to)

            # BED name: Query|HMM|GeneType|Conservedness|Coverage
            bed_name = f"{query_name}|{target_hmm}|{gene_class}|{conservedness}|{coverage}"
            bed.writerow([info["chr"], bed_start, bed_end, bed_name, score, info["strand"]])

    print(f"BED file with gene type, conservedness, and coverage written: {args.output}")

if __name__ == "__main__":
    main()

