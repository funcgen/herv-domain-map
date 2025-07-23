#!/usr/bin/env python3
"""
Extract BED entries of internal (-int / _internal) ERV insertions
from a RepeatMasker .out file (LTR/ERV class).
"""

import re
import argparse
from tqdm import tqdm

def is_internal(name: str) -> bool:
    name = name.lower()
    return "-int" in name or "internal" in name

def extract_internal_bed(rm_out, output_bed):
    entries = []

    with open(rm_out) as fh:
        for line in fh:
            # Skip comment lines or headers (typically start with SW or score)
            if line.strip() == "" or line.startswith("SW") or line.startswith("score"):
                continue

            fields = re.split(r"\s+", line.strip())
            if len(fields) < 11:
                continue  # skip malformed lines

            repeat_name  = fields[9]
            repeat_class = fields[10]

            if not (is_internal(repeat_name) and repeat_class.startswith("LTR/ERV")):
                continue

            divergence = fields[1]
            chrom      = fields[4]
            start      = int(fields[5]) - 1  # BED is 0-based
            end        = int(fields[6])
            strand     = '-' if fields[8].upper() == 'C' else '+'

            name = (
                f"{repeat_name}_pos_{chrom}_{start+1}_{end}"
                f"_strand_{strand}_divergence_{divergence}"
            )
            score = "0"

            entries.append([chrom, str(start), str(end), name, score, strand])


    # Write BED
    with open(output_bed, "w") as out:
        for e in entries:
            out.write("\t".join(e) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract BED of ERV internal sequences from RepeatMasker output.")
    parser.add_argument("--repeatmasker", "-r", required=True, help="RepeatMasker .out file")
    parser.add_argument("--output", "-o", required=True, help="Output BED file")
    args = parser.parse_args()

    extract_internal_bed(args.repeatmasker, args.output)

if __name__ == "__main__":
    main()
