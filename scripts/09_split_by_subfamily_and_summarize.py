#!/usr/bin/env python3

import argparse
import os
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Split BED by ERV subfamily and generate summary TSVs.")
    parser.add_argument("--input", "-i", required=True, help="Input BED file")
    parser.add_argument("--output_dir", "-o", required=True, help="Output directory for subfamily files and summaries")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Read BED entries
    subfamily_dict = defaultdict(list)

    with open(args.input, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue
            chrom, start, end, name, score, strand = row[:6]
            start = int(start)
            end = int(end)

            # Extract subfamily name (everything before first hyphen)
            subfamily = name.split("-")[0]

            # Parse fields from the name column
            name_parts = name.split("|")
            try:
                domain_class = name_parts[3]
                hmm = name_parts[2]
                confidence = name_parts[4]
                coverage = name_parts[5]
            except IndexError:
                domain_class = "NA"
                hmm = "NA"
                confidence = "NA"
                coverage = "NA"

            subfamily_dict[subfamily].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "coverage": coverage,
                "strand": strand,
                "domain_class": domain_class,
                "hmm": hmm,
                "confidence": confidence
            })

    # Write outputs
    for subfamily, entries in subfamily_dict.items():
        # BED file
        bed_path = os.path.join(args.output_dir, f"{subfamily}.bed")
        with open(bed_path, "w", newline="") as bed_f:
            writer = csv.writer(bed_f, delimiter="\t")
            for e in entries:
                writer.writerow([e["chrom"], e["start"], e["end"], e["name"], e["coverage"], e["strand"]])

        # Summary TSV
        summary_path = os.path.join(args.output_dir, f"{subfamily}_summary.tsv")
        with open(summary_path, "w", newline="") as sum_f:
            writer = csv.writer(sum_f, delimiter="\t")
            writer.writerow(["Subfamily", "Domain Class", "Chromosome", "Start", "End", "Coverage", "Confidence", "HMM"])
            for e in entries:
                writer.writerow([subfamily, e["domain_class"], e["chrom"], e["start"], e["end"],
                                 e["coverage"], e["confidence"], e["hmm"]])

    print(f"✅ Processed {len(subfamily_dict)} subfamilies into {args.output_dir}")

if __name__ == "__main__":
    main()
