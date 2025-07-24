#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Generate ERV subfamily functionality summary with quality labels from BED file.")
    parser.add_argument("--input", "-i", required=True, help="Input BED file (filtered, no overlaps)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV summary file")
    return parser.parse_args()

# Define a ranking for conservedness labels
LABEL_RANKS = {"Full-length": 1, "High": 2, "Moderate": 3, "Low": 4}

def get_better_label(label1, label2):
    if label1 is None:
        return label2
    if label2 is None:
        return label1
    return label1 if LABEL_RANKS.get(label1, 5) < LABEL_RANKS.get(label2, 5) else label2

def main():
    args = parse_args()

    summary = defaultdict(lambda: {
        "GAG": False, "GAG_quality": None,
        "POL": False, "POL_quality": None,
        "ENV": False, "ENV_quality": None,
        "ACCESSORY": False, "ACCESSORY_quality": None,
        "DOMAINS": set()
    })

    with open(args.input, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue
            chrom, start, end, name, score, strand = row[:6]
            subfamily = name.split("-")[0]

            parts = name.split("|")
            try:
                domain_class = parts[3].upper()
                domain_name = parts[2]
                conservedness = parts[4]
            except IndexError:
                domain_class = "NA"
                domain_name = "NA"
                conservedness = None

            # Update presence and quality
            if domain_class in ["GAG", "POL", "ENV", "ACCESSORY"]:
                summary[subfamily][domain_class] = True
                current_best = summary[subfamily][f"{domain_class}_quality"]
                summary[subfamily][f"{domain_class}_quality"] = get_better_label(current_best, conservedness)

            summary[subfamily]["DOMAINS"].add(domain_name)

    # Write summary TSV
    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        header = [
            "Subfamily", 
            "GAG", "GAG Quality",
            "POL", "POL Quality",
            "ENV", "ENV Quality",
            "ACCESSORY", "ACCESSORY Quality",
            "Potentially Functional",
            "Domains Detected"
        ]
        writer.writerow(header)

        for subfamily, info in sorted(summary.items()):
            potential_flag = "Yes" if info["GAG"] and info["POL"] and info["ENV"] else "No"

            row = [
                subfamily,
                "Yes" if info["GAG"] else "No", info["GAG_quality"] or "-",
                "Yes" if info["POL"] else "No", info["POL_quality"] or "-",
                "Yes" if info["ENV"] else "No", info["ENV_quality"] or "-",
                "Yes" if info["ACCESSORY"] else "No", info["ACCESSORY_quality"] or "-",
                potential_flag,
                "; ".join(sorted(info["DOMAINS"])) if info["DOMAINS"] else "-"
            ]
            writer.writerow(row)

    print(f"Summary with quality and potential functionality generated: {args.output}")

if __name__ == "__main__":
    main()
