#!/usr/bin/env python3

import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="Filter overlapping BED entries by domain class, strand, and coverage (no pandas).")
    parser.add_argument("--input", "-i", required=True, help="Input BED file")
    parser.add_argument("--output", "-o", required=True, help="Filtered BED output file")
    parser.add_argument("--log", "-l", required=True, help="Log file with removed entries")
    parser.add_argument("--overlap", "-x", type=float, default=80.0, help="Overlap threshold in percent (default: 80%%)")
    return parser.parse_args()

def compute_overlap(a_start, a_end, b_start, b_end):
    overlap_start = max(a_start, b_start)
    overlap_end = min(a_end, b_end)
    return max(0, overlap_end - overlap_start)

def main():
    args = parse_args()
    overlap_thresh = args.overlap / 100.0

    # Read BED entries
    entries = []
    with open(args.input, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 6:
                continue  # Skip incomplete rows
            chrom, start, end, name, score, strand = row[:6]
            start = int(start)
            end = int(end)
            domain_class = name.split("|")[3]
            entries.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "score": score,
                "strand": strand,
                "domain_class": domain_class
            })

    # Group by (chrom, strand, domain_class)
    groups = {}
    for entry in entries:
        key = (entry["chrom"], entry["strand"], entry["domain_class"])
        groups.setdefault(key, []).append(entry)

    kept = []
    removed = []

    # Process each group
    for key, group in groups.items():
        group.sort(key=lambda x: (x["start"], x["end"]))
        keep_flags = [True] * len(group)

        for i in range(len(group)):
            if not keep_flags[i]:
                continue
            a = group[i]
            a_len = a["end"] - a["start"]

            for j in range(i + 1, len(group)):
                if not keep_flags[j]:
                    continue
                b = group[j]
                b_len = b["end"] - b["start"]

                overlap_len = compute_overlap(a["start"], a["end"], b["start"], b["end"])
                if overlap_len == 0:
                    continue

                smaller_len = min(a_len, b_len)
                overlap_pct = overlap_len / smaller_len

                if overlap_pct >= overlap_thresh:
                    if a_len >= b_len:
                        keep_flags[j] = False
                        removed.append(b)
                    else:
                        keep_flags[i] = False
                        removed.append(a)
                        break

        for k, entry in enumerate(group):
            if keep_flags[k]:
                kept.append(entry)

    # Write output
    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        for entry in kept:
            writer.writerow([entry["chrom"], entry["start"], entry["end"], entry["name"], entry["score"], entry["strand"]])

    # Write log
    with open(args.log, "w", newline="") as log_f:
        writer = csv.writer(log_f, delimiter="\t")
        if removed:
            writer.writerow(["chrom", "start", "end", "name", "score", "strand"])  # header
            for entry in removed:
                writer.writerow([entry["chrom"], entry["start"], entry["end"], entry["name"], entry["score"], entry["strand"]])
        else:
            log_f.write("No entries removed.\n")

    print(f"Filtering complete. Kept: {len(kept)} | Removed: {len(removed)}")

if __name__ == "__main__":
    main()
