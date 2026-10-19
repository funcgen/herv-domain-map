#!/usr/bin/env python3

import argparse
from pathlib import Path

def parse_hmmalign_sto(sto_file):
    alignment = {}
    with open(sto_file) as f:
        for line in f:
            if line.startswith("#=GR") and "PP" in line:
                parts = line.strip().split()
                seqname = parts[1]
                pp = parts[3]
                alignment.setdefault(seqname, {"seq": "", "pp": ""})
                alignment[seqname]["pp"] += pp
            elif not line.startswith("#") and line.strip() and not line.startswith("//"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    seqname = parts[0]
                    seq = parts[1]
                    alignment.setdefault(seqname, {"seq": "", "pp": ""})
                    alignment[seqname]["seq"] += seq
    return alignment

def score_alignment(pp_string):
    # Count the original PP symbol set 5-9 and *. These are rounded posterior
    # bins, not exact P >= 0.5 and not evidence of protein functionality.
    high_conf = sum(1 for c in pp_string if c in ["*", "9", "8", "7", "6", "5"])
    total = len(pp_string)
    return high_conf / total if total > 0 else 0

def main():
    parser = argparse.ArgumentParser(description="Parse hmmalign alignments and score posterior support.")
    parser.add_argument("--input", required=True, help="Directory with .sto files")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    sto_files = list(Path(args.input).glob("*.sto"))
    if not sto_files:
        print("No .sto files found in input directory.")
        return

    with open(args.output, "w") as out_f:
        out_f.write("Sequence\tProfile\tSupportedAlignmentColumnFraction\tColumnSupportCategory\tResiduePosteriorSupportedFraction\n")
        for file in sto_files:
            alignment = parse_hmmalign_sto(file)
            for seqname, data in alignment.items():
                pp = data["pp"]
                coverage = score_alignment(pp)
                assert len(data['seq']) == len(pp)
                residue_pp = ''.join(c for aa,c in zip(data['seq'],pp) if aa not in '.-')
                residue_support = score_alignment(residue_pp)

                # Extract the exact profile from the right-hand metadata fields
                parts = seqname.split("|")
                domain_type = parts[-2] if len(parts) >= 4 else "Unknown"

                # Categories describe alignment-column posterior support only
                if coverage >= 0.85:
                    status = "High"
                elif coverage >= 0.70:
                    status = "Moderate"
                elif coverage >= 0.50:
                    status = "Partial"
                elif coverage >= 0.30:
                    status = "Weak"
                elif coverage > 0:
                    status = "Very Weak"
                else:
                    status = "No Alignment"

                out_f.write(f"{seqname}\t{domain_type}\t{coverage:.3f}\t{status}\t{residue_support:.3f}\n")

    print(f"Report written to {args.output}")

if __name__ == "__main__":
    main()
