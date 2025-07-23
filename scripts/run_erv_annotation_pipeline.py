#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def run_cmd(cmd, desc):
    print(f"\n {desc}\n{' '.join(cmd)}")
    result = subprocess.run(cmd, shell=False)
    if result.returncode != 0:
        print(f"Step failed: {desc}")
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(description="ERV Protein Annotation Pipeline")
    parser.add_argument("--repeatmasker", "-r", required=True, help="RepeatMasker .out file")
    parser.add_argument("--genome", "-g", required=True, help="Reference genome FASTA (indexed with pyfaidx)")
    parser.add_argument(
        "--orf-mode", choices=["canonical", "viral"], default="canonical",
        help="ORF detection mode: 'canonical' (AUG-start) or 'viral' (stop-to-stop)"
    )
    parser.add_argument("--domain-classes", "-d", required=True, help="HMM domain to gene class TSV")
    parser.add_argument("--hmm-db", required=True, help="HMM database file for hmmscan")
    parser.add_argument("--output-prefix", "-o", required=True, help="Output prefix (e.g., results/ERV)")

    parser.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum allowed E-value for hmmscan hits (default: 1e-5)")
    parser.add_argument("--min-coverage", type=float, default=0, help="Minimum HMM profile coverage (fraction between 0–1) to keep hits (default: 0)")
    parser.add_argument("--max-domain-evalue", type=float, default=1e-5, help="Maximum allowed E-value at the domain level (default: 1e-5)")

    parser.add_argument("--orf-minsize", type=int, default=180, help="Minimum ORF size in nucleotides for getorf (default: 180)")
    args = parser.parse_args()

    # Intermediate files
    raw_bed = f"{args.output_prefix}_raw.bed"
    merged_bed = f"{args.output_prefix}_merged.bed"
    internal_fasta = f"{args.output_prefix}_internal.fasta"
    orfs_fasta = f"{args.output_prefix}_orfs.fasta"
    orfs_fixed_fasta = f"{args.output_prefix}_orfs_fixed.fasta"
    hmmscan_tbl = f"{args.output_prefix}_hmmscan.tbl"
    filtered_tsv = f"{args.output_prefix}_filtered.tsv"
    bed_output = f"{args.output_prefix}_domains.bed"
    bed_fasta = f"{args.output_prefix}_domains.fasta"

    # Step 1: Generate BED file from RepeatMasker .out
    run_cmd([
        "python", "01_generate_ERVs_bed.py",
        "-r", args.repeatmasker,
        "-o", raw_bed
    ], "Generate raw BED from RepeatMasker")

    # Step 2: Merge overlapping ERV-internal entries
    run_cmd([
        "python", "02_merge_ERVs_annotations.py",
        raw_bed,
        merged_bed
    ], "Merge overlapping ERV internal regions")

    # Step 3: Extract sequences from merged BED
    run_cmd([
        "python", "03_extract_fasta_from_bed.py",
        "-b", merged_bed,
        "-g", args.genome,
        "-o", internal_fasta
    ], "Extract internal sequences from merged BED")

    # Step 4: ORF Prediction
    orf_cmd = [
        "getorf",
        "-sequence", internal_fasta,
        "-outseq", orfs_fasta,
        "-minsize", str(args.orf_minsize),
        "-reverse", "N"
    ]

    if args.orf_mode == "canonical":
        orf_cmd.extend(["-find", "1", "-methionine", "Y"])
    else:
        orf_cmd.extend(["-find", "0", "-methionine", "N"])

    run_cmd(orf_cmd, f"Run EMBOSS getorf on internal sequences ({args.orf_mode} mode)")

    # Step 5: Fix ORF Headers
    run_cmd([
        "python", "04_fix_orfs_headers.py",
        "-i", orfs_fasta,
        "-o", orfs_fixed_fasta
    ], "Fix ORF Headers for HMMER")

    # Step 6: HMM Scan (with logging)
    hmmscan_log = f"{args.output_prefix}_hmmscan.log"
    with open(hmmscan_log, "w") as log_file:
        print(f"\n Run hmmscan on ORFs\nLogging to {hmmscan_log}")
        result = subprocess.run([
            "hmmscan", "--domtblout", hmmscan_tbl, args.hmm_db, orfs_fixed_fasta
        ], stdout=log_file, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            print("HMMscan failed")
            sys.exit(result.returncode)

    # Step 7: Filter HMM Hits
    run_cmd([
        "python", "05_filter_hmmscan_hits.py",
        "-i", hmmscan_tbl,
        "-o", filtered_tsv,
        "--domain-classes", args.domain_classes,
        "--max-evalue", str(args.max_evalue),
        "--max-domain-evalue", str(args.max_domain_evalue),
        "--min-coverage", str(args.min_coverage)
    ], "Filter HMMscan Hits")

    # Step 8: Map Domains to BED
    run_cmd([
        "python", "06_map_domains_to_bed.py",
        "-i", filtered_tsv,
        "-o", bed_output
    ], "Map Domains to BED")

    # Step 9: Extract BED Sequences
    run_cmd([
        "bedtools", "getfasta",
        "-fi", args.genome,
        "-bed", bed_output,
        "-fo", bed_fasta,
        "-name"
    ], "Extract Sequences from BED (bedtools getfasta)")

    print(f"\nPipeline completed!\n Final BED: {bed_output}\n Domain FASTA: {bed_fasta}")

if __name__ == "__main__":
    main()
