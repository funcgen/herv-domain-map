#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os

def run_cmd(cmd, desc, log_file=None):
    print(f"\n{desc}\n{' '.join(cmd)}")
    if log_file:
        with open(log_file, "w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    else:
        result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"Step failed: {desc}")
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(description="ERV Domain Analysis Pipeline")
    parser.add_argument("--filtered-tsv", required=True, help="Filtered HMM hits TSV (from previous pipeline)")
    parser.add_argument("--hmm-db", required=True, help="HMM database file used for hmmscan and alignment")
    parser.add_argument("--orfs-fasta", required=True, help="ORF protein FASTA file")
    parser.add_argument("--bed", required=True, help="BED file (domains mapped to genome)")
    parser.add_argument("--sto-dir", required=True, help="Directory containing .sto alignments")
    parser.add_argument("--output-prefix", required=True, help="Output prefix for results (e.g., results/ERV_analysis)")
    parser.add_argument("--interproscan-dir", required=True, help="Path to InterProScan installation folder")
    parser.add_argument("--interproscan-output", required=False, default=None, help="Optional output folder for InterProScan (default: <output-prefix>_interproscan)")
    parser.add_argument("--skip-interproscan", action="store_true", help="Skip InterProScan step")
    parser.add_argument("--cpus", type=int, default=32, help="Number of CPUs to use for InterProScan (default: 32)")
    parser.add_argument("--skip-interproscan-parse", action="store_true", help="Skip parsing InterProScan XML output")


    args = parser.parse_args()

    # Define outputs
    domain_fasta = f"{args.output_prefix}_domains.fasta"
    alignment_summary = f"{args.output_prefix}_alignment_summary.tsv"
    subfamily_dir = f"{args.output_prefix}_subfamily"
    functionality_summary = f"{args.output_prefix}_functionality_summary.tsv"
    interproscan_output = args.interproscan_output or f"{args.output_prefix}_interproscan"
    interproscan_summary = f"{args.output_prefix}_conserved_residues.tsv"

    os.makedirs(subfamily_dir, exist_ok=True)

    # Step 1: Extract Domain Sequences
    run_cmd([
        "python", "07_extract_domain_seq_from_orf.py",
        "-t", args.filtered_tsv,
        "-f", args.orfs_fasta,
        "-o", domain_fasta
    ], "Extract Domain Sequences for Alignment")

    # Step 2: Align Domain Sequences
    align_log = f"{args.output_prefix}_align.log"
    run_cmd([
        "bash", "9_align_aa_sequences.sh", domain_fasta, args.hmm_db, args.sto_dir
    ], "Align Domain Sequences (hmmalign)", log_file=align_log)

    # Step 3: Analyze Alignments
    run_cmd([
        "python", "08_analyse_alignments.py",
        "--input", args.sto_dir,
        "--output", alignment_summary
    ], "Analyze Alignment Coverage")

    # Step 4: Split by Subfamily and Summarize
    run_cmd([
        "python", "09_split_by_subfamily_and_summarize.py",
        "--input", args.bed,
        "--output_dir", subfamily_dir
    ], "Split by Subfamily and Generate Summaries")

    # Step 5: Generate ERV Functionality Summary
    run_cmd([
        "python", "10_generate_erv_functionality_summary.py",
        "--input", args.bed,
        "--output", functionality_summary
    ], "Generate ERV Functionality Summary")

    # Step 6: Run InterProScan (Optional)
    if not args.skip_interproscan:
        run_cmd([
            "bash", "11_run_interproscan.sh",
            "-i", domain_fasta,
            "-o", interproscan_output,
            "-d", args.interproscan_dir,
            "-c", str(args.cpus)
        ], "Run InterProScan")

    # Step 7: Parse InterProScan XML for conserved residues
    interproscan_xml = os.path.join(interproscan_output, os.path.basename(domain_fasta) + ".xml")
    if not args.skip_interproscan and not args.skip_interproscan_parse:
        run_cmd([
            "python", "12_find_conserved_residues.py",
            "-i", interproscan_xml,
            "-o", interproscan_summary
        ], "Parse InterProScan XML to extract conserved residues")



    print("\nERV Domain Analysis Pipeline Completed!")


if __name__ == "__main__":
    main()
