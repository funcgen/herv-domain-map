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
    parser = argparse.ArgumentParser(description="Filter ENV domains from domain FASTA and run Phobius")
    parser.add_argument("--domain-fasta", required=True, help="Input FASTA file with domain protein sequences")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files (e.g., results/ENV)")
    parser.add_argument("--phobius-dir", required=True, help="Path to Phobius installation directory")
    parser.add_argument("--mode", choices=["short", "gff"], default=None, help="Phobius output format: short or gff")
    args = parser.parse_args()

    # Output files
    env_fasta = f"{args.output_prefix}_ENV_domains.fasta"
    phobius_out = f"{args.output_prefix}_phobius.txt"
    log_file = f"{args.output_prefix}_phobius.log"

    # Step 1: Extract ENV domains
    run_cmd([
        "python", "13_extract_env_sequences.py",
        "-i", args.domain_fasta,
        "-o", env_fasta
    ], "Filter ENV domain sequences from input FASTA")

    # Step 2: Run Phobius on ENV domains
    phobius_script = "./14_phobius_run_tool.sh"
    cmd = [phobius_script, "-i", env_fasta, "-o", phobius_out, "-d", args.phobius_dir]
    if args.mode == "short":
        cmd.append("--short")
    elif args.mode == "gff":
        cmd.append("--gff")

    run_cmd(cmd, "Run Phobius on ENV domain sequences", log_file=log_file)

    print("\nENV Domain Phobius Pipeline Completed!")

if __name__ == "__main__":
    main()
