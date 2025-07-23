#!/bin/bash

# Exit on error
set -e

# ------------------------------
# 14_run_interproscan.sh
# Run InterProScan on protein sequences
# ------------------------------

# Usage message
usage() {
  echo "Usage: $0 -i <input_fasta> -o <output_dir> -d <interproscan_dir> [-c <cpus>]"
  echo "  -i: Input protein FASTA file"
  echo "  -o: Output directory"
  echo "  -d: Path to InterProScan installation directory"
  echo "  -c: Number of CPUs to use (default: 8)"
  exit 1
}

# Default values
CPUS=8

# Parse arguments
while getopts "i:o:d:c:" opt; do
  case ${opt} in
    i ) INPUT_FASTA=$OPTARG ;;
    o ) OUTPUT_DIR=$OPTARG ;;
    d ) INTERPROSCAN_DIR=$OPTARG ;;
    c ) CPUS=$OPTARG ;;
    * ) usage ;;
  esac
done

# Check required arguments
if [[ -z "$INPUT_FASTA" || -z "$OUTPUT_DIR" || -z "$INTERPROSCAN_DIR" ]]; then
  usage
fi

# Check input file
if [[ ! -f "$INPUT_FASTA" ]]; then
  echo "❌ ERROR: Input FASTA file not found: $INPUT_FASTA"
  exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run InterProScan
echo "➡️ Running InterProScan on $INPUT_FASTA"
"$INTERPROSCAN_DIR/interproscan.sh" \
  -i "$INPUT_FASTA" \
  -dp \
  -cpu "$CPUS" \
  --output-dir "$OUTPUT_DIR"

echo "InterProScan completed: Results in $OUTPUT_DIR"
