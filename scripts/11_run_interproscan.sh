#!/bin/bash

# Exit on error
set -euo pipefail
INPUT_FASTA=""; OUTPUT_DIR=""; INTERPROSCAN_DIR=""

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
  echo "ERROR: Input FASTA file not found: $INPUT_FASTA"
  exit 2
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run InterProScan
echo "➡Running InterProScan on $INPUT_FASTA"
"$INTERPROSCAN_DIR/interproscan.sh" \
  -i "$INPUT_FASTA" \
  -f XML,JSON,TSV,GFF3 \
  -dp \
  -cpu "$CPUS" \
  --output-dir "$OUTPUT_DIR"

# Verify JSON completeness against both source sequences and XML before reporting success.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_BASENAME="$(basename -- "$INPUT_FASTA")"
"${PYTHON_BIN:-python3}" "$SCRIPT_DIR/16_validate_interproscan_json.py" \
  "$OUTPUT_DIR/$OUTPUT_BASENAME.json" \
  --fasta "$INPUT_FASTA" \
  --xml "$OUTPUT_DIR/$OUTPUT_BASENAME.xml" \
  --repair-omitted-root \
  --report "$OUTPUT_DIR/$OUTPUT_BASENAME.json.validation.json"

echo "InterProScan completed: Results in $OUTPUT_DIR"
