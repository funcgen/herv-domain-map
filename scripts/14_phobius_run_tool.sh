#!/bin/bash

# Tool to run Phobius on a FASTA file
# Usage: ./phobius_run_tool.sh -i input.fasta -o output.txt [-d /path/to/phobius] [--short|--gff]

# Default values
PHOBIUS_DIR="$(dirname "$0")"
MODE=""
INPUT=""
OUTPUT=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -d|--dir)
            PHOBIUS_DIR="$2"
            shift 2
            ;;
        --short)
            MODE="-short"
            shift
            ;;
        --gff)
            MODE="-gff"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 -i input.fasta -o output.txt [-d /path/to/phobius] [--short|--gff]"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h for help."
            exit 1
            ;;
    esac
done

# Validate input
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Error: Input and output must be specified."
    echo "Usage: $0 -i input.fasta -o output.txt [-d /path/to/phobius] [--short|--gff]"
    exit 1
fi

if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file '$INPUT' does not exist."
    exit 1
fi

# Validate input
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Error: Input and output must be specified."
    echo "Usage: $0 -i input.fasta -o output.txt [-d /path/to/phobius] [--short|--gff]"
    exit 1
fi

# Convert to absolute paths
INPUT="$(readlink -f "$INPUT")"
OUTPUT="$(readlink -f "$OUTPUT")"


# Change to Phobius directory to use local decodeanhmm, model, options
cd "$PHOBIUS_DIR" || {
    echo "Error: Cannot cd to Phobius directory '$PHOBIUS_DIR'"
    exit 1
}

# Run Phobius
echo "Running Phobius on $INPUT ..."
perl phobius.pl $MODE "$INPUT" > "$OUTPUT"

if [[ $? -eq 0 ]]; then
    echo "Phobius finished successfully."
    echo "Results written to: $OUTPUT"
else
    echo "Phobius failed to execute."
    exit 1
fi
