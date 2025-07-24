#!/usr/bin/env python3
import argparse

def fix_header(header):
    """
    Convert 'name [start - end]' into 'name|start-end'
    """
    if '[' in header and ']' in header:
        try:
            name, coords = header.split('[', 1)
            coords = coords.strip('] \n').replace(' ', '')
            coords = coords.replace('-', '-')
            return f"{name.strip()}|{coords}"
        except Exception as e:
            print(f"Warning: Could not parse header: {header.strip()}")
            return header.strip()
    else:
        return header.strip()

def main():
    parser = argparse.ArgumentParser(description="Fix ORF FASTA headers for hmmscan input.")
    parser.add_argument("-i", "--input", required=True, help="Input ORF FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with fixed headers")
    args = parser.parse_args()

    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(">" + fix_header(line[1:]) + "\n")
            else:
                fout.write(line)

    print(f"Fixed headers written to: {args.output}")

if __name__ == "__main__":
    main()
