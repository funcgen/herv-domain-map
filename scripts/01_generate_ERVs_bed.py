#!/usr/bin/env python3
"""
Extract BED entries of internal ERV insertions
from a RepeatMasker .out file.

Internal regions are defined by:
- repeat_class in ERV / MaLR / Gypsy LTR classes
- repeat_name ending in -int / _int (case-insensitive) or _I (internal)
"""

import re
import argparse
from tqdm import tqdm


def is_internal_class(rep_class: str) -> bool:
    """
    Recognize RepeatMasker classes that contain HERV internal regions
    based on class patterns (no whitelist of names).
    """
    rep_class = rep_class.upper()

    return (
        rep_class.startswith("LTR/ERV")        # LTR/ERV1, LTR/ERVK, LTR/ERVL...
        or rep_class.startswith("LTR/ERVL-MALR")  # MaLR (MLT, THE1...)
        or rep_class.startswith("LTR/GYPSY")      # MamGyp internal
    )


def is_internal(name: str) -> bool:
    """
    Detect internal ERV regions based on RepeatMasker name patterns:

    - '-int' or '_int' at the end (MER4-int, HERVK13-int, MLT1J-int, ...)
    - '_I' at the end (HERV1_I, HERV4_I, ERV3-16A3_I, THE1_I)

    """
    if not name:
        return False

    n = name.strip()

    # 1. Standard internal naming: '-int' or '_int' (case-insensitive)
    if re.search(r'(?:-|_)int$', n, flags=re.IGNORECASE):
        return True

    # 2. Internal-only families using the '_I' suffix
    #    (HERV1_I, HERV4_I, ERV3-16A3_I, THE1_I, ...)
    #    Note: this does NOT catch '-I' because of the underscore.
    if re.search(r'_I$', n, flags=re.IGNORECASE):
        return True

    return False


def extract_internal_bed(rm_out, output_bed):
    entries = []

    with open(rm_out) as fh:
        for line in fh:
            # Skip comment lines or headers (typically start with SW or score)
            if line.strip() == "" or line.startswith("SW") or line.startswith("score"):
                continue

            fields = re.split(r"\s+", line.strip())
            if len(fields) < 11:
                continue  # skip malformed lines

            repeat_name  = fields[9]
            repeat_class = fields[10]

            # Filter by class (ERV / MaLR / Gypsy) and internal name pattern
            if not (is_internal_class(repeat_class) and is_internal(repeat_name)):
                continue

            divergence = fields[1]
            chrom      = fields[4]
            start      = int(fields[5]) - 1  # BED is 0-based
            end        = int(fields[6])
            strand     = '-' if fields[8].upper() == 'C' else '+'

            name = (
                f"{repeat_name}_pos_{chrom}_{start+1}_{end}"
                f"_strand_{strand}_divergence_{divergence}"
            )
            score = "0"

            entries.append([chrom, str(start), str(end), name, score, strand])

    # Write BED
    with open(output_bed, "w") as out:
        for e in entries:
            out.write("\t".join(e) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Extract BED of ERV internal sequences from RepeatMasker output.")
    parser.add_argument("--repeatmasker", "-r", required=True, help="RepeatMasker .out file")
    parser.add_argument("--output", "-o", required=True, help="Output BED file")
    args = parser.parse_args()

    extract_internal_bed(args.repeatmasker, args.output)


if __name__ == "__main__":
    main()
