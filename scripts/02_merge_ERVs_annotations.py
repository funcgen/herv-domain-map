import sys

# You can change the allowed gap here (in base pairs)
ALLOWED_GAP = 100

def parse_bed_line(line):
    """
    Parse a BED line into a dictionary with fields including:
    - chrom, start, end, name, score, strand, subfamily (parsed from name)
    """
    fields = line.strip().split('\t')
    chrom, start, end, name, score, strand = fields
    subfamily = name.split('-int')[0] if '-int' in name else name  # Extract subfamily prefix
    return {
        'chrom': chrom,
        'start': int(start),
        'end': int(end),
        'name': name,
        'score': score,
        'strand': strand,
        'subfamily': subfamily
    }

def format_bed_line(entry):
    """
    Format a merged BED entry back to string, assigning a new standardized name.
    The start is incremented by 1 for 1-based naming in name field.
    """
    new_name = f"{entry['subfamily']}-int_pos_{entry['chrom']}_{entry['start']+1}_{entry['end']}_strand_{entry['strand']}"
    return f"{entry['chrom']}\t{entry['start']}\t{entry['end']}\t{new_name}\t{entry['score']}\t{entry['strand']}"

def merge_entries(entries):
    """
    Merge overlapping or closely spaced entries (within ALLOWED_GAP) of the same group:
    - Assumes entries are from same chromosome, strand, and subfamily.
    """
    if not entries:
        return []

    entries.sort(key=lambda x: x['start'])
    merged = [entries[0]]

    for current in entries[1:]:
        last = merged[-1]

        # Check if current entry is close enough to merge
        if (current['start'] <= last['end'] + ALLOWED_GAP and
            current['chrom'] == last['chrom'] and
            current['strand'] == last['strand'] and
            current['subfamily'] == last['subfamily']):
            
            # Extend the end coordinate of the last merged entry
            last['end'] = max(last['end'], current['end'])
        else:
            merged.append(current)

    return merged

def main(input_file, output_file):
    """
    Main function to:
    1. Parse BED input
    2. Group by (chrom, strand, subfamily)
    3. Merge nearby regions within each group
    4. Write to output BED
    """
    with open(input_file, 'r') as f:
        entries = [parse_bed_line(line) for line in f if line.strip()]

    grouped = {}
    for entry in entries:
        key = (entry['chrom'], entry['strand'], entry['subfamily'])
        grouped.setdefault(key, []).append(entry)

    merged_entries = []
    for group in grouped.values():
        merged_entries.extend(merge_entries(group))

    with open(output_file, 'w') as f:
        for entry in merged_entries:
            f.write(format_bed_line(entry) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_erv_bed.py <input_bed_file> <output_bed_file>")
    else:
        main(sys.argv[1], sys.argv[2])


