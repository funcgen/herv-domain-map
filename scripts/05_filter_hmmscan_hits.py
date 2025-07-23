#!/usr/bin/env python3
import argparse
import re
import pandas as pd

def load_domain_classification(tsv_path):
    mapping = {}
    with open(tsv_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            name = fields[0].lower()
            gene_class = fields[1].lower()
            mapping[name] = gene_class
    return mapping

def classify_gene(hmm_name, domain_dict):
    hmm_name = hmm_name.lower()
    return domain_dict.get(hmm_name, "other")

def conservedness_label(cov):
    if cov >= 0.95:
        return "Full-length"
    elif cov >= 0.80:
        return "High"
    elif cov >= 0.60:
        return "Moderate"
    else:
        return "Low"

def extract_locus(query):
    """Remove '_<number>|...' from query"""
    return re.sub(r"_[0-9]+\|.*$", "", query)

def extract_subfamily(query):
    """Extract subfamily from query, stopping at '-int'"""
    return re.sub(r"-int.*", "", query)

def extract_domain_type(domain):
    """Get domain type from domain name (before first '_')"""
    return domain.split("_")[0]

def main():
    parser = argparse.ArgumentParser(description="Filter HMMscan results by E-value and coverage, and keep only the best hit per Locus and DomainType.")
    parser.add_argument("-i", "--input", required=True, help="Input .tbl file from hmmscan --tblout")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file (filtered and deduplicated)")
    parser.add_argument("--domain-classes", required=True, help="TSV file mapping HMM names to gene classes")
    parser.add_argument("--max-evalue", type=float, default=1e-5, help="Max full-sequence E-value (default 1e-5)")
    parser.add_argument("--max-domain-evalue", type=float, default=1e-5, help="Max domain E-value (default 1e-5)")
    parser.add_argument("--min-coverage", type=float, default=0, help="Min coverage (default 0)")
    args = parser.parse_args()

    domain_dict = load_domain_classification(args.domain_classes)

    filtered_rows = []

    with open(args.input) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = re.split(r'\s+', line.strip())
            try:
                target_name = fields[0]
                hmm_length = int(fields[2])
                query_name = fields[3]
                full_seq_evalue = float(fields[6])
                full_seq_score = float(fields[7])
                domain_evalue = float(fields[11])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
            except (IndexError, ValueError):
                continue

            aligned_length = hmm_to - hmm_from + 1
            coverage = round(aligned_length / hmm_length, 3)

            if (
                full_seq_evalue <= args.max_evalue and
                domain_evalue <= args.max_domain_evalue and
                coverage >= args.min_coverage
            ):
                gene = classify_gene(target_name, domain_dict)
                conservedness = conservedness_label(coverage)
                filtered_rows.append({
                    "Query": query_name,
                    "Gene": gene,
                    "Domain": target_name,
                    "Score": full_seq_score,
                    "E-value": full_seq_evalue,
                    "Domain_E-value": domain_evalue,
                    "Coverage": coverage,
                    "Conservedness": conservedness,
                    "Line": line.strip()
                })


    if not filtered_rows:
        print("No hits passed filtering criteria.")
        return

    # Convert to DataFrame
    df = pd.DataFrame(filtered_rows)

    # Extract additional columns
    df["Locus"] = df["Query"].apply(extract_locus)
    df["Subfamily"] = df["Query"].apply(extract_subfamily)
    df["DomainType"] = df["Domain"].apply(extract_domain_type)
    df["Name"] = df["Locus"] + "_" + df["DomainType"]

    # Keep only best score per (Locus, DomainType)
    df_best = (
        df.sort_values("Score", ascending=False)
          .drop_duplicates(subset=["Name", "DomainType"])
          .reset_index(drop=True)
    )

    # Reorder columns for output
    df_best = df_best[[
    "Query", "Gene", "Domain", "DomainType", "Score", "E-value",
    "Domain_E-value", "Coverage", "Conservedness", "Subfamily", "Locus", "Name", "Line"
]]

    # Write final output
    df_best.to_csv(args.output, sep="\t", index=False)

    print(f"Final filtered file with one hit per Locus and DomainType written to: {args.output}")

if __name__ == "__main__":
    main()

