#!/usr/bin/env python3

import argparse
from lxml import etree
import pandas as pd

def parse_interproscan_xml_grouped(xml_path):
    ns = {'ipr': 'https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas'}
    tree = etree.parse(xml_path)
    root = tree.getroot()
    records = []

    for idx, protein in enumerate(root.findall("ipr:protein", namespaces=ns)):
        xrefs = protein.findall("ipr:xref", namespaces=ns)
        xref_ids = [x.attrib.get("id", f"protein_{idx}") for x in xrefs]
        
        for rps_match in protein.findall("ipr:matches/ipr:rpsblast-match", namespaces=ns):
            domain = rps_match.find("ipr:signature", namespaces=ns)
            domain_desc = domain.attrib.get("desc", "unknown") if domain is not None else "unknown"

            for site in rps_match.findall(".//ipr:rpsblast-site", namespaces=ns):
                site_desc = site.attrib.get("description", "unknown")
                residue_list = [
                    f"{loc.attrib['residue']}{loc.attrib['start']}"
                    for loc in site.findall("ipr:site-locations/ipr:site-location", namespaces=ns)
                    if loc.attrib.get("residue") and loc.attrib.get("start")
                ]

                for protein_id in xref_ids:
                    if residue_list:
                        records.append({
                            "Protein": protein_id,
                            "Domain": domain_desc,
                            "Site Description": site_desc,
                            "Conserved Residues": ";".join(residue_list)
                        })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="Extract conserved residues from InterProScan XML")
    parser.add_argument("-i", "--input", required=True, help="Input InterProScan XML file")
    parser.add_argument("-o", "--output", default="conserved_residues_compact.tsv", help="Output TSV file name")
    args = parser.parse_args()

    df = parse_interproscan_xml_grouped(args.input)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[✓] Output written to: {args.output} ({len(df)} rows)")

if __name__ == "__main__":
    main()
