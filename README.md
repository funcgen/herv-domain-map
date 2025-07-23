# HERV Domain Annotation and Conservation

This repository contains tools and pipelines for identifying, annotating, and evaluating the conservation of protein-coding domains in **human endogenous retroviruses (HERVs)**. The project focuses on detecting key retroviral domains such as **Gag**, **Pol**, **Env**, and **accessory proteins** in predicted HERV open reading frames (ORFs), using HMM-based models and motif detection.

## 🧬 Overview

Although most HERVs are considered non-functional, a subset retain conserved protein-coding domains that may influence host cellular processes or pathology. This project aims to systematically annotate these domains and assess their sequence conservation and potential functional relevance.

## 🚀 Features

- Domain annotation using **HMMER** with profiles from **GyDB**.
- ORF-domain mapping and filtering (best hit per gene class per ORF)
- BED output for domain coordinates
- Coverage scoring and motif search (e.g., YXDD, PPXY, cysteine motifs)
- Summary tables of domain presence per HERV subfamily
- Optional ENV extraction and **Phobius** analysis for transmembrane domain prediction

## 🛠️ Dependencies

- Python 3.8+
- `pandas`, `biopython`, `argparse`
- [HMMER 3.x](http://hmmer.org/)
- [Phobius](https://phobius.sbc.su.se/) (optional)
- BEDTools (for some genomic operations)

## 📁 Repository Structure

```bash
.
├── data/                    # Input data: ORFs, HMM profiles, etc.
├── scripts/
│   ├── 01_filter_hmmscan_hits.py
│   ├── 02_map_domains_to_bed.py
│   ├── 03_score_alignment_and_motifs.py
│   ├── 16_extract_env_sequences.py
│   └── 17_phobius_run_tool.sh
├── results/                 # Output BED files, summary tables
├── README.md
└── requirements.txt
```
## 📋 Example Usage
 ```bash
# Step 1: Filter HMMER hits
python scripts/01_filter_hmmscan_hits.py --input raw_hits.tbl --output filtered_hits.tsv

# Step 2: Map domains to genome
python scripts/02_map_domains_to_bed.py --input filtered_hits.tsv --output domains.bed

# Step 3: Score alignment and detect motifs
python scripts/03_score_alignment_and_motifs.py --input alignments/ --output summary.tsv
```

🧠 Citation
If you use this resource in your research, please cite:

[Authors], "A Comprehensive Annotation of Conserved Protein Domains in Human Endogenous Retroviruses", [Journal, Year] (coming soon)

📎 License
MIT License. See LICENSE file for details.
