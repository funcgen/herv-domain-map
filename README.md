# HERV Domain Annotation and Conservation

This repository contains tools and pipelines for identifying, annotating, and evaluating the conservation of protein-coding domains in **human endogenous retroviruses (HERVs)**. The project focuses on detecting key retroviral domains such as **Gag**, **Pol**, **Env**, and **accessory proteins** in predicted HERV open reading frames (ORFs), using HMM-based models and motif detection.
<p align="center">
 <img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/607eafc3-0462-4868-bd88-e257619f24f5" />
</p>


## 🧬 Overview

Although most HERVs are considered non-functional, a subset retain conserved protein-coding domains that may influence host cellular processes or pathology. This project aims to systematically annotate these domains and assess their structural conservation and potential functional relevance.

## 🚀 Features

- Domain annotation using **HMMER** with profiles from **GyDB**.
- ORF-domain mapping and filtering (best hit per domain class)
- BED output for domain coordinates
- Coverage scoring and motif search
- Summary tables of domain presence per HERV subfamily
- Env-domains **Phobius** analysis for transmembrane domain prediction

## 🛠️ Dependencies

Conda environment:
``` 
name: ervscan
channels:
  - conda-forge
  - bioconda
  - defaults
  - r
dependencies:
  - bedtools
  - python=3.8
  - samtools
  - bioconda::pyfaidx
  - bioconda::hmmer
  - bioconda::emboss
  - bioconda::seqkit
  - pandas
  - bioconda::meme
prefix: /home/groups/funcgen/tmontser/miniconda3/envs/ervscan
``` 
Also: 
- [Phobius](https://phobius.sbc.su.se/)
- [InterProScan](https://www.ebi.ac.uk/interpro/download/InterProScan/)

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

### Step 1: Run the Domain Annotation Pipeline
This step detects retroviral domains in ORFs using HMMER and GyDB profiles and filters domain hits.

```bash
python run_erv_annotation_pipeline.py \
  -r ../GRCh38.primary_assembly.genome.fa.out \
  -g ../GRCh38.primary_assembly.genome.fa \
  --orf-mode viral \
  -d hmm_profiles/gydb_domains_classification.tsv \
  --min-coverage 0 \
  --hmm-db hmm_profiles/combined_gydb.hmm \
  --output-prefix results/ERV_GyDB_v4
```
### Step 2: Run the Domain Analysis Pipeline
This step performs a basic analysis of the results, and optionally integrates InterProScan annotations.

```bash
python run_erv_analysis_pipeline.py \
  --filtered-tsv results/ERV_GyDB_v4_filtered.tsv \
  --hmm-db hmm_profiles/combined_gydb.hmm \
  --orfs-fasta results/ERV_GyDB_v4_orfs_fixed.fasta \
  --bed results/ERV_GyDB_v4_domains.bed \
  --sto-dir results/alignments \
  --output-prefix results/ERV_GyDB_analysis \
  --interproscan-dir /home/groups/funcgen/tmontser/software/interproscan/interproscan-5.75-106.0 \
  --cpus 32
```

### Step 3: Run Env Transmembrane Prediction (Phobius)
This step analyzes Env domain sequences for transmembrane features using Phobius.

```bash
python run_erv_env_phobius_pipeline.py \
  --domain-fasta results/ERV_GyDB_analysis_domains.fasta \
  --output-prefix results/ENV \
  --phobius-dir /home/groups/funcgen/tmontser/software/phobius/phobius101_linux/phobius
```

## 📦 Dataset

This repository contains the scripts and pipelines used to generate and analyze the dataset described in our study.

The full dataset—including annotated ORFs, domain predictions, alignment scores, and functional annotations—is available on Zenodo:

🔗 [https://doi.org/10.5281/zenodo.16318928](https://doi.org/10.5281/zenodo.16318928)

Please cite our article (see below) and the corresponding dataset DOI if you reuse the data in your own work.  


## 🧠 Citation  
If you use this resource in your research, please cite:  

Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). *A comprehensive annotation of conserved protein domains in human endogenous retroviruses*. bioRxiv. https://doi.org/XXX

## 📎 License  
MIT License. See LICENSE file for details.
