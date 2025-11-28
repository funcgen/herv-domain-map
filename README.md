# HERV Domain Annotation and Conservation

This repository contains scripts to run the pipeline for identifying, annotating, and evaluating the conservation of protein-coding domains in **human endogenous retroviruses (HERVs)**. The project focuses on detecting key retroviral domains such as **Gag**, **Pol**, **Env**, and **accessory proteins** in predicted HERV open reading frames (ORFs), using HMM-based models and motif detection.
<p align="center">
 <img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/607eafc3-0462-4868-bd88-e257619f24f5" />
</p>


## 🧬 Overview

Although most HERVs are considered non-functional, a subset retain conserved protein-coding domains that may influence host cellular processes or pathology. This project aims to systematically annotate these domains and assess their structural conservation and potential functional relevance.

## 🚀 Features

- Domain annotation using **HMMER** with profiles from **GyDB**
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
├── scripts/
 ├── 01_generate_ERVs_bed.py
 ├── 02_merge_ERVs_annotations.py
 ├── 03_extract_fasta_from_bed.py
 ├── 04_fix_orfs_headers.py
 ├── 05_filter_hmmscan_hits.py
 ├── 06_map_domains_to_bed.py
 ├── 07_extract_domain_seq_from_orf.py
 ├── 08_analyse_alignments.py
 ├── 09_split_by_subfamily_and_summarize.py
 ├── 10_generate_erv_functionality_summary.py
 ├── 11_run_interproscan.sh
 ├── 12_find_conserved_residues.py
 ├── 13_extract_env_sequences.py
 ├── 14_phobius_run_tool.sh
 ├── run_erv_annotation_pipeline.py
 ├── run_erv_analysis_pipeline.py
 ├── run_erv_env_phobius_pipeline.py
├── README.md
├── LICENSE
```
## 📋 Example Usage

### Step 1: Run the Domain Annotation Pipeline
This step detects retroviral domains in ORFs using HMMER and GyDB profiles and filters domain hits.
For this step, the pipeline requires the files gydb_domains_classification.tsv and combined_gydb.hmm. The former is a TSV file that maps each HMM profile to its corresponding functional class (GAG, POL, ENV, Accessory, or Other), while the latter is the combined HMM profile database used by hmmscan. Both files are available in the Zenodo repository:  
- [gydb_domains_classification.tsv](https://zenodo.org/api/records/17662456/draft/files/gydb_domains_classification.tsv/content)
- [combined_gydb.hmm](https://zenodo.org/api/records/17662456/draft/files/combined_gydb.hmm/content)

```bash
python run_erv_annotation_pipeline.py \
  -r ../GRCh38.primary_assembly.genome.fa.out \
  -g ../GRCh38.primary_assembly.genome.fa \
  --orf-mode viral \
  -d /path/to/gydb_domains_classification.tsv \
  --min-coverage 0 \
  --hmm-db /path/to/combined_gydb.hmm \
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
  --interproscan-dir /path/to/interproscan \
  --cpus 32
```
InterProScan can be downloaded by following the official installation instructions available in the [InterProScan documentation](https://interproscan-docsdev.readthedocs.io/en/latest/HowToDownload.html).


### Step 3: Run Env Transmembrane Prediction (Phobius)
This step analyzes Env domain sequences for transmembrane features using Phobius.

```bash
python run_erv_env_phobius_pipeline.py \
  --domain-fasta results/ERV_GyDB_analysis_domains.fasta \
  --output-prefix results/ENV \
  --phobius-dir /path/to/phobius
```
Phobius can be downloaded from the official website at this [link](https://software.sbc.su.se/phobius.html).


### Step 4: Analysis of the Results (optional)
To reproduce the downstream analysis and generate the figures used in our manuscript, you can run the R script: `analysis_hervs_article.R`


## 📦 Dataset

This repository contains the scripts used to generate and analyze the dataset described in our study.

The full dataset—including annotated ORFs, domain predictions, alignment scores, and functional annotations—is available on Zenodo:

🔗 [https://doi.org/10.5281/zenodo.16318928](https://doi.org/10.5281/zenodo.16318928)

Please cite our article (see below) and the corresponding dataset DOI if you reuse the data in your own work.  


## 🧠 Citation  
If you use this resource in your research, please cite:  

Montserrat-Ayuso, T., & Esteve-Codina, A. (2025). *A comprehensive annotation of conserved protein domains in human endogenous retroviruses*. bioRxiv. https://doi.org/10.1101/2025.07.25.666750  

## 📎 License  
MIT License. See LICENSE file for details.
