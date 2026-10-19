#!/bin/bash
#SBATCH --job-name=herv_pipeline
#SBATCH --output=logs/analysis_complete_pipeline_%j.out
#SBATCH --error=logs/analysis_complete_pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --qos=normal
#SBATCH --time=06:00:00  # 2 days max runtime
#SBATCH --partition=general

# Example for a NEW genome reconstruction. For a domain-only database release,
# follow the frozen-ORF rescan commands in the release README instead.
set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Activate environment
echo "Activating environment."
source ~/.bashrc
conda activate te_preprocessing


echo "Running annotation pipeline."
python "$SCRIPT_DIR"/run_erv_annotation_pipeline.py \
  -r ../GRCh38.primary_assembly.genome.fa.out \
  -g ../GRCh38.primary_assembly.genome.fa \
  --orf-mode viral \
  -d hmm_profiles/gydb_domains_classification.tsv \
  --min-coverage 0 \
  --species human \
  --hmm-Z 316 \
  --profile-thresholds "$SCRIPT_DIR/../config/ERVL_profile_thresholds.tsv" \
  --primary-scope "$SCRIPT_DIR/../config/ERVL_primary_scope.json" \
  --hmm-db hmm_profiles/combined_gydb.hmm \
  --short-gap 200 \
  --long-gap 0 \
  --pad 10 \
  --require-rmsk \
  --merge-debug-tsv results/ERV_GyDB_v6_merge_debug.tsv \
  --output-prefix results/ERV_GyDB_v6


echo "Activating necessary modules."

module load Java/11.0.2
module load PCRE2

echo "Running analysis pipeline."
python "$SCRIPT_DIR"/run_erv_analysis_pipeline.py \
--filtered-tsv results/ERV_GyDB_v6_filtered.tsv \
--hmm-db hmm_profiles/combined_gydb.hmm \
--orfs-fasta results/ERV_GyDB_v6_orfs_fixed.fasta \
--bed results/ERV_GyDB_v6_domains.bed \
--sto-dir results/alignments \
--output-prefix results/ERV_GyDB_analysis \
--interproscan-dir /home/groups/funcgen/tmontser/software/interproscan/interproscan-5.75-106.0 \
--cpus 32


echo "Running Phobius pipeline."
python "$SCRIPT_DIR"/run_erv_env_phobius_pipeline.py \
--domain-fasta results/ERV_GyDB_analysis_domains.fasta \
--output-prefix results/ENV \
--phobius-dir /home/groups/funcgen/tmontser/software/phobius/phobius101_linux/phobius


echo "Done."
