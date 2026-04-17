#!/usr/bin/env bash
set -euo pipefail

# =========================
# User-configurable values
# =========================
LOCATION="FossilLake"          # CHANGE AS NEEDED!
SAMPLE="30_1302373217_FL"      # CHANGE AS NEEDED!
THRESHOLD="0.01"               # Donut plot grouping threshold
TRUNC_LEN="250"                # DADA2 truncation length

# Input/output names
FASTQ_FILE="16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz"
WORKDIR="16S_rRNA_seq_${SAMPLE}"
MANIFEST="manifest.tsv"

# Classifier
CLASSIFIER_FILE="silva-138-99-515-806-nb-classifier.qza"
CLASSIFIER_URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"

# =========================
# Preflight checks
# =========================
command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH."; exit 1; }
command -v biom >/dev/null 2>&1 || { echo "Error: biom not found in PATH."; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python not found in PATH."; exit 1; }
command -v realpath >/dev/null 2>&1 || { echo "Error: realpath not found in PATH."; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "Error: wget not found in PATH."; exit 1; }

[[ -f "$FASTQ_FILE" ]] || { echo "Error: input FASTQ not found: $FASTQ_FILE"; exit 1; }
[[ -f "make_donut_graph.py" ]] || { echo "Error: make_donut_graph.py not found in current directory."; exit 1; }

# =========================
# Preparatory steps
# =========================
mkdir -p "$WORKDIR"
mv "$FASTQ_FILE" "$WORKDIR/"
cp make_donut_graph.py "$WORKDIR/"
cd "$WORKDIR"

printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"
printf "%s\t%s\n" "$LOCATION" "$(realpath "$FASTQ_FILE")" >> "$MANIFEST"

# =========================
# Import FASTQ, assess quality, denoise
# =========================
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path demux-single-end.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data demux-single-end.qza \
  --o-visualization demux-summary.qzv

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-single-end.qza \
  --p-trunc-len "$TRUNC_LEN" \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

# =========================
# Summarize feature table
# =========================
qiime feature-table summarize \
  --i-table table.qza \
  --o-feature-frequencies feature-frequencies.qza \
  --o-sample-frequencies sample-frequencies.qza \
  --o-summary table-summary.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# =========================
# Assign taxonomy using SILVA
# =========================
wget -O "$CLASSIFIER_FILE" "$CLASSIFIER_URL"

qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER_FILE" \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# =========================
# Generate taxonomy outputs
# =========================
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file "$MANIFEST" \
  --o-visualization taxa-bar-plots.qzv

qiime tools export \
  --input-path taxonomy.qza \
  --output-path exported-taxonomy

qiime tools export \
  --input-path table.qza \
  --output-path exported-table

# =========================
# Combine sequences with taxonomy
# =========================
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --o-visualization rep-seqs-with-taxonomy.qzv

# =========================
# Make genus-level abundance table
# =========================
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

qiime feature-table summarize \
  --i-table genus-table.qza \
  --o-feature-frequencies genus-feature-frequencies.qza \
  --o-sample-frequencies genus-sample-frequencies.qza \
  --o-summary genus-table-summary.qzv

qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus-table

biom convert \
  -i exported-genus-table/feature-table.biom \
  -o genus-table.tsv \
  --to-tsv

# =========================
# Generate donut graph
# =========================
python make_donut_graph.py \
  --manifest "$MANIFEST" \
  --table genus-table.tsv \
  --threshold "$THRESHOLD"

echo "Done. Results are in: $(pwd)"
