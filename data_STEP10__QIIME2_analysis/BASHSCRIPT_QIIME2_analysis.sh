#!/usr/bin/env bash
set -euo pipefail

# Load conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1
 
# =========================
# User-configurable values
# =========================
LOCATION="HorseThief_Reservoir" # CHANGE AS NEEDED!
SAMPLE="30_1326789214_HTF"      # CHANGE AS NEEDED!
THRESHOLD="0.01"               # Donut plot grouping threshold
TRUNC_LEN="250"                # DADA2 truncation length

# Input/output names
ORIG_DIR="$(pwd)"
FASTQ_FILE="16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz"
WORKDIR="16S_rRNA_seq_${SAMPLE}"
MANIFEST="manifest.tsv"

# Classifier
CLASSIFIER_FILE_READS="silva-138-99-seqs-515-806.qza"
CLASSIFIER_FILE_TAX="silva-138-99-tax-515-806.qza"

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
mv "$CLASSIFIER_FILE_READS" "$WORKDIR/"
mv "$CLASSIFIER_FILE_TAX" "$WORKDIR/"
mv make_donut_graph.py "$WORKDIR/"
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
qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs.qza \
  --i-reference-reads  "$CLASSIFIER_FILE_READS" \
  --i-reference-taxonomy "$CLASSIFIER_FILE_TAX" \
  --p-perc-identity 0.97 \
  --p-query-cov 0.8 \
  --o-classification taxonomy.qza \
  --o-search-results blast-results.qza

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

# Convert BIOM table to TSV (prefixed filename)
GENUS_TSV="${WORKDIR}_genus-table.tsv"

biom convert \
  -i exported-genus-table/feature-table.biom \
  -o "$GENUS_TSV" \
  --to-tsv

# =========================
# Generate donut graph
# =========================
python make_donut_graph.py \
  --manifest "$MANIFEST" \
  --table "$GENUS_TSV" \
  --threshold "$THRESHOLD"

# Rename donut outputs (prepend WORKDIR)
PNG_OUT="${WORKDIR}_donut.png"
SVG_OUT="${WORKDIR}_donut.svg"

mv *_donut.png "$PNG_OUT"
mv *_donut.svg "$SVG_OUT"

# Copy outputs back to original directory
cp "$GENUS_TSV" "$ORIG_DIR/"
cp "$PNG_OUT" "$ORIG_DIR/"
cp "$SVG_OUT" "$ORIG_DIR/"
