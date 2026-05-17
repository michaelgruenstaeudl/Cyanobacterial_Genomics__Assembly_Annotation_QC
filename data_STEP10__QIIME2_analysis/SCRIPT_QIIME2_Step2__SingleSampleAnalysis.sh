#!/usr/bin/env bash
set -euo pipefail

SAMPLE="${1:?Usage: bash BASHSCRIPT_QIIME2_single_sample.sh SAMPLE_NAME}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1

THRESHOLD="0.01"
TRUNC_LEN="250"
MIN_FREQUENCY="100"

ORIG_DIR="$(pwd)"
FASTQ_FILE="16S_rRNA_seq_${SAMPLE}.assembled.fastq.gz"
WORKDIR="16S_rRNA_seq_${SAMPLE}"
MANIFEST="manifest_${SAMPLE}.tsv"
SHARED_ASV_IDS="shared_ASV_ids.tsv"

CLASSIFIER_FILE_READS="silva-138-99-seqs-515-806.qza"
CLASSIFIER_FILE_TAX="silva-138-99-tax-515-806.qza"

command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH."; exit 1; }
command -v biom >/dev/null 2>&1 || { echo "Error: biom not found in PATH."; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python not found in PATH."; exit 1; }
command -v realpath >/dev/null 2>&1 || { echo "Error: realpath not found in PATH."; exit 1; }

[[ -f "$FASTQ_FILE" ]] || { echo "Error: input FASTQ not found: $FASTQ_FILE"; exit 1; }
[[ -f "$SHARED_ASV_IDS" ]] || { echo "Error: shared ASV whitelist not found: $SHARED_ASV_IDS"; exit 1; }
[[ -f "$CLASSIFIER_FILE_READS" ]] || { echo "Error: missing $CLASSIFIER_FILE_READS"; exit 1; }
[[ -f "$CLASSIFIER_FILE_TAX" ]] || { echo "Error: missing $CLASSIFIER_FILE_TAX"; exit 1; }
[[ -f "make_donut_graph.py" ]] || { echo "Error: make_donut_graph.py not found."; exit 1; }

mkdir -p "$WORKDIR"
cp "$FASTQ_FILE" "$WORKDIR/"
cp "$SHARED_ASV_IDS" "$WORKDIR/"
cp "$CLASSIFIER_FILE_READS" "$WORKDIR/"
cp "$CLASSIFIER_FILE_TAX" "$WORKDIR/"
cp make_donut_graph.py "$WORKDIR/"
cd "$WORKDIR"

printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"
printf "%s\t%s\n" "$SAMPLE" "$(realpath "$FASTQ_FILE")" >> "$MANIFEST"

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime quality-filter q-score \
  --i-demux demux.qza \
  --p-min-quality 30 \
  --p-quality-window 3 \
  --p-min-length-fraction 0.75 \
  --p-max-ambiguous 0 \
  --o-filtered-sequences demux-q30.qza \
  --o-filter-stats q30-filter-stats.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux-q30.qza \
  --p-trunc-len "$TRUNC_LEN" \
  --p-trunc-q 30 \
  --p-max-ee 1.0 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats.qza \
  --o-base-transition-stats base-transition-stats.qza

qiime feature-table filter-features \
  --i-table table-dada2.qza \
  --m-metadata-file "$SHARED_ASV_IDS" \
  --o-filtered-table table-shared.qza

qiime feature-table filter-seqs \
  --i-data rep-seqs-dada2.qza \
  --i-table table-shared.qza \
  --o-filtered-data rep-seqs-shared.qza

qiime feature-table filter-features \
  --i-table table-shared.qza \
  --p-min-frequency "$MIN_FREQUENCY" \
  --o-filtered-table table-filtered.qza

qiime feature-table filter-seqs \
  --i-data rep-seqs-shared.qza \
  --i-table table-filtered.qza \
  --o-filtered-data rep-seqs-filtered.qza

qiime feature-classifier classify-consensus-blast \
  --i-query rep-seqs-filtered.qza \
  --i-reference-reads "$CLASSIFIER_FILE_READS" \
  --i-reference-taxonomy "$CLASSIFIER_FILE_TAX" \
  --p-perc-identity 0.97 \
  --p-query-cov 0.8 \
  --o-classification taxonomy.qza \
  --o-search-results blast-results.qza

qiime taxa filter-table \
  --i-table table-filtered.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude Unassigned \
  --p-include g__ \
  --o-filtered-table table-genus-assigned.qza

qiime taxa collapse \
  --i-table table-genus-assigned.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus-table

GENUS_TSV="${WORKDIR}_genus-table.tsv"

biom convert \
  -i exported-genus-table/feature-table.biom \
  -o "$GENUS_TSV" \
  --to-tsv

python make_donut_graph.py \
  --manifest "$MANIFEST" \
  --table "$GENUS_TSV" \
  --threshold "$THRESHOLD"

mv *_donut.png "${WORKDIR}_donut.png"
mv *_donut.svg "${WORKDIR}_donut.svg"

cp "$GENUS_TSV" "$ORIG_DIR/"
cp "${WORKDIR}_donut.png" "$ORIG_DIR/"
cp "${WORKDIR}_donut.svg" "$ORIG_DIR/"
