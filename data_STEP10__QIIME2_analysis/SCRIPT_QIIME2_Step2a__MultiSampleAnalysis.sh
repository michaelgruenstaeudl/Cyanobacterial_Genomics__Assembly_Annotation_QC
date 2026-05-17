#!/usr/bin/env bash
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1

TRUNC_LEN="250"
MIN_FREQUENCY="100"
MIN_SAMPLES="3"

MANIFEST="manifest_all.tsv"
WORKDIR="16S_rRNA_seq_multisample_filter"

mkdir -p "$WORKDIR"

printf "sample-id\tabsolute-filepath\n" > "$MANIFEST"

for FASTQ in 16S_rRNA_seq_*.assembled.fastq.gz; do
  SAMPLE_ID="${FASTQ#16S_rRNA_seq_}"
  SAMPLE_ID="${SAMPLE_ID%.assembled.fastq.gz}"
  printf "%s\t%s\n" "$SAMPLE_ID" "$(realpath "$FASTQ")" >> "$MANIFEST"
done

cp "$MANIFEST" "$WORKDIR/"
cd "$WORKDIR"

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
  --p-min-frequency "$MIN_FREQUENCY" \
  --p-min-samples "$MIN_SAMPLES" \
  --o-filtered-table table-filtered.qza

qiime feature-table filter-seqs \
  --i-data rep-seqs-dada2.qza \
  --i-table table-filtered.qza \
  --o-filtered-data rep-seqs-filtered.qza

qiime tools export \
  --input-path rep-seqs-filtered.qza \
  --output-path exported-shared-asvs

#cp exported-shared-asvs/dna-sequences.fasta ../shared_ASVs.fasta

awk '/^>/ {sub(/^>/, "", $1); print $1}' exported-shared-asvs/dna-sequences.fasta > shared_ASV_ids.tmp
printf "feature-id\n" > shared_ASV_ids.tsv
cat shared_ASV_ids.tmp >> shared_ASV_ids.tsv
rm shared_ASV_ids.tmp

cp shared_ASV_ids.tsv ../shared_ASV_ids.tsv
