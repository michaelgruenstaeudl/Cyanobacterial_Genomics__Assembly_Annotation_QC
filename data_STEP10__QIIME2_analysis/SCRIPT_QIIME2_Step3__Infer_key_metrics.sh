#!/usr/bin/env bash
set -uo pipefail

SAMPLE="${1:?Usage: bash SCRIPT_QIIME2_Step3__Infer_key_metrics.sh SAMPLE_NAME}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1

ORIG_DIR="$(pwd)"
WORKDIR="16S_rRNA_seq_${SAMPLE}"

command -v qiime >/dev/null 2>&1 || { echo "Error: qiime not found in PATH."; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "Error: awk not found in PATH."; exit 1; }
command -v biom >/dev/null 2>&1 || { echo "Error: biom not found in PATH."; exit 1; }

LOG_FILE="${ORIG_DIR}/${WORKDIR}_stats.log"

# Initialize log file immediately so it always exists
: > "$LOG_FILE"
echo "Sample: ${SAMPLE}" >> "$LOG_FILE"

[[ -d "$WORKDIR" ]] || { echo "Error: work directory not found: $WORKDIR" | tee -a "$LOG_FILE"; exit 1; }
cd "$WORKDIR"

# Prefer the most downstream available table artifact from Step 3 output
if [[ -f "table-filtered.qza" ]]; then
  TABLE_ARTIFACT="table-filtered.qza"
elif [[ -f "table-dada2.qza" ]]; then
  TABLE_ARTIFACT="table-dada2.qza"
elif [[ -f "table.qza" ]]; then
  TABLE_ARTIFACT="table.qza"
else
  echo "Error: none of table-filtered.qza, table-dada2.qza, table.qza found in $WORKDIR" | tee -a "$LOG_FILE"
  exit 1
fi

TAXONOMY_ARTIFACT="taxonomy.qza"
[[ -f "$TAXONOMY_ARTIFACT" ]] || { echo "Error: taxonomy artifact not found: $TAXONOMY_ARTIFACT" | tee -a "$LOG_FILE"; exit 1; }
echo "Using table artifact: $TABLE_ARTIFACT" >> "$LOG_FILE"

##### Number of raw reads per sample
# Summarize the number of raw reads
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux-summary.qzv

qiime tools export \
  --input-path demux-summary.qzv \
  --output-path exported-demux-summary
RAW_FILE=$(find exported-demux-summary -type f | grep -i "sample" | grep -E '\.tsv$|\.csv$' | head -n 1)

# Extract all sample IDs and their read counts
awk '
BEGIN{FS="[,\t]"}
NR==1 {
  for(i=1;i<=NF;i++){
    if(tolower($i) ~ /sample/) sid=i
    if(tolower($i) ~ /count|sequence/) cnt=i
  }
}
NR>1 {
  printf "Raw reads (%s): %s\n", $sid, $cnt
}
' "$RAW_FILE" >> $LOG_FILE

##### Number of ASVs (features)
# Log ASV inference method because DADA2 uses exact sequence variants rather than OTU clustering by similarity threshold
echo "ASV method: DADA2 (no similarity threshold, exact sequence inference)" >> $LOG_FILE

# Summarize feature table
qiime feature-table summarize \
  --i-table "$TABLE_ARTIFACT" \
  --o-feature-frequencies table-feature-frequencies.qza \
  --o-sample-frequencies table-sample-frequencies.qza \
  --o-summary table-summary.qzv

# Extract total number of ASVs from the exported BIOM table
qiime tools export \
  --input-path "$TABLE_ARTIFACT" \
  --output-path exported-table
ASV_COUNT=$(biom summarize-table -i exported-table/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of ASVs: ${ASV_COUNT}" >> $LOG_FILE

##### Number of genera
# Collapse ASVs to genus level
qiime taxa collapse \
  --i-table "$TABLE_ARTIFACT" \
  --i-taxonomy "$TAXONOMY_ARTIFACT" \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

qiime feature-table summarize \
  --i-table genus-table.qza \
  --o-feature-frequencies genus-feature-frequencies.qza \
  --o-sample-frequencies genus-sample-frequencies.qza \
  --o-summary genus-table-summary.qzv

# Extract total number of genera from the exported collapsed BIOM table
qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus
GENUS_COUNT=$(biom summarize-table -i exported-genus/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of genera: ${GENUS_COUNT}" >> $LOG_FILE

##### Number of families
# Collapse ASVs to family level
qiime taxa collapse \
  --i-table "$TABLE_ARTIFACT" \
  --i-taxonomy "$TAXONOMY_ARTIFACT" \
  --p-level 5 \
  --o-collapsed-table family-table.qza

qiime feature-table summarize \
  --i-table family-table.qza \
  --o-feature-frequencies family-feature-frequencies.qza \
  --o-sample-frequencies family-sample-frequencies.qza \
  --o-summary family-table-summary.qzv

# Extract total number of families from the exported collapsed BIOM table
qiime tools export \
  --input-path family-table.qza \
  --output-path exported-family
FAMILY_COUNT=$(biom summarize-table -i exported-family/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of families: ${FAMILY_COUNT}" >> $LOG_FILE

##### Summary alpha diversity metrics (with values)
# Compute alpha diversity metrics
qiime diversity alpha \
  --i-table "$TABLE_ARTIFACT" \
  --p-metric shannon \
  --o-alpha-diversity shannon-vector.qza
qiime diversity alpha \
  --i-table "$TABLE_ARTIFACT" \
  --p-metric pielou_e \
  --o-alpha-diversity pielou-e-vector.qza
qiime diversity alpha \
  --i-table "$TABLE_ARTIFACT" \
  --p-metric simpson \
  --o-alpha-diversity simpson-vector.qza

# Export alpha diversity vectors
qiime tools export \
  --input-path shannon-vector.qza \
  --output-path exported-shannon
qiime tools export \
  --input-path pielou-e-vector.qza \
  --output-path exported-pielou-e
qiime tools export \
  --input-path simpson-vector.qza \
  --output-path exported-simpson

# Extract and log Shannon diversity values (rounded to 3 decimals)
SHANNON_FILE=$(find exported-shannon -type f | grep -E '\.tsv$' | head -n 1)
awk 'BEGIN{FS="\t"} NR>1 {printf "Shannon diversity (%s): %.3f\n", $1, $2}' "$SHANNON_FILE" >> $LOG_FILE

# Extract and log Pielou evenness values (rounded to 3 decimals)
PIELOU_FILE=$(find exported-pielou-e -type f | grep -E '\.tsv$' | head -n 1)
awk 'BEGIN{FS="\t"} NR>1 {printf "Pielou evenness (%s): %.3f\n", $1, $2}' "$PIELOU_FILE" >> $LOG_FILE

# Extract and log Simpson diversity values (rounded to 3 decimals)
SIMPSON_FILE=$(find exported-simpson -type f | grep -E '\.tsv$' | head -n 1)
awk 'BEGIN{FS="\t"} NR>1 {printf "Simpson diversity (%s): %.3f\n", $1, $2}' "$SIMPSON_FILE" >> $LOG_FILE
