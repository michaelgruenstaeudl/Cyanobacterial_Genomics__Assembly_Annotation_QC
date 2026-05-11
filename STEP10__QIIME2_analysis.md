### Read processing of 16S rRNA amplicon sequences

#### Installation of the necessary tools
```bash
# Installation of seqkit (for correcting paired-end reads)
conda install -c bioconda seqkit

# Installation of PEAR (for pairing paired-end reads)
conda install -c bioconda pear
```

#### STEP 1. Filtering and pairing the paired-end reads
```bash

SAMPLE=30_1326789214_HTF

# Prepare the paired-end sequence reads
seqkit sana ${SAMPLE}_R1_001.fastq.gz -o ${SAMPLE}_R1_001.cleaned.fastq 2>${SAMPLE}_R1_001.cleaned.fastq.log
gzip ${SAMPLE}_R1_001.cleaned.fastq
seqkit sana ${SAMPLE}_R2_001.fastq.gz -o ${SAMPLE}_R2_001.cleaned.fastq 2>${SAMPLE}_R2_001.cleaned.fastq.log
gzip ${SAMPLE}_R2_001.cleaned.fastq
seqkit pair -1 ${SAMPLE}_R1_001.cleaned.fastq.gz -2 ${SAMPLE}_R2_001.cleaned.fastq.gz 2>${SAMPLE}_R2_001.cleaned.paired.log

# Paired-end sequence reads
pear -f ${SAMPLE}_R1_001.cleaned.paired.fastq.gz -r ${SAMPLE}_R2_001.cleaned.paired.fastq.gz -o 16S_rRNA_seq_${SAMPLE}
for i in 16S*_${SAMPLE}*.fastq; do gzip $i; done
```

---

### Comprehensive metagenomic analysis of a 16S rRNA amplicon sequencing sample using QIIME2

```bash
# Load QIIME2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2026.1

# Define log-file
LOG_FILE="16S_rRNA_seq_30_1326789214_HTF_METRICS.log"
```

#### STEP 2. Conduct the QIIME2 analysis including sequence classification

```bash
bash BASHSCRIPT_QIIME2_analysis.sh
```

---

#### STEP 3. Inference of key metrics

##### Number of raw reads per sample
```bash
# Summarize the number of raw reads
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
```

##### Number of ASVs (features)
```bash
# Log ASV inference method because DADA2 uses exact sequence variants rather than OTU clustering by similarity threshold
echo "ASV method: DADA2 (no similarity threshold, exact sequence inference)" >> $LOG_FILE

# Summarize feature table
qiime feature-table summarize \
  --i-table table.qza \
  --o-summary table-summary.qzv

# Extract total number of ASVs from the exported BIOM table
qiime tools export \
  --input-path table.qza \
  --output-path exported-table
ASV_COUNT=$(biom summarize-table -i exported-table/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of ASVs: ${ASV_COUNT}" >> $LOG_FILE
```

##### Number of genera
```bash
# Collapse ASVs to genus level and 
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

qiime feature-table summarize \
  --i-table genus-table.qza \
  --o-summary genus-table-summary.qzv

# Extract total number of genera from the exported collapsed BIOM table
qiime tools export \
  --input-path genus-table.qza \
  --output-path exported-genus
GENUS_COUNT=$(biom summarize-table -i exported-genus/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of genera: ${GENUS_COUNT}" >> $LOG_FILE
```

##### Number of families
```bash
# Collapse ASVs to family level
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table family-table.qza
qiime feature-table summarize \
  --i-table family-table.qza \
  --o-summary family-table-summary.qzv

# Extract total number of families from the exported collapsed BIOM table
qiime tools export \
  --input-path family-table.qza \
  --output-path exported-family
FAMILY_COUNT=$(biom summarize-table -i exported-family/feature-table.biom | awk -F': ' '/^Num observations:|^Number of observations:/ {print $2; exit}' | tr -d '[:space:]')
echo "Number of families: ${FAMILY_COUNT}" >> $LOG_FILE
```

##### Summary alpha diversity metrics (with values)
```bash
# Compute alpha diversity metrics
qiime diversity alpha \
  --i-table table.qza \
  --p-metric shannon \
  --o-alpha-diversity shannon-vector.qza
qiime diversity alpha \
  --i-table table.qza \
  --p-metric pielou_e \
  --o-alpha-diversity pielou-e-vector.qza
qiime diversity alpha \
  --i-table table.qza \
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
```

---

#### STEP 4. Visualization across samples

```bash
# Visualizing diversity at the taxonomic level of class
python stacked_taxonomy_bars.py . -r c -o HorseThief_metagenomics_rank-class -t 1 --metadata metadata_cyanobacteria_samples.tsv

# Visualizing diversity at the taxonomic level of family but only for cyanobacteria
python stacked_taxonomy_bars.py . -r f -o HorseThief_metagenomics_rank-family_only_cyanos -t 1 --filter-rank c --filter-name Cyanobacteriia --metadata metadata_cyanobacteria_samples.tsv
```
