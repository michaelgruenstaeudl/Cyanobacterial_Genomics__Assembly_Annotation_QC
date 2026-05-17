
#### STEP 1. Read processing of 16S rRNA amplicon sequences
```bash
# Installation of seqkit (for correcting paired-end reads)
conda install -c bioconda seqkit
# Installation of PEAR (for pairing paired-end reads)
conda install -c bioconda pear

for SAMPLE in 30_1305444119_HT2 30_1326789214_HTF 30_1302373217_HT3; do

  # Filtering the paired-end reads
  seqkit sana ${SAMPLE}_R1_001.fastq.gz -o ${SAMPLE}_R1_001.cleaned.fastq 2>${SAMPLE}_R1_001.cleaned.fastq.log
  gzip ${SAMPLE}_R1_001.cleaned.fastq
  seqkit sana ${SAMPLE}_R2_001.fastq.gz -o ${SAMPLE}_R2_001.cleaned.fastq 2>${SAMPLE}_R2_001.cleaned.fastq.log
  gzip ${SAMPLE}_R2_001.cleaned.fastq
  seqkit pair -1 ${SAMPLE}_R1_001.cleaned.fastq.gz -2 ${SAMPLE}_R2_001.cleaned.fastq.gz 2>${SAMPLE}_R2_001.cleaned.paired.log

  # Pairing the paired-end reads
  pear -f ${SAMPLE}_R1_001.cleaned.paired.fastq.gz -r ${SAMPLE}_R2_001.cleaned.paired.fastq.gz -o 16S_rRNA_seq_${SAMPLE}
  for i in 16S*_${SAMPLE}*.fastq; do gzip $i; done

done
```

---

#### STEP 2. Comprehensive metagenomic analysis of a 16S rRNA amplicon sequencing sample using QIIME2

```bash
# Define input and output
LOCATION="HorseThief_Reservoir"

# Conduct the QIIME2 analysis on all samples
# Purpose: Generating ASVs across samples for estimating optimal quality filtering and minimal frequency levels
bash SCRIPT_QIIME2_Step2a__MultiSampleAnalysis.sh

# Conduct the QIIME2 analysis (incl. sequence classification) on one sample at a time
for i in 30_1305444119_HT2 30_1326789214_HTF 30_1302373217_HT3; do
  bash SCRIPT_QIIME2_Step2b__SingleSampleAnalysis.sh $i;
done
```

---

#### STEP 3. Inference of key metrics

```bash
# Inference of key metrics
for i in 30_1305444119_HT2 30_1326789214_HTF 30_1302373217_HT3; do
  bash SCRIPT_QIIME2_Step3__Infer_key_metrics.sh $i;
done
```
