#!/bin/bash
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --job-name=Kraken2_16GB
#SBATCH --mail-type=all
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00

# Load software
source /homes/mgruenstaeudl/.bashrc
conda activate kraken2_env

# Run Bactopia
INPATH=/homes/mgruenstaeudl/HorseThiefReservoir_Fall2024

## Illumina data
kraken2 \
  --db ${INPATH}/01b_Kraken2_filtering/k2_standard_16_GB_20260626 \
  --paired \
  --threads 16 \
  --confidence 0.06 \
  --report illumina.kraken.report \
  --output illumina.kraken.out \
  ${INPATH}/01a_processed_reads/Illumina_filt_R1_paired.fastq.gz \
  ${INPATH}/01a_processed_reads/Illumina_filt_R2_paired.fastq.gz

extract_kraken_reads.py \
  -k illumina.kraken.out \
  -r illumina.kraken.report \
  -s1 ${INPATH}/01a_processed_reads/Illumina_filt_R1_paired.fastq.gz \
  -s2 ${INPATH}/01a_processed_reads/Illumina_filt_R2_paired.fastq.gz \
  -o illumina_R1_cyanosOnly.fastq \
  -o2 illumina_R2_cyanosOnly.fastq \
  --taxid 0 1117 \
  --include-children \
  --fastq-output

gzip illumina_R1_cyanosOnly.fastq
gzip illumina_R2_cyanosOnly.fastq

## Nanopore data
kraken2 \
  --db ${INPATH}/01b_Kraken2_filtering/k2_standard_16_GB_20260626 \
  --threads 16 \
  --confidence 0.02 \
  --report nanopore.kraken.report \
  --output nanopore.kraken.out \
  ${INPATH}/01a_processed_reads/Nanopore_filtered.q75.fastq.gz

extract_kraken_reads.py \
  -k nanopore.kraken.out \
  -r nanopore.kraken.report \
  -s ${INPATH}/01a_processed_reads/Nanopore_filtered.q75.fastq.gz \
  -o nanopore_cyanosOnly.fastq \
  --taxid 0 1117 \
  --include-children \
  --fastq-output

gzip nanopore_cyanosOnly.fastq
