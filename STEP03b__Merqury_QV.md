### Inference of QV for genome assembly via Mercury

#### Installation of Mercury on Beocat
```bash
# Load compiler
module load GCC
# Installation of Mercury
conda config --set solver classic
conda config --set channel_priority strict
env -u LD_LIBRARY_PATH conda create -n merqury -c conda-forge -c bioconda merqury --solver=classic
# Activation of Mercury
conda activate merqury
```

#### Running Merqury

##### Set up an SBATCH file for evaluation of a single assembly: `run_merqury_qv.sh`
```bash

#!/bin/bash
#SBATCH --job-name=merqury_qv
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz
ASM=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/BacterialGenome_Bactopia.fasta

OUT=merqury_qv
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
mkdir -p "logs/$OUT"
PREFIX="${OUT}_$(date +%Y-%m-%d)"

# 1) Build read k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "$OUT/${PREFIX}.reads.meryl"

# 2) Build assembly k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$ASM" \
  output "$OUT/${PREFIX}.asm.meryl"

# 3) Compute assembly-only k-mers (those in asm but not in reads)
meryl difference \
  "$OUT/${PREFIX}.asm.meryl" \
  "$OUT/${PREFIX}.reads.meryl" \
  output "$OUT/${PREFIX}.asm.only.meryl"

# 4) Compute QV
qv.sh \
  "$OUT/${PREFIX}.asm.meryl" \
  "$OUT/${PREFIX}.asm.only.meryl" \
  > "$OUT/${PREFIX}_QV.txt"

echo "Wrote:"
echo "  $OUT/${PREFIX}_QV.txt"
```
