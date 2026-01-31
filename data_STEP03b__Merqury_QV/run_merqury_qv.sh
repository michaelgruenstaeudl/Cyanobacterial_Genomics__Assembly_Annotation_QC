#!/bin/bash
#SBATCH --job-name=merqury_qv
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz

ASM1=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Merqury/input/BacterialGenome_Bactopia.fasta
ASM2=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Merqury/input/PlasmidGenome_Bactopia.fasta

OUT=merqury_qv
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

# Keep Merqury-related logs inside $OUT (even if something calls ./logs/)
mkdir -p logs

PREFIX="$(basename "$OUT")_$(date +%Y-%m-%d)"

# STEP 1. Build read k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "${PREFIX}.reads.meryl"

# STEP 2. Compute QV for two assemblies (qv.sh writes ${PREFIX}.qv)
QV_SH="$CONDA_PREFIX/share/merqury/eval/qv.sh"

"$QV_SH" \
  "${PREFIX}.reads.meryl" \
  "$ASM1" \
  "$ASM2" \
  "${PREFIX}"
