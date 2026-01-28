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

##### Set up an SBATCH file for evaluation of a single assembly: `run_mercury_qv.sh`
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
R1=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R1_paired.fastq.gz
ASM=/homes/mgruenstaeudl/data/Limnothrix/04_backmapping/01a_Illumina_input/FinalAssembly_Bactopia.fasta

OUT=merqury_qv
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

ln -sf "$R1" .
ln -sf "$R2" .
ln -sf "$ASM" .

R1_BN="$(basename "$R1")"
R2_BN="$(basename "$R2")"
ASM_BN="$(basename "$ASM")"

# Stable prefix from R1
base="${R1_BN%.fastq.gz}"
base="${base%.fq.gz}"
base="${base%_R1*}"
base="${base%_1*}"

# 1) Build read k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1_BN" "$R2_BN" \
  output "${base}.meryl"

# 2) Build assembly k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$ASM_BN" \
  output "asm.meryl"

# 3) Compute assembly-only k-mers (those in asm but not in reads)
meryl difference "asm.meryl" "${base}.meryl" output "asm.only.meryl"

# 4) Compute QV
# Merqury's qv.sh prints QV and related stats; capture output to a file.
qv.sh "asm.meryl" "asm.only.meryl" > "${OUT}_QV.txt"

echo "Wrote:"
echo "  $(pwd)/${OUT}_QV.txt"
```
