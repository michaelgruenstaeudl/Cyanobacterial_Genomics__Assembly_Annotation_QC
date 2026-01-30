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

##### Set up an SBATCH file for evaluation of assembly quality: `run_merqury_qv.sh`
```bash
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
ASM=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Merqury/input/BacterialGenome_Bactopia.fasta

OUT=merqury_qv
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

# Keep Merqury-related logs inside $OUT (even if something calls ./logs/)
mkdir -p logs

PREFIX="$(basename "$OUT")_$(date +%Y-%m-%d)"

# 1) Build read k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "${PREFIX}.reads.meryl"

# 2) Build assembly k-mer database
meryl count k="$K" threads="$THREADS" memory=8 \
  "$ASM" \
  output "${PREFIX}.asm.meryl"

# 3) Compute assembly-only k-mers (those in asm but not in reads)
meryl difference \
  "${PREFIX}.asm.meryl" \
  "${PREFIX}.reads.meryl" \
  output "${PREFIX}.asm.only.meryl"

# 4) Compute QV
qv.sh \
  "${PREFIX}.asm.meryl" \
  "${PREFIX}.asm.only.meryl" \
  > "${PREFIX}_QV.txt"

echo "Wrote:"
echo "  $OUT/${PREFIX}_QV.txt"
```

#--- FILE HYGIENE ------------------------------------------------------

# This script runs from OUTSIDE $OUT (do not cd).
OUT=merqury_qv

shopt -s nullglob

# Compress all meryl databases created by the QV workflow:
#   *.reads.meryl  *.asm.meryl  *.asm.only.meryl
for mdb in "$OUT"/*.meryl; do
  if [[ -d "$mdb" ]]; then
    tarball="${mdb}.tar.gz"

    echo "  -> Compressing $(basename "$mdb")"
    tar -czf "$tarball" -C "$OUT" "$(basename "$mdb")"

    # Verify archive exists and is non-empty
    if [[ -s "$tarball" ]]; then
      rm -rf "$mdb"
      echo "     Removed original $(basename "$mdb")"
    else
      echo "     ERROR: Failed to create $(basename "$tarball"); keeping original" >&2
    fi
  fi
done

shopt -u nullglob
