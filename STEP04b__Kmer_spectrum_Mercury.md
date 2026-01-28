### K-mer spectrum analysis for bacterial genome assembly

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

##### Set up an SBATCH file for evaluation of a single assembly: `run_mercury_single_asm.sh`
```bash
#!/bin/bash
#SBATCH --job-name=merqury_single
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R2_paired.fastq.gz
ASM=/homes/mgruenstaeudl/data/Limnothrix/04_backmapping/01a_Illumina_input/FinalAssembly_Bactopia.fasta
OUT=merqury_single_assembly
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

# Stable sample prefix from R1
base="${R1_BN%.fastq.gz}"
base="${base%.fq.gz}"
base="${base%_R1*}"
base="${base%_1*}"

# Build the meryl database from paired-end reads (correct syntax for your meryl)
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1_BN" "$R2_BN" \
  output "${base}.meryl"

# Run Merqury (single assembly)
merqury.sh "${base}.meryl" "$ASM_BN" > "${OUT}_merqury.out" 2> "${OUT}_merqury.err"

```

##### Set up an SBATCH file for evaluation of two assemblies: `run_mercury_compare_two_asm.sh`
```bash
#!/bin/bash
#SBATCH --job-name=merqury_compare_two
#SBATCH --mail-user=m_gruenstaeudl@fhsu.edu
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/denovo_assembly/Limnothrix/Illumina_filt_R2_paired.fastq.gz

ASM1=/homes/mgruenstaeudl/data/Limnothrix/04_backmapping/01a_Illumina_input/FinalAssembly_Bactopia.fasta
ASM2=/homes/mgruenstaeudl/data/Limnothrix/04_backmapping/01a_Illumina_input/Plasmid_Bactopia.fasta

OUT=merqury_compare
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
cd "$OUT"

ln -sf "$R1" .
ln -sf "$R2" .
ln -sf "$ASM1" .
ln -sf "$ASM2" .

R1_BN="$(basename "$R1")"
R2_BN="$(basename "$R2")"
ASM1_BN="$(basename "$ASM1")"
ASM2_BN="$(basename "$ASM2")"

# Stable sample prefix from R1
base="${R1_BN%.fastq.gz}"
base="${base%.fq.gz}"
base="${base%_R1*}"
base="${base%_1*}"

# Build the meryl database from paired-end reads (syntax compatible with your meryl)
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1_BN" "$R2_BN" \
  output "${base}.meryl"

# Run Merqury comparison (two assemblies)
merqury.sh "${base}.meryl" "$ASM1_BN" "$ASM2_BN" > "${OUT}_merqury.out" 2> "${OUT}_merqury.err"

```
