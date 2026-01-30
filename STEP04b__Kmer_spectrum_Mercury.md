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
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10

set -euo pipefail

eval "$(/homes/mgruenstaeudl/miniconda3/bin/conda shell.bash hook)"
conda activate merqury

#--- INPUT -------------------------------------------------------------
R1=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
R2=/homes/mgruenstaeudl/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz
ASM=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/FinalAssembly_Bactopia.fasta
OUT=merqury_single_assembly
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
mkdir -p "logs/$OUT"  # IMPORTANT for Merqury!
mkdir -p figs/
PREFIX="${OUT}_$(date +%Y-%m-%d)"

# Build the meryl database from paired-end reads
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "$OUT/${PREFIX}.meryl"

# Run Merqury (single assembly)
merqury.sh "$OUT/${PREFIX}.meryl" "$ASM" "$OUT/${PREFIX}"

#--- VISUALIZE ---------------------------------------------------------
Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f merqury_single_assembly/merqury_single_assembly_2026-01-29.FinalAssembly_Bactopia.spectra-cn.hist \
  -o figs/merqury_single_assembly_2026-01-29.FinalAssembly_Bactopia.spectra-cn \
  -x 10 -y 6 -m 200 -t line -p

```

##### Set up an SBATCH file for evaluation of two assemblies: `run_mercury_compare_two_asm.sh`
```bash
#!/bin/bash
#SBATCH --job-name=merqury_compare_two
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

ASM1=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/FinalAssembly_Bactopia.fasta
ASM2=/homes/mgruenstaeudl/data/Limnothrix/05b_kmer_spectrum_Mercury/input/PlasmidGenome_Bactopia.fasta

OUT=merqury_compare
K=21
THREADS="${SLURM_CPUS_PER_TASK:-10}"

#--- RUN ---------------------------------------------------------------
mkdir -p "$OUT"
mkdir -p "logs/$OUT"  # IMPORTANT for Merqury!
mkdir -p figs/
PREFIX="${OUT}_$(date +%Y-%m-%d)"

# Build the meryl database from paired-end reads
meryl count k="$K" threads="$THREADS" memory=8 \
  "$R1" "$R2" \
  output "$OUT/${PREFIX}.meryl"

# Run Merqury comparison (two assemblies)
merqury.sh "$OUT/${PREFIX}.meryl" "$ASM1" "$ASM2" "$OUT/${PREFIX}"

#--- VISUALIZE ---------------------------------------------------------
Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f merqury_compare/merqury_compare_2026-01-29.spectra-cn.hist \
  -o figs/merqury_compare_2026-01-29.spectra-cn \
  -x 10 -y 6 -m 200 -t line -p

Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f merqury_compare/merqury_compare_2026-01-29.FinalAssembly_Bactopia.spectra-cn.hist \
  -o figs/merqury_compare_2026-01-29.FinalAssembly_Bactopia.spectra-cn \
  -x 10 -y 6 -m 200 -t line -p

Rscript "$MERQURY/plot/plot_spectra_cn.R" \
  -f merqury_compare/merqury_compare_2026-01-29.Plasmid_Bactopia.spectra-cn.hist \
  -o figs/merqury_compare_2026-01-29.Plasmid_Bactopia.spectra-cn \
  -x 10 -y 6 -m 200 -n 100000 -t line -p
```
