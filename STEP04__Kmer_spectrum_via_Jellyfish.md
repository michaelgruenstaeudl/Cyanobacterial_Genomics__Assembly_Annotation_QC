### K-mer spectrum analysis for bacterial genome assembly

#### Preparation
```bash
module load Jellyfish
```

#### K-mer counting with Jellyfish on Illumina data

```bash

# Parameters
K=21
THREADS=8
HASH="1G"

# Output directory
ILM_OUT=~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_Illumina_k${K}
mkdir -p "$ILM_OUT"
cd "$ILM_OUT"

# Input reads (edit as needed; FASTQ or FASTQ.GZ both work with zcat -f)
ILLUMINA_READS=(
  ~/data/Limnothrix/02_processed_reads/Illumina_filt_R1_paired.fastq.gz
  ~/data/Limnothrix/02_processed_reads/Illumina_filt_R2_paired.fastq.gz
)

# Count canonical k-mers (-C) and write Jellyfish database
# zcat -f works for both .fastq and .fastq.gz
jellyfish count -C -m "$K" -s "$HASH" -t "$THREADS" \
  /dev/fd/0 -o "illumina_reads_k${K}.jf" < <(zcat -f "${ILLUMINA_READS[@]}")

# Make k-mer multiplicity histogram
jellyfish histo -t "$THREADS" "illumina_reads_k${K}.jf" > "illumina_reads_k${K}.histo"
```

#### K-mer counting with Jellyfish on ONT data

```bash

# Parameters
K=17
THREADS=8
HASH="1G"

# Output directory
ONT_OUT=~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_ONT_k${K}
mkdir -p "$ONT_OUT"
cd "$ONT_OUT"

# Input reads (edit as needed; FASTQ or FASTQ.GZ)
ONT_READS=(
  ~/data/Limnothrix/02_processed_reads/Nanopore_filtered.q75.fastq.gz
)

# Count canonical k-mers (-C) and write Jellyfish database
jellyfish count -C -m "$K" -s "$HASH" -t "$THREADS" \
  /dev/fd/0 -o "ont_reads_k${K}.jf" < <(zcat -f "${ONT_READS[@]}")

# Make k-mer multiplicity histogram
jellyfish histo -t "$THREADS" "ont_reads_k${K}.jf" > "ont_reads_k${K}.histo"
```


#### Visualization of Jellyfish via Python's Seaborn

```python
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_histo(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None, names=["multiplicity", "n_kmers"])


def plot_spectrum(
    df: pd.DataFrame,
    title: str,
    out_svg: str,
    x_max: int = 200,
    log_y: bool = True,
):
    sns.set_theme(style="whitegrid")

    df_plot = df[df["multiplicity"] <= x_max].copy()

    plt.figure(figsize=(11, 6))
    ax = sns.barplot(
        data=df_plot,
        x="multiplicity",
        y="n_kmers",
    )

    ax.set_title(title)
    ax.set_xlabel("k-mer multiplicity")
    ax.set_ylabel("Number of distinct k-mers")

    if log_y:
        ax.set_yscale("log")

    step = 5 if x_max <= 100 else 10
    ax.set_xticks(list(range(0, x_max + 1, step)))

    plt.tight_layout()
    plt.savefig(out_svg, format="svg")
    plt.close()


# ---------------- Illumina ----------------
K_ILLUMINA = 21
ILM_BASE = os.path.expanduser(
    f"~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_Illumina_k{K_ILLUMINA}"
)

illumina_histo = os.path.join(
    ILM_BASE, f"illumina_reads_k{K_ILLUMINA}.histo"
)
illumina_svg = os.path.join(
    ILM_BASE, f"illumina_reads_k{K_ILLUMINA}_spectrum.svg"
)

df_illumina = load_histo(illumina_histo)
plot_spectrum(
    df_illumina,
    title=f"Illumina read k-mer spectrum (k={K_ILLUMINA})",
    out_svg=illumina_svg,
    x_max=200,
    log_y=True,
)


# ---------------- Oxford Nanopore ----------------
K_ONT = 17
ONT_BASE = os.path.expanduser(
    f"~/data/Limnothrix/05_kmer_spectrum/KmerSpectrum_ONT_k{K_ONT}"
)

ont_histo = os.path.join(
    ONT_BASE, f"ont_reads_k{K_ONT}.histo"
)
ont_svg = os.path.join(
    ONT_BASE, f"ont_reads_k{K_ONT}_spectrum.svg"
)

df_ont = load_histo(ont_histo)
plot_spectrum(
    df_ont,
    title=f"Oxford Nanopore read k-mer spectrum (k={K_ONT})",
    out_svg=ont_svg,
    x_max=200,
    log_y=True,
)

print("Wrote:")
print(" ", illumina_svg)
print(" ", ont_svg)

```