### Visualization of sequencing coverage via Circleator

#### Backmapping of reads to complete genome assembly
Best done on a HPC cluster, as this step is computationally intensive.

##### Backmapping of Illumina reads
```bash

MYSAMPLE="FinalAssembly_Bactopia"

REF="${MYSAMPLE}.fasta"
R1="${MYSAMPLE}_R1_paired.fastq.gz"
R2="${MYSAMPLE}_R2_paired.fastq.gz"
THREADS=8

# Loading modules if conducted on HPC cluster
# module load BWA-MEM2
# If bwa-mem2 not installed as module, download:
# curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf -
module load SAMtools

# Checking that required tools available
command -v bwa-mem2 >/dev/null 2>&1 || { echo "ERROR: bwa-mem2 not found"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }

# Index reference
bwa-mem2 index "$REF"

# Map, sort, index
bwa-mem2 mem -t "$THREADS" \
  -R "@RG\tID:${MYSAMPLE}\tSM:${MYSAMPLE}\tPL:ILLUMINA" \
  "$REF" "$R1" "$R2" \
| samtools sort -@ "$THREADS" -o "${MYSAMPLE}.backmap.sorted.bam"

samtools index "${MYSAMPLE}.backmap.sorted.bam"

# Alignment summaries
samtools flagstat "${MYSAMPLE}.backmap.sorted.bam" > "${MYSAMPLE}.backmap.flagstat.txt"
samtools stats "${MYSAMPLE}.backmap.sorted.bam" > "${MYSAMPLE}.backmap.stats.txt"

# Per-base depth
samtools depth -a "${MYSAMPLE}.backmap.sorted.bam" > "${MYSAMPLE}.depth.txt"

# Genome-wide average coverage depth
awk '{sum += $3; count += 1} END {if (count > 0) print sum / count}' \
  "${MYSAMPLE}.depth.txt" > "${MYSAMPLE}.average_depth.txt"

```

##### Backmapping of ONT reads
```bash

MYSAMPLE="FinalAssembly_Bactopia"

REF="${MYSAMPLE}.fasta"
ONT_FASTQ="ont_reads.fastq.gz"
THREADS=8

# Define output files
BAM="${MYSAMPLE}.ont.backmap.sorted.bam"
DEPTH_TXT="${MYSAMPLE}.ont.depth.txt"
AVG_TXT="${MYSAMPLE}.ont.average_depth.txt"

# Loading modules if conducted on HPC cluster

# Checking that required tools available
command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not found"; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "ERROR: samtools not found"; exit 1; }

# Index the reference FASTA (useful for downstream tools)
samtools faidx "$REF"

# Map ONT reads to the reference and sort alignments
minimap2 -t "$THREADS" -ax map-ont "$REF" "$ONT_FASTQ" \
  | samtools sort -@ "$THREADS" -o "$BAM"

# Index the sorted BAM file
samtools index "$BAM"

# Basic alignment quality control statistics
samtools flagstat "$BAM" > "${MYSAMPLE}.ont.flagstat.txt"
samtools stats "$BAM"    > "${MYSAMPLE}.ont.stats.txt"

# Calculate per-base coverage depth
samtools depth -a "$BAM" > "$DEPTH_TXT"

# Calculate genome-wide average coverage depth
awk '{sum+=$3; n++} END {if(n>0) printf "%.6f\n", sum/n; else print "NA"}' "$DEPTH_TXT" > "$AVG_TXT"

```

#### Visualization of coverage across circular genome

##### Installation of Circleator and dependecies

###### On Debian
```bash
# Dependencies
sudo apt install perl bioperl libbatik-java vcftools samtools

# Additional perl modules
sudo cpan
install CPAN
reload cpan
install JSON
install Log::Log4perl
install SVG
install Text::CSV
install Bio::Perl
install Bio::FeatureIO::gff
install Module::Build
exit

# Installing Circleator systemwide
curl -L -o Circleator-1.0.2.tar.gz \
  https://github.com/jonathancrabtree/Circleator/archive/refs/tags/1.0.2.tar.gz

tar xzvf Circleator-1.0.2.tar.gz
cd Circleator-1.0.2

perl Build.PL
./Build
./Build test
sudo ./Build install
```

##### Configuration file for Circleator
Circleator draws figures based on the information in its configuration file

```
## =========================
## Circleator config: 4.5 Mb bacterial genome
## =========================
## "Genes": "#4daf4a", # medium green
## "tRNAs": "#984ea3", # purple
## "rRNAs": "#e41a1c", # red
## "Illumina Coverage": "#3aa0d5", # blue
    
coords outerf=1.07,tick-interval=50000,label-interval=250000,label-units=kb,label-precision=0

small-cgap

## Forward strand
genes-fwd  heightf=0.05,color1=#4daf4a
tRNAs-fwd  heightf=0.05,innerf=same,color1=#984ea3
rRNAs-fwd  heightf=0.05,innerf=same,color1=#e41a1c

tiny-cgap

## Reverse strand
genes-rev  heightf=0.05,color1=#4daf4a
tRNAs-rev  heightf=0.05,innerf=same,color1=#984ea3
rRNAs-rev  heightf=0.05,innerf=same,color1=#e41a1c

small-cgap
new cov graph 0.14 graph-function=BAMCoverage,bam-file=FinalAssembly_Bactopia.backmap.sorted.bam,bam-seqid=FinalAssembly_Bactopia,graph-min=0,graph-max=data_max,window-size=5000,heightf=0.25,opacity=0.85,color1=#3aa0d5
```

##### Running Circleator
Note: Make sure that the genome name in the GenBank flatfile (i.e., second column of lines LOCUS and ACCESSION) are identical to the genome name of the BAM file (i.e., find via: `samtools view -H FinalAssembly_Bactopia.backmap.sorted.bam | grep '^@SQ' | head
samtools idxstats FinalAssembly_Bactopia.backmap.sorted.bam`).

```bash
MYSAMPLE=FinalAssembly_Bactopia
circleator --data=${MYSAMPLE}_sameHeaderAsBAM.gbk --config=${MYSAMPLE}.circleator.conf > ${MYSAMPLE}.svg
```

##### Correcting the raw SVG
Circleator produces an SVG in the old SVG 1.0 format that is not rendered correctly in today's SVG editors. Hence, the output of Circleator must be rasterized first using [Apache Batik](https://xmlgraphics.apache.org/batik/download.html).

For that, download the [Batik binary](https://www.apache.org/dyn/closer.cgi?filename=/xmlgraphics/batik/binaries/batik-bin-1.19.tar.gz&action=download) into the same directopry as the output of Circleator and unzip it.

```bash
# Set BATIK_HOME and confirm you are in the right place
export BATIK_HOME="$PWD/batik-1.19"
ls "$BATIK_HOME"/lib | head

# Convert SVG to PDF using Batik with full classpath
java -cp "$BATIK_HOME/lib/*:$BATIK_HOME/extensions/*:$BATIK_HOME/batik-rasterizer-1.19.jar" \
  org.apache.batik.apps.rasterizer.Main \
  -m application/pdf \
  -d FinalAssembly_Bactopia.pdf \
  FinalAssembly_Bactopia.svg
```

##### Generate a color legend using Python
Circleator does not produce any legends for the figures it produces. Hence, the legends must be generated separately by the user.

```bash
python3 Circleator_legend_maker.sh
```