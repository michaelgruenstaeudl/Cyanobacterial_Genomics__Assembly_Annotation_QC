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
## Outer coordinate circle: reduce clutter
##  - ticks every 50 kb
##  - labels every 250 kb
coords outerf=1.07,tick-interval=50000,label-interval=250000,label-units=kb,label-precision=0

## Optional GC% track (smoothed)
#small-cgap
#%GC0-100 innerf=1.03,heightf=0.07,opacity=0.6,graph-min=0,graph-max=100,window-size=1000

### forward and reverse strand genes
small-cgap
genes-fwd heightf=0.05
tRNAs-fwd heightf=0.05,innerf=same
rRNAs-fwd heightf=0.05,innerf=same
tiny-cgap
genes-rev heightf=0.05
tRNAs-rev heightf=0.05,innerf=same
rRNAs-rev heightf=0.05,innerf=same

## Coverage track (smoothed; NOT nucleotide-by-nucleotide)
## Choose a window-size that matches what you want to see:
##  - 1000 = more detail
##  - 2000 = very smooth
##  - 5000 = easier for file size and display -- SELECTED
small-cgap
new cov graph 0.14 graph-function=BAMCoverage,bam-file=FinalAssembly_Bactopia.backmap.sorted.bam,bam-seqid=FinalAssembly_Bactopia,graph-min=0,graph-max=data_max,window-size=5000,heightf=0.25,opacity=0.85,color1=#3aa0d5
```

##### Running Circleator

```bash
MYSAMPLE=FinalAssembly_Bactopia
circleator --data=${MYSAMPLE}_corrected.gbk --config=${MYSAMPLE}.circleator.conf > ${MYSAMPLE}.svg
```