### Read Filtering

#### Installation of tools for quality filtering
```
# Installation of filtlong
cd ~
mkdir -p git
mkdir -p bin
cd git
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j
cp bin/filtlong ~/bin/filtlong
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

#### Quality filtering of sequence reads

##### Illumina reads
```
# Combining Illumina reads
LOC_ILLUMINA_DAT1_R1=~/data/HorseThief/HorseThief_Illumina_Jan2025_run1_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT1_R2=~/data/HorseThief/HorseThief_Illumina_Jan2025_run1_R2.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R1=~/data/HorseThief/HorseThief_Illumina_Jan2025_run2_R1.fastq.gz  # Change!
LOC_ILLUMINA_DAT2_R2=~/data/HorseThief/HorseThief_Illumina_Jan2025_run2_R2.fastq.gz  # Change!

cat $LOC_ILLUMINA_DAT1_R1 $LOC_ILLUMINA_DAT2_R1 > Illumina_all_R1.fastq.gz
cat $LOC_ILLUMINA_DAT1_R2 $LOC_ILLUMINA_DAT2_R2 > Illumina_all_R2.fastq.gz

# Filtering the Illumina reads using \texttt{Trimmomatic
module load Trimmomatic

# Note: This will take 20-30 min!
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 \
  Illumina_all_R1.fastq.gz Illumina_all_R2.fastq.gz \
  Illumina_filt_R1_paired.fastq.gz Illumina_filt_R1_unpaired.fastq.gz \
  Illumina_filt_R2_paired.fastq.gz Illumina_filt_R2_unpaired.fastq.gz \
  ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  
rm Illumina_all_R1.fastq.gz Illumina_all_R2.fastq.gz
```

##### NanoPore reads
```
# Combining Nanopore reads
LOC_NANOPORE_DAT1=~/data/HorseThief/Nanopore_2025_06_12__fastq_pass/  # Change!
LOC_NANOPORE_DAT2=~/data/HorseThief/Nanopore_2025_06_24__fastq_pass/  # Change!

zcat $LOC_NANOPORE_DAT1/*.fastq.gz >  Nanopore_all.fastq
zcat $LOC_NANOPORE_DAT2/*.fastq.gz >> Nanopore_all.fastq

# Filtering the Nanopore reads using filtlong
filtlong --min_length 1 --keep_percent 75 Nanopore_all.fastq | 
	gzip > Nanopore_filtered.q75.fastq.gz

rm Nanopore_all.fastq
```

---

#### Installation of tools for taxonomic filtering
```
# Installation of Kraken2 and KrakenTools
conda create -n kraken2_env \
  -c conda-forge \
  -c bioconda \
  kraken2 krakentools
conda activate kraken2_env

# Download database (80 GB large!!!) and move into folder
mkdir -p k2_standard_20260626
wget  https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz
tar -xzf k2_standard_20260626.tar.gz \
  -C k2_standard_20260626
```

#### Taxonomic filtering of sequence reads

##### Illumina reads
```
kraken2 \
  --db k2_standard_20260626 \
  --paired \
  --threads 8 \
  --confidence 0.10 \
  --report illumina.kraken.report \
  --output illumina.kraken.out \
  ../01a_processed_reads/Illumina_filt_R1_paired.fastq.gz \
  ../01a_processed_reads/Illumina_filt_R2_paired.fastq.gz

extract_kraken_reads.py \
  -k illumina.kraken.out \
  -r illumina.kraken.report \
  -s1 ../01a_processed_reads/Illumina_filt_R1_paired.fastq.gz \
  -s2 ../01a_processed_reads/Illumina_filt_R2_paired.fastq.gz \
  -o illumina_R1_cyanosOnly.fastq \
  -o2 illumina_R2_cyanosOnly.fastq \
  --taxid 1117 \
  --include-children \
  --fastq-output
```

##### Oxford Nanopore
```
kraken2 \
  --db k2_standard_20260626 \
  --threads 16 \
  --confidence 0.10 \
  --report nanopore.kraken.report \
  --output nanopore.kraken.out \
  ../01a_processed_reads/Nanopore_filtered.q75.fastq.gz

extract_kraken_reads.py \
  -k nanopore.kraken.out \
  -r nanopore.kraken.report \
  -s ../01a_processed_reads/Nanopore_filtered.q75.fastq.gz \
  -o nanopore_cyanosOnly.fastq \
  --taxid 1117 \
  --include-children \
  --fastq-output
```