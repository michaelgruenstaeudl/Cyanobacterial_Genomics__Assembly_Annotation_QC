### Converting assemblies to complete bacterial genomes

#### Using Circlator to evaluate overlap of linear contig
```
# Installation
pip install circlator
sudo apt install mummer bwa samtools spades prodigal
conda install bioconda::canu
circlator progcheck
 
# Run Circlator
circlator fixstart Limnothrix_polished_genome.fa fixstart_out

# Check the output
head fixstart_out.fasta

# Confirm overlap of contig ends (500 bp)
tail -c 1000 fixstart_out.fasta | head -c 500 > end.fa && head -c 500 fixstart_out.fasta > start.fa && diff start.fa end.fa && echo "Checked 500 bp overlap at contig ends"


# Run Circlator on mitochondrial genome
circlator all --genes_fa Dloop.fasta MG936619_too_long_and_wrong_start.fasta SRR6664769_1.fastq too_long_and_wrong_start  
## Doesn't work with this dataset because number of reads too small

circlator fixstart --genes_fa Dloop.fasta MG936619_too_long_and_wrong_start.fasta fixing_start_of_too_long  ## Doesn't work 
```

#### 1. Using rotate to rotate mitogenome to correct start
https://github.com/richarddurbin/rotate
See examples in: https://wellcomeopenresearch.org/articles/8-401

# rotate to anchor string, allowing for 4 mismatches
./rotate -s TACGACCTCGATGTTGGATCA -m 4 mammalia.fa > mammalia.rotated.fa


#### 2a. If rotated mitogenome too long, trim off extra end
Can be done via custom script:
1. self-align a mitochondrial FASTA using MUMmer
2. detect a high-identity overlap between the beginning and end of the sequence
3. trim the duplicated region
4. output a circularized FASTA.

```
#!/usr/bin/env bash

# Usage:
# trim_circular_overlap.sh mito.fasta trimmed.fasta

set -e

FASTA=$1
OUT=$2

PREFIX=mito_self

# Self alignment
nucmer --maxmatch -p $PREFIX $FASTA $FASTA > /dev/null

# Get coordinates
show-coords -rcl $PREFIX.delta > ${PREFIX}.coords

# Find terminal overlap
OVERLAP=$(awk '
NR>5 {
    if ($1 < 1000 && $4 > $7-1000 && $7-$4 > 200) {
        print $7-$4
        exit
    }
}' ${PREFIX}.coords)

if [ -z "$OVERLAP" ]; then
    echo "No terminal overlap detected."
    cp $FASTA $OUT
    exit 0
fi

echo "Detected overlap: $OVERLAP bp"

LEN=$(grep -v ">" $FASTA | tr -d '\n' | wc -c)
TRIM=$((LEN - OVERLAP))

echo "Original length: $LEN"
echo "Trimmed length:  $TRIM"

# Extract trimmed sequence
seq=$(grep -v ">" $FASTA | tr -d '\n')

echo ">circularized_mito" > $OUT
echo ${seq:0:$TRIM} >> $OUT
```

#### 2b. If rotated mitogenome too short, extend via Novoplasty
Use "partial_mito.fasta" as "seed input":
```
Project:
-----------------------
Project name          = mito_extend
Type                  = mito
Genome Range          = 14000-20000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = yes
Seed Input            = partial_mito.fasta
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = reads_R1.fastq.gz
Reverse reads         = reads_R2.fastq.gz

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = yes
```

#### If 2b. successful, do 2a.

