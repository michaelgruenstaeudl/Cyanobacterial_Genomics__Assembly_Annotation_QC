### Converting assemblies to complete bacterial genomes

Note: The following code doesn't work properly at the moment!

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
