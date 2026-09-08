### Genome Assembly (Hybrid Assembly)

#### Installation of the necessary tools
```
# Installation of bactopia
mamba create -n bactopia -c conda-forge -c bioconda bactopia
```

#### Hybrid DeNovo Assembly using Bactopia

##### Assembly primarily through Illumina reads
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Unicycler to assemble the short reads first, then bridging gaps with long reads.
# Good for low-coverage noisy long-reads.

bactopia \
   --sample HorseThief_Unicycler \
   --r1 Illumina_filt_R1_paired.fastq.gz \
   --r2 Illumina_filt_R2_paired.fastq.gz \
   --ont Nanopore_filtered.q75.fastq.gz \
   --hybrid
   --phix /homes/mgruenstaeudl/references/phix174.fasta   # Removal of any PhiX spike-ins
```

##### Assembly primarily through Oxford Nanopore reads
```
conda activate bactopia

# Under this option, Bactopia creates a hybrid assembly using Dragonflye to assemble the long-reads first, then polishing the assembly with the short-reads.
# Good for high-coverage high-quality long-reads.

bactopia \
  --sample Limnothrix_Unicycler \
  --r1 ${INPATH}/Illumina_filt_R1_paired.fastq.gz \
  --r2 ${INPATH}/Illumina_filt_R2_paired.fastq.gz \
  --ont ${INPATH}/Nanopore_filtered.q75.fastq.gz \
  --hybrid \
  --genome_size 4600000 \
  --coverage 0 \
  --use_bbmap \
  --phix /homes/mgruenstaeudl/references/phix174.fasta \
  --phix_k 31 \
  --maq 20 \
  --minlength 75 \
  --trimq 15 \
  --min_contig_len 500 \
  --min_contig_cov 2 \
  --unicycler_mode conservative \
  --nfconfig /homes/mgruenstaeudl/hybrid_assembly/Limnothrix_denovo_assembly.conf
```