# Cyanobacterial Genome Assembly

Scripts and workflows for the **assembly, polishing, quality assessment, and visualization** of complete cyanobacterial genomes.

---

## Overview

This repository contains step-by-step protocols and scripts for generating high-quality cyanobacterial genome assemblies from sequencing data. The workflow covers hybrid assembly, contig circularization, quality evaluation, and multiple visualization approaches commonly used in comparative and structural genomics.

---

## Workflow Structure

Each chapter corresponds to one logical step in the genome assembly and evaluation pipeline. Chapters are designed to be followed sequentially, but individual steps can also be used independently.

---

## Chapters

### 01. Genome Assembly
- [Hybrid genome assembly](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP01__Read_filtering_and_genome_assembly.md)

---

### 02. Contig Processing
- [Contig circularization](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP02__Circularization_of_contig.md)

---

### 03. Assembly Quality Assessment
- [Quality check via QUAST](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03a__Quality_eval_via_QUAST.md)

- [Inference of Merqury QV values](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03b__Merqury_QV.md)


###### What the .qv file format means

The Merqury QV output file has five columns:

```
<label>   <err_kmers>   <total_kmers>   <QV>   <error_rate>
```

Column descriptions

| Column      | Meaning                                          |
| ----------- | ------------------------------------------------ |
| label       | Assembly name (or combined assemblies)           |
| err_kmers   | Number of assembly k-mers not found in the reads |
| total_kmers | Total number of assembly k-mers                  |
| QV          | Phred-scaled quality value                       |
| error_rate  | Estimated base error rate                        |



###### RESULTS FOR BACTERIAL GENOME

```text
<label>                         <err_kmers>   <total_kmers>   <QV>     <error_rate>
BacterialGenome_Bactopia        0             4536393        +inf    0
PlasmidGenome_Bactopia          0             4018           +inf    0
both                            0             4540411        +inf    0
```

**Notes:**
* A `+inf` QV indicates zero observed error k-mers relative to the reads.
* An `error_rate` of `0` reflects no estimated base errors under this metric.


###### **TO DO:** Do both quality checks (QUAST, Mercury) also for plasmid genome

---

### 04. k-mer Spectrum Analysis
- [Visualization of k-mer spectra via Jellyfish](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP04a__Kmer_spectrum_Jellyfish.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP04a__Kmer_spectrum_Jellyfish/kmer_spectrum_Illumina_k21_ONT_k17_combined.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Extend the x-axis in this figure.

- [Visualization of k-mer spectra via Merqury](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP04b__Kmer_spectrum_Merqury.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP04b_Kmer_spectrum_Merqury/output/merqury_compare_2026-01-30.BothGenomes.spectra-cn.fl.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Create the same figure based on the Nanopore data; maybe create a facet plot for both.

---

### 05. Whole-Genome Alignment
- [Visualization of MUMmer4 results as dotplots](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP05a__MUMmer4_dotplots.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP05__MUMmer4_dotplots/Limnothrix_sp_BLA16_vs_BactopiaAssembly.dotplot.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

- [Inference of inversions](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP05b__Inversions.md)

```text
Limnothrix_sp_BL_A_16_CP166615	82994	113322	30329
Limnothrix_sp_BL_A_16_CP166615	592628	736045	143418
Limnothrix_sp_BL_A_16_CP166615	1756173	1817785	61613
Limnothrix_sp_BL_A_16_CP166615	2674925	2733309	58385
Limnothrix_sp_BL_A_16_CP166615	3085082	3131342	46261
Limnothrix_sp_BL_A_16_CP166615	4038055	4191007	152953
```

---

### 06. Coverage Visualization
- [Visualization of sequencing coverage via Circleator](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP06__Coverage_Viz_Circleator.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP06__Coverage_Viz_Circleator/02_output/FinalAssembly_Bactopia_Circleator_plusLegend.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Do coverage visualization also for plasmid genome

---

### 07. Synteny and Structural Analysis
- Show inversions within the assembly using Circos

- [Show synteny and collinearity between *Limnothrix* B-16 and the assembly using Circos](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP07b__Circos__synteny_collinearity_across_genomes.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP07b__Circos__synteny_collinearity_across_genomes/process_and_output/circos/circos.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">


---

### 08. Gene-Level Visualization
- [Visualization of gene location via GenoVi](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08a__Visualization_via_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08a__Visualization_via_GenoVi/FinalAssembly_Bactopia__output_from_GenoVi.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Do GenoVi visualization also for plasmid genome

- [Barchart-style gene tally compatible with GenoVi visualization](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08b__Gene_tally_of_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08b__Gene_tally_of_GenoVi/genovi_COG_Classification_COG_barplot.png" style="display:block; margin-left:auto; margin-right:auto; width:100%;">


### 09. Tabular summary of genome characteristics and annotation
- [Generate a tabular summary of the characteristics and annotations of an input genome](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP09__Generate_tabular_genome_summary.py)

---

## Notes

- Image files are stored in the corresponding `data_STEPXX__*` directories.
- All workflows assume Linux environments and standard bioinformatics toolchains.
- Individual steps can be adapted to other bacterial genomes with minimal modification.

---
