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
- [Quality check via QUAST](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP03__Quality_eval_via_QUAST.md)
- Quality control analysis: inference of Mercury QV values

---

### 04. k-mer Spectrum Analysis
- [Visualization of k-mer spectra via Jellyfish](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP04a__Kmer_spectrum_Jellyfish.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP04a__Kmer_spectrum_Jellyfish/kmer_spectrum_Illumina_k21_ONT_k17_combined.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

---

### 05. Whole-Genome Alignment
- [Visualization of MUMmer4 results as dotplots](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP05__MUMmer4_dotplots.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP05__MUMmer4_dotplots/Limnothrix_sp_BLA16_vs_BactopiaAssembly.dotplot.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

---

### 06. Coverage Visualization
- [Visualization of sequencing coverage via Circleator](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP06__Coverage_Viz_Circleator.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP06__Coverage_Viz_Circleator/02_output/FinalAssembly_Bactopia_Circleator_plusLegend.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

---

### 07. Synteny and Structural Analysis
- Show inversions within the assembly using Circos

- Show synteny and collinearity between *Limnothrix* B-16 and the assembly using Circos


---

### 08. Gene-Level Visualization
- [Visualization of gene location via GenoVi](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08a__Visualization_via_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08a__Visualization_via_GenoVi/FinalAssembly_Bactopia__output_from_GenoVi.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

- [Barchart-style gene tally compatible with GenoVi visualization](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP08b__Gene_tally_of_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP08b__Gene_tally_of_GenoVi/genovi_COG_Classification_COG_barplot.png" style="display:block; margin-left:auto; margin-right:auto; width:100%;">

---

## Notes

- Image files are stored in the corresponding `data_STEPXX__*` directories.
- All workflows assume Linux environments and standard bioinformatics toolchains.
- Individual steps can be adapted to other bacterial genomes with minimal modification.

---
