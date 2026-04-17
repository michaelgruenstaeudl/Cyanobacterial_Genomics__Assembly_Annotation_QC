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

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP04b__Kmer_spectrum_Merqury/output/merqury_compare_2026-01-30.BothGenomes.spectra-cn.fl.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

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

### 08. Evaluation of genome annotations
- [Evaluate if the reading frames of the genes of the genome are intact](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP08__Annotation_evaluation/PYSCRIPT_Evaluate_reading_frames_of_genes.py)


- [Compare gene set of two input genomes by gene name and start-position proximity](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP08__Annotation_evaluation/PYSCRIPT_Compare_genes_by_name_and_position.py)

```python
python PYSCRIPT_Compare_genes_by_name_and_position.py Limnothrix_sp_HT2024_Bactopia.gb Limnothrix_sp_HT2024_bacass.gb --max-start-diff 500
```

- [Standardize the annotations of a bacterial genome](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP08__Annotation_evaluation/PYSCRIPT_Standardize_annotations_of_bacterial_genome.py)
This script ensures that every `CDS` and every `gene` annotation contain at least a `gene`-tag as well as a `product`-tag. The `gene`-tag contains the four-letter gene abbreviation. The full behaviour of the script is as follows:
```python
python PYSCRIPT_Standardize_annotations_of_bacterial_genome.py input.gb output.gb
```

#### Behavior Table
| Situation                                                                                         | Action                                                                                                                   | Log style                                                             |
| ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| CDS already has valid `gene` and `product`                                                        | Copy missing values to the paired `gene` feature; CDS remains authoritative                                              | White if nothing changes, yellow if synchronization changes something |
| Only the `gene` feature has valid `gene` and `product`                                            | Copy missing values to the paired `CDS` feature                                                                          | Yellow if this changes the CDS, otherwise white                       |
| `gene` and `CDS` disagree                                                                         | Prefer CDS values and record the conflict in the report                                                                  | Yellow if resolved successfully, red if still unresolved              |
| Valid `product` exists but `gene` is missing                                                      | Try local mapping first; otherwise query UniProt **cyanobacteria** by product to infer the most common gene abbreviation | Yellow on successful resolution, red on failure                       |
| Valid `gene` exists but `product` is missing                                                      | Try local mapping first; otherwise query UniProt **cyanobacteria** by gene to infer the most common product description  | Yellow on successful resolution, red on failure                       |
| `standard_name` is present                                                                        | Use it as supporting information during local resolution and as the displayed name in logs                               | Shown as `standard_name` instead of `locus_tag`                       |
| `standard_name` is `hypothetical protein CDS` or `hypothetical protein gene`                      | Do **not** query UniProt and do **not** log the annotation; still standardize `/product` to `hypothetical protein`       | No log output                                                         |
| Product is already `hypothetical protein` for one of those hypothetical-standard-name annotations | Leave it as `hypothetical protein` or rewrite it to the same standard form                                               | No log output                                                         |
| Nothing reliable can be inferred                                                                  | Fall back to unresolved values such as `unknown_gene` or remaining missing information, and record it in the report      | Red                                                                   |
| Annotation does not change                                                                        | Keep existing values as they are                                                                                         | White                                                                 |

#### Priority Order

| Priority | Rule                                                                                                                                                                                              |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1        | Prefer existing CDS qualifiers over gene qualifiers                                                                                                                                               |
| 2        | Copy missing values across paired `gene` and `CDS` features                                                                                                                                       |
| 3        | Apply the local mapping table                                                                                                                                                                     |
| 4        | Query UniProt cyanobacteria by product if `gene` is missing                                                                                                                                       |
| 5        | Query UniProt cyanobacteria by gene if `product` is missing                                                                                                                                       |
| 6        | Skip logging and UniProt lookup for annotations whose `standard_name` is `hypothetical protein CDS` or `hypothetical protein gene`, but still standardize their product to `hypothetical protein` |
| 7        | Record conflicts and unresolved cases in the report                                                                                                                                               |


---

### 09. Gene-Level Visualization
- [Visualization of gene location via GenoVi](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP09a__Visualization_via_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP09a__Visualization_via_GenoVi/FinalAssembly_Bactopia__output_from_GenoVi.png" style="display:block; margin-left:auto; margin-right:auto; width:50%;">

###### **TO DO:** Do GenoVi visualization also for plasmid genome

- [Barchart-style gene tally compatible with GenoVi visualization](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/STEP09b__Gene_tally_of_GenoVi.md)

<img src="https://raw.githubusercontent.com/michaelgruenstaeudl/Cyanobacterial_Genome_Assembly/main/data_STEP09b__Gene_tally_of_GenoVi/genovi_COG_Classification_COG_barplot.png" style="display:block; margin-left:auto; margin-right:auto; width:100%;">


---

### 10. QIIME2 analysis
- [QIIME2 analysis](https://github.com/michaelgruenstaeudl/CyanobacterialGenomeAssemblyAndAnnotation/blob/main/data_STEP10__QIIME2_analysis/)

---

## Notes

- Image files are stored in the corresponding `data_STEPXX__*` directories.
- All workflows assume Linux environments and standard bioinformatics toolchains.
- Individual steps can be adapted to other bacterial genomes with minimal modification.

---
