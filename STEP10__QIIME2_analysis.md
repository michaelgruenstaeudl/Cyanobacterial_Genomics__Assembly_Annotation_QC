### Conduct QIIME2 analysis and generate visualization across samples

#### Conduct QIIME2 analysis
```bash
bash BASHSCRIPT_QIIME2_analysis.sh
```

#### Generate visualization across samples
```bash
# Visualizing diversity at the taxonomic level of class
python stacked_taxonomy_bars.py . -r c -o HorseThief_metagenomics_rank-class -t 1 --metadata metadata_cyanobacteria_samples.tsv

# Visualizing diversity at the taxonomic level of family but only for cyanobacteria
python stacked_taxonomy_bars.py . -r f -o HorseThief_metagenomics_rank-family_only_cyanos -t 1 --filter-rank c --filter-name Cyanobacteriia --metadata metadata_cyanobacteria_samples.tsv
```
