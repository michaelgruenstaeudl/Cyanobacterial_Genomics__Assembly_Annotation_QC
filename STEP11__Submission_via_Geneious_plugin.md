### Submission of annotated genome to GenBank using the Geneious plugin

#### How to install PGAP locally
##### 1. Ensure you have docker installed and have it in your sudo group
Having it in your sudo group is especially important so 'sudo' does not need to be written by pgap.py

##### 2. Download the PGAP launcher script
```bash
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
```
This downloads the official PGAP launcher script from NCBI’s GitHub repository and saves it locally as pgap.py.

##### 3. Make the script executable
```bash
chmod +x pgap.py
```
This allows the script to be run directly from the command line.

##### 4. Download/update PGAP databases and containers
```bash
./pgap.py --update
```
This installs or updates the PGAP runtime assets, including the required databases and Docker/Singularity container images.

##### 5. Run the built-in test genome
```bash
./pgap.py -r -o mg37_results \
  -g $HOME/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta \
  -s "Mycoplasmoides genitalium"
```
This runs PGAP on the example Mycoplasmoides genitalium genome and writes the results to mg37_results. A successful run verifies that PGAP, its databases, and the container environment are working.

#### Annotations via PGAP
This step runs PGAP using the supplied YAML configuration files (`*_GENOME.yaml`) for the chromosome and plasmid assemblies. It generates standardized prokaryotic genome annotations, including CDS, RNA features, product names, and submission-ready output files based on the input sequences and metadata.

```bash
pgap.py -r Limnothrix_sp_HT2024_chromosome_GENOME.yaml

pgap.py -r Limnothrix_sp_HT2024_plasmid_GENOME.yaml
```

#### Step4__PYSCRIPT_Merge_GenBank_tags.py
This script merges annotation tags and selected feature metadata from a user-defined genome GenBank file (`*_GENOME.gb`) into a PGAP-generated GenBank file (`*_PGAP.gbk`), producing a consolidated output file (`*_MERGED.gb`). It is used to preserve curated user annotations or naming choices while retaining PGAP-standardized structural annotation, and writes a log file documenting the merge process.

```bash
python Step4__PYSCRIPT_Merge_GenBank_tags.py Limnothrix_sp_HT2024_plasmid_PGAP.gbk Limnothrix_sp_HT2024_plasmid_GENOME.gb Limnothrix_sp_HT2024_plasmid_MERGED.gb > Limnothrix_sp_HT2024_plasmid_MERGED.log

python Step4__PYSCRIPT_Merge_GenBank_tags.py Limnothrix_sp_HT2024_chromosome_PGAP.gbk Limnothrix_sp_HT2024_chromosome_GENOME.gb Limnothrix_sp_HT2024_chromosome_MERGED.gb > Limnothrix_sp_HT2024_chromosome_MERGED.log
```

#### Script5_final_corrections_for_GB_submission.py
This script takes the merged GenBank annotation files (`*_MERGED.gb`) and produces cleaned final versions (`*_FINAL.gb`) suitable for GenBank submission via Geneious' GenBank plugin. It standardizes feature annotations, removes formatting or qualifier issues, and applies final compliance corrections to a genome record.

```bash
python Script5_final_corrections_for_GB_submission.py Limnothrix_sp_HT2024_chromosome_MERGED.gb Limnothrix_sp_HT2024_chromosome_FINAL.gb

python Script5_final_corrections_for_GB_submission.py Limnothrix_sp_HT2024_plasmid_MERGED.gb Limnothrix_sp_HT2024_plasmid_FINAL.gb
```

Important: Add a **Genbank submission metadata** block to the genome record within Geneious before conversion to sqn-format via the submission plugin.

#### fix_protein_ids.py
This script edits the ASN.1 submission files (`.asn`) produced via Geneious's GenBank submission plugin to replace invalid protein identifiers with GenBank-compatible gnl|dbname|unique_id style protein IDs. It ensures each CDS-associated protein record has a unique valid identifier and can optionally apply a specified `locus_tag_prefix`.

```bash
python Script5_final_corrections_for_GB_submission.py Limnothrix_sp_HT2024_chromosome_MERGED.gb Limnothrix_sp_HT2024_chromosome_FINAL.gb

python Script5_final_corrections_for_GB_submission.py Limnothrix_sp_HT2024_plasmid_MERGED.gb Limnothrix_sp_HT2024_plasmid_FINAL.gb
```
