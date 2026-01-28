#### Installation
```bash
# Install system dependencies for Python and QUAST
sudo apt update
sudo apt install -y python3-full python3-venv git

# Create and activate a Python virtual environment
python3 -m venv quast_env
source quast_env/bin/activate

# Upgrade pip and setuptools in the venv
pip install --upgrade pip setuptools

# Navigate to working directory
cd ~/cyanobacteria_analysis

# Clone QUAST if not already present
if [ ! -d "quast" ]; then
    git clone https://github.com/ablab/quast.git
fi

# Navigate to QUAST directory
cd quast
```

#### Run QUAST on both assemblies
```bash
python quast.py ~/cyanobacteria_analysis/FinalAssembly_Bactopia_corrected.fasta -o ~/cyanobacteria_analysis/quast_Bactopia_report --gene-finding && \
python quast.py ~/cyanobacteria_analysis/FinalAssembly_bacass_corrected.fasta -o ~/cyanobacteria_analysis/quast_Bacass_report --gene-finding
```