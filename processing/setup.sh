#!/bin/bash
set -eo pipefail
echo ">>> miniforge..."
if [ ! -d "$HOME/miniforge3" ]; then
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
fi
source "$HOME/miniforge3/etc/profile.d/conda.sh"
if ! conda env list | grep -q ibsncgs; then
    conda create -y -n ibsncgs python=3.10
    conda activate ibsncgs
    conda install -y -c bioconda -c conda-forge sra-tools fastp bowtie2 kraken2 bracken pandas scikit-learn shap
else
    conda activate ibsncgs
fi

echo ">>> Base Kraken2 (~8 Go)..."
mkdir -p ~/dbs/kraken2/k2_standard_08gb_20240605
if [ ! -f ~/dbs/kraken2/k2_standard_08gb_20240605/hash.k2d ]; then
    cd ~/dbs/kraken2
    wget -q --tries=5 --timeout=120 https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz
    tar -xzf k2_standard_08gb_20240605.tar.gz -C k2_standard_08gb_20240605
    rm -f k2_standard_08gb_20240605.tar.gz
fi

echo ">>> Index Bowtie2 GRCh38 (~4 Go)..."
mkdir -p ~/dbs/bowtie2_grch38
if [ -z "$(find ~/dbs/bowtie2_grch38 -name '*.1.bt2' ! -name '*.rev.1.bt2' 2>/dev/null | head -n1)" ]; then
    cd ~/dbs/bowtie2_grch38
    wget -q --tries=5 --timeout=120 https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
    unzip -q -o GRCh38_noalt_as.zip
    rm -f GRCh38_noalt_as.zip
fi
echo ">>> Index Bowtie2 détecté : $(find ~/dbs/bowtie2_grch38 -name '*.1.bt2' ! -name '*.rev.1.bt2' | head -n1)"
[ -n "$(find ~/dbs/bowtie2_grch38 -name '*.1.bt2' ! -name '*.rev.1.bt2' | head -n1)" ] || { echo "!!! Index Bowtie2 absent après setup"; exit 1; }
echo ">>> Setup terminé."
df -h ~ | tail -n1
