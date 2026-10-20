#!/bin/bash
# Usage: 01_download_samples.sh <SRR> <outdir>
set -eo pipefail
ACCESSION="$1"
OUTDIR="${2:-reads}"
mkdir -p "$OUTDIR"

echo ">>> Liens ENA pour $ACCESSION..."
# La réponse ENA contient une colonne run_accession puis fastq_ftp : on prend la DERNIERE colonne de la ligne de données
FASTQ_FTP=$(curl -s --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=fastq_ftp&format=tsv" \
    | awk -F'\t' 'NR==2{print $NF}' | tr -d '\r')

if [ -z "$FASTQ_FTP" ]; then
    echo "!!! Aucun lien ENA pour $ACCESSION - fallback sra-tools"
    prefetch "$ACCESSION" -O "$OUTDIR"
    fasterq-dump "$OUTDIR/$ACCESSION" --split-files --outdir "$OUTDIR"
    gzip -f "$OUTDIR/${ACCESSION}_1.fastq" "$OUTDIR/${ACCESSION}_2.fastq"
    exit 0
fi

IFS=';' read -ra LINKS <<< "$FASTQ_FTP"
for link in "${LINKS[@]}"; do
    fname=$(basename "$link")
    echo ">>> Téléchargement de $fname..."
    wget -q -c --tries=8 --timeout=60 --waitretry=10 "http://${link}" -O "$OUTDIR/$fname"
done
ls -lh "$OUTDIR"/${ACCESSION}*
