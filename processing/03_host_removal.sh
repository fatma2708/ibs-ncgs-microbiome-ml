#!/bin/bash
set -eo pipefail
ACCESSION="$1"; TRIM_DIR="${2:-trimmed}"; CLEAN_DIR="${3:-clean}"
# Index auto-détecté : on cherche *.1.bt2 sous ~/dbs/bowtie2_grch38 (la structure du zip peut varier)
BT2_INDEX="${BT2_INDEX:-$(find "$HOME/dbs/bowtie2_grch38" -name '*.1.bt2' ! -name '*.rev.1.bt2' 2>/dev/null | head -n1 | sed 's/\.1\.bt2$//')}"
mkdir -p "$CLEAN_DIR" reports/bowtie2
if [ -z "$BT2_INDEX" ] || [ ! -f "${BT2_INDEX}.1.bt2" ]; then
    echo "!!! Index Bowtie2 introuvable: ${BT2_INDEX}.1.bt2"; find "$HOME/dbs/bowtie2_grch38" -maxdepth 3 | head -20; exit 1; fi
echo ">>> Bowtie2 (retrait humain) $ACCESSION..."
bowtie2 -x "$BT2_INDEX" \
    -1 "${TRIM_DIR}/${ACCESSION}_1.trim.fastq.gz" -2 "${TRIM_DIR}/${ACCESSION}_2.trim.fastq.gz" \
    --un-conc-gz "${CLEAN_DIR}/${ACCESSION}_hostremoved_%.fastq.gz" \
    -S /dev/null -p 4 2> "reports/bowtie2/${ACCESSION}.bowtie2.log" || {
    echo "!!! Bowtie2 a échoué. Fin du log :"; tail -n 40 "reports/bowtie2/${ACCESSION}.bowtie2.log"; exit 1; }
if [ ! -f "${CLEAN_DIR}/${ACCESSION}_hostremoved_1.fastq.gz" ]; then
    echo "!!! Sortie Bowtie2 introuvable. Contenu de ${CLEAN_DIR} :"; ls -la "$CLEAN_DIR"
    echo "--- log :"; tail -n 20 "reports/bowtie2/${ACCESSION}.bowtie2.log"; exit 1; fi
mv "${CLEAN_DIR}/${ACCESSION}_hostremoved_1.fastq.gz" "${CLEAN_DIR}/${ACCESSION}_R1.clean.fastq.gz"
mv "${CLEAN_DIR}/${ACCESSION}_hostremoved_2.fastq.gz" "${CLEAN_DIR}/${ACCESSION}_R2.clean.fastq.gz"
echo ">>> $(grep 'overall alignment rate' reports/bowtie2/${ACCESSION}.bowtie2.log || echo 'taux N/A')"
