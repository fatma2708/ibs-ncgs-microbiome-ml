#!/bin/bash
set -eo pipefail
ACCESSION="$1"; READS_DIR="${2:-reads}"; TRIM_DIR="${3:-trimmed}"
mkdir -p "$TRIM_DIR" reports/fastp
R1=$(ls "${READS_DIR}/${ACCESSION}"*_1.fastq.gz 2>/dev/null | head -n1 || true)
R2=$(ls "${READS_DIR}/${ACCESSION}"*_2.fastq.gz 2>/dev/null | head -n1 || true)
if [ -z "$R1" ] || [ -z "$R2" ]; then echo "!!! R1/R2 introuvables pour $ACCESSION"; exit 1; fi
echo ">>> fastp $ACCESSION..."
fastp -i "$R1" -I "$R2" \
    -o "${TRIM_DIR}/${ACCESSION}_1.trim.fastq.gz" -O "${TRIM_DIR}/${ACCESSION}_2.trim.fastq.gz" \
    --detect_adapter_for_pe --thread 4 \
    --json "reports/fastp/${ACCESSION}.fastp.json" --html "reports/fastp/${ACCESSION}.fastp.html"
