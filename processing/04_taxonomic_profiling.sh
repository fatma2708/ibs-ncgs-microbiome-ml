#!/bin/bash
set -eo pipefail
ACCESSION="$1"; CLEAN_DIR="${2:-clean}"
KRAKEN_DB="${KRAKEN_DB:-$HOME/dbs/kraken2/k2_standard_08gb_20240605}"
mkdir -p kraken2_out bracken_out
echo ">>> Kraken2 $ACCESSION..."
kraken2 --db "$KRAKEN_DB" --paired \
    "${CLEAN_DIR}/${ACCESSION}_R1.clean.fastq.gz" "${CLEAN_DIR}/${ACCESSION}_R2.clean.fastq.gz" \
    --report "kraken2_out/${ACCESSION}.kraken2.report" \
    --output "kraken2_out/${ACCESSION}.kraken2.output" --threads 4
echo ">>> Bracken $ACCESSION..."
bracken -d "$KRAKEN_DB" -i "kraken2_out/${ACCESSION}.kraken2.report" \
    -o "bracken_out/${ACCESSION}.bracken" -r 150 -l S -t 10
