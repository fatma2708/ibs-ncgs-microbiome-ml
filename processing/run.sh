#!/bin/bash
# Usage: run.sh [IBS|NCGS]   -- reprenable, tolérant aux erreurs par échantillon
set -eo pipefail
FILTER="${1:-}"
MAXSEC="${MAXSEC:-39000}"     # ~10h50 : on ne démarre plus d'échantillon après ce délai (limite Kaggle 12h)
START=$(date +%s)
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate ibsncgs
mkdir -p bracken_out kraken2_out reports

# Reprise : récupérer les .bracken d'une version précédente ajoutée en Input
for d in $(find /kaggle/input -type d -name bracken_out 2>/dev/null); do
    cp -n "$d"/*.bracken bracken_out/ 2>/dev/null || true
done

save_results() {
    echo ">>> Sauvegarde des résultats..."
    rm -f kraken2_out/*.output
    zip -qr "results_${FILTER:-ALL}.zip" bracken_out kraken2_out reports phase1_samples.csv failed.txt 2>/dev/null || true
    echo ">>> .bracken présents: $(ls bracken_out/*.bracken 2>/dev/null | wc -l)"
    [ -f failed.txt ] && echo ">>> Échecs: $(sort -u failed.txt | tr '\n' ' ')" || true
}
trap save_results EXIT

process_sample() {
    local RUN="$1"
    bash 01_download_samples.sh "$RUN" reads    || return 1
    bash 02_qc_trimming.sh "$RUN" reads trimmed || return 1
    rm -f reads/${RUN}*
    bash 03_host_removal.sh "$RUN" trimmed clean || return 1
    rm -f trimmed/${RUN}*
    bash 04_taxonomic_profiling.sh "$RUN" clean || return 1
    rm -f clean/${RUN}* kraken2_out/${RUN}.kraken2.output
    return 0
}

CONSEC=0
while IFS=',' read -r RUN GRP REST <&3; do
    RUN=$(echo "$RUN" | tr -d '\r'); GRP=$(echo "$GRP" | tr -d '\r')
    if [ -n "$FILTER" ] && [ "$GRP" != "$FILTER" ]; then continue; fi
    if [ -f "bracken_out/${RUN}.bracken" ]; then echo ">>> $RUN déjà traité, on saute."; continue; fi
    if [ $(( $(date +%s) - START )) -gt "$MAXSEC" ]; then
        echo ">>> Limite de temps atteinte : arrêt propre avant $RUN. Relance pour continuer."; break
    fi
    echo "=============================================="
    echo "  $RUN ($GRP)  -  $(date '+%H:%M:%S')"
    echo "=============================================="
    if process_sample "$RUN" </dev/null; then
        echo ">>> $RUN OK  -  $(date '+%H:%M:%S')"; CONSEC=0
    else
        CONSEC=$((CONSEC+1))
        echo "!!! $RUN ECHEC - on passe au suivant"
        echo "$RUN" >> failed.txt
        rm -f reads/${RUN}* trimmed/${RUN}* clean/${RUN}*
    fi
    df -h /kaggle/working 2>/dev/null | tail -n1 || true
    if [ "$CONSEC" -ge 3 ]; then echo "!!! 3 échecs consécutifs : arrêt (erreur systématique)."; exit 1; fi
done 3< <(tail -n +2 phase1_samples.csv)
echo ">>> FIN DE BOUCLE"
