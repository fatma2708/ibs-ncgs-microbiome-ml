# IBS vs NCGS — taxonomie + métabolomique (réanalyse exploratoire)

Classification binaire IBS vs NCGS (sensibilité au gluten non cœliaque) à partir de données publiques
(BioProject PRJNA981497, Dixit et al. 2026, *mSphere*) : profils taxonomiques (Kraken2/Bracken, retraités depuis les lectures brutes)
et métabolomique (Data S3). CLR, validation croisée imbriquée, test de permutation, fusion tardive, SHAP.

> **Limites :** n = 28 (14/14), une seule cohorte, pas de validation externe. Résultats exploratoires, pas un outil diagnostique.
> L'appariement métabolomique ↔ échantillons de séquençage est **déduit** (k-ième `Mg_` trié = `IBS_k`/`NCGS_k`), non documenté par les auteurs.

## Structure
- `processing/` — retraitement des lectures (fastp → Bowtie2 GRCh38 → Kraken2 → Bracken), reprenable, conçu pour 12 h de CPU Kaggle.
- `src/run_ml.py` — modèles et évaluation.
- `workflow/Snakefile` — étape ML en une commande. `docker/` — environnements.
- `tests/` + `.github/workflows/ci.yml` — test sur données synthétiques.

## Utilisation
```bash
docker build -f docker/Dockerfile -t ibs-ncgs-ml .
# Data S3 : télécharger msphere.00856-25-s0003.xlsx depuis l'article, le placer dans data/s3.xlsx
docker run --rm -v $PWD:/app -w /app --entrypoint snakemake ibs-ncgs-ml -c1 --snakefile workflow/Snakefile
# test rapide : docker run --rm --entrypoint bash ibs-ncgs-ml tests/test_pipeline.sh
```
Les profils `bracken_out/*.bracken` s'obtiennent avec `processing/` (bases Kraken2 `k2_standard_08gb_20240605` et GRCh38 à fournir).
