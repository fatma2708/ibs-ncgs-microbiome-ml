# IBS vs NCGS — taxonomie + métabolomique (réanalyse exploratoire)

Classification binaire IBS vs NCGS (sensibilité au gluten non cœliaque) à partir de données publiques
(BioProject PRJNA981497, Dixit et al. 2026, *mSphere*) : profils taxonomiques (Kraken2/Bracken, retraités depuis les lectures brutes)
et métabolomique (Data S3). CLR, validation croisée imbriquée, test de permutation, fusion tardive, SHAP.

> **Limites :** n = 28 (14/14), une seule cohorte, pas de validation externe. Résultats exploratoires, pas un outil diagnostique.
> L'appariement métabolomique ↔ échantillons de séquençage est **déduit** (k-ième `Mg_` trié = `IBS_k`/`NCGS_k`), non documenté par les auteurs.

## Résultats (n = 28 : 14 IBS, 14 NCGS)
| | AUC (IC 95 %) | permutation p |
|---|---|---|
| Taxonomie (régression logistique, CLR) | 0,80 (0,59–0,95) | 0,02 |
| Métabolomique LC-MS | 0,40 (0,19–0,64) | 0,60 |
| Fusion tardive | 0,79 (0,58–0,94) | 0,02 |

Robustesse : 0,79–0,81 sur 4 configurations (espèce/genre, seuils de prévalence) ; forêt aléatoire 0,77 ; diversité de Shannon seule 0,64 ; profondeur de séquençage seule 0,56.
9 des 30 variables principales recoupent l'analyse DESeq2 publiée (≈3 attendues par hasard). Les Enterobacteriaceae en bloc ne diffèrent pas.
**Contrôle de biais de sélection :** un rapport *Klebsiella*/*Rothia* choisi sur toutes les données donne AUC 0,95 (optimiste). Quand les deux taxons sont choisis *à l'intérieur* de chaque pli d'entraînement (`src/nested_ratio.py`), l'AUC tombe à 0,64 (espèce, p = 0,13) / 0,68 (genre, p = 0,08). *Klebsiella* est le taxon le plus souvent sélectionné comme enrichi dans l'IBS ; le couple n'est pas robuste. Hypothèse exploratoire retenue : *Klebsiella* plus élevé dans l'IBS, à tester sur une cohorte indépendante.
Figures dans `figures/`, sorties brutes dans `results/`.

## Structure
- `processing/` — retraitement des lectures (fastp → Bowtie2 GRCh38 → Kraken2 → Bracken), reprenable, conçu pour 12 h de CPU Kaggle.
- `src/run_ml.py` — modèles et évaluation ; `src/extra_analyses.py` — courbe d'apprentissage, IC bootstrap 95 %, contrôle de la profondeur ; `src/make_figures.py` — figures ; `src/compare_deseq2.py` — recoupement avec les taxons différentiels publiés (Data S2) ; `src/compare_models.py` — LR vs forêt aléatoire vs référence (Shannon) vs signature publiée (biaisée) ; `src/nested_ratio.py` — évaluation honnête de la sélection de taxons ; `src/sensitivity.py` — 4 configurations pré-déclarées (espèce/genre, seuils de prévalence).
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

## Audit de l'évaluation (`src/audit_evaluation.py`, `results/audit_evaluation.json`)
- **Pas de fuite détectée** : filtre de prévalence, standardisation et choix de C appris sur les lignes d'entraînement seulement ; CLR par échantillon (test de garde inclus).
- **Définition de l'AUC** : moyenne par répétition = 0,79 (SD 0,03) ; AUC des prédictions moyennées sur 20 répétitions = 0,81. À citer comme « AUC ≈ 0,8 ».
- **IC 95 %** : DeLong 0,64–0,98 ; bootstrap 0,61–0,96. Les deux traitent les scores hors-pli comme fixes, donc l'incertitude réelle est un peu plus grande.
- **Permutations** : nulle centrée sur 0,49 (SD 0,14, donc bruit important à n = 28). p = 0,015 (une exécution) et 0,030 (moyenne de 5 répétitions, estimand apparié).
- **Calibration du test** : sur 40 jeux à étiquettes mélangées, p ≤ 0,10 dans 10 % des cas et p ≤ 0,05 dans 2,5 % (p moyen 0,48) : pas de faux positifs en excès.

## Validation externe
Aucune cohorte publique indépendante IBS vs NCGS avec métagénomes shotgun n'a été trouvée (voir `DATASET.md`). Le modèle gelé (`src/freeze_model.py`, `models/frozen_model_v1.joblib`) et le script `src/validate_external.py` (règles de décision prédéfinies) sont prêts pour une cohorte future. Appliquer le modèle gelé aux données de découverte donne AUC = 1,0 : c'est un test de mécanique, pas une validation.

## Données
Voir [`DATASET.md`](DATASET.md) et `data/dataset_table.csv` (42 runs dans le BioProject, 28 utilisés = tous les IBS et NCGS ; NCGS_PG exclu car mêmes patients après régime).

## Suivi de l'hypothèse Klebsiella (`src/klebsiella_followup.py`, exploratoire, p non corrigés)
- Le signal vient surtout de *K. pneumoniae* (médiane 0,66 % vs 0,04 %, p = 0,010), avec le même sens pour *K. quasipneumoniae* et *K. variicola*. *K. oxytoca*, *grimontii* et *michiganensis* ne diffèrent pas.
- Pas de différence pour *Escherichia*, *Enterobacter*, *Citrobacter* : ce n'est pas une prolifération générale d'entérobactéries.
- Un panel plus large de producteurs de butyrate ne diffère pas (p = 0,21) : seul *Anaerostipes* ressort. La corrélation *Anaerostipes*–*Klebsiella* est négative globalement (ρ = −0,36) mais pas au sein des groupes : pas de preuve d'un lien direct.
