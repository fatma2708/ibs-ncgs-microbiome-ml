# IBS vs NCGS: gut taxonomy + metabolomics (exploratory reanalysis)

Binary classification of IBS vs NCGS (non-celiac gluten sensitivity) from public data
(BioProject PRJNA981497, Dixit et al. 2026, *mSphere*): taxonomic profiles (Kraken2/Bracken, **reprocessed from the raw reads**)
and LC-MS metabolomics (Data S3). CLR transform, nested cross-validation, permutation test, late fusion, SHAP.

> **Limits:** n = 28 (14/14), a single cohort, no external validation. Exploratory results, not a diagnostic tool.
> The metabolomics ↔ sequencing-sample mapping is **inferred** (k-th sorted `Mg_` = `IBS_k`/`NCGS_k`) and not documented by the authors.
> It affects only the metabolomics and fusion views; the taxonomy result uses no mapping.

## Results (n = 28: 14 IBS, 14 NCGS)
| | AUC (95% CI) | permutation p |
|---|---|---|
| Taxonomy (logistic regression, CLR) | 0.80 (0.59–0.95) | 0.02 |
| LC-MS metabolomics | 0.40 (0.19–0.64) | 0.60 |
| Late fusion | 0.79 (0.58–0.94) | 0.02 |

Headline: **AUC ≈ 0.8** for taxonomy (see the audit section for exact definitions and a DeLong CI of 0.64–0.98).

- **Robustness:** 0.79–0.81 across 4 pre-declared configurations (species/genus, prevalence thresholds); random forest 0.77; Shannon diversity alone 0.64; sequencing depth alone 0.56.
- **Agreement with the paper:** 9 of the top 30 features overlap the published DESeq2 analysis (about 3 expected by chance). Enterobacteriaceae as a whole do not differ.
- **Selection-bias control:** a *Klebsiella*/*Rothia* ratio chosen on all the data gives AUC 0.95 (optimistic). When the two taxa are chosen *inside* each training fold (`src/nested_ratio.py`), AUC drops to 0.64 (species, p = 0.13) / 0.68 (genus, p = 0.08). The pair is not robust.
- **Differential abundance (FDR-corrected):** *Anaerostipes* is lower in IBS (q = 0.026, genus level), along with *Rothia*, *Wujia*, *Actinomyces* and *Granulicatella*.
- **Hypothesis retained (exploratory):** *Klebsiella* (mostly *K. pneumoniae*) is higher in IBS than NCGS (about 9× at genus level, p = 0.013 uncorrected, not significant after FDR). To be tested in an independent cohort.

Figures are in `figures/`, raw outputs in `results/`.

## Structure
- `processing/`: read reprocessing (fastp → Bowtie2 GRCh38 → Kraken2 → Bracken), resumable, designed for Kaggle's 12 h CPU limit.
- `src/run_ml.py`: models and evaluation. `src/extra_analyses.py`: learning curve, 95% bootstrap CI, depth-confounding check. `src/make_figures.py`: figures. `src/compare_deseq2.py`: overlap with the published differential taxa (Data S2). `src/compare_models.py`: LR vs random forest vs Shannon baseline vs published signature (biased). `src/nested_ratio.py`: honest evaluation of taxon selection. `src/sensitivity.py`: 4 pre-declared configurations.
- `src/audit_evaluation.py`: audit of the evaluation procedure. `src/freeze_model.py` and `src/validate_external.py`: frozen model and external-validation script. `src/fdr_differential.py` and `src/klebsiella_followup.py`: differential abundance and Klebsiella follow-up.
- `workflow/Snakefile`: the ML step in one command. `docker/`: environments.
- `tests/` + `.github/workflows/ci.yml`: test on synthetic data.
- `DATASET.md`, `data/dataset_table.csv`: dataset documentation.

## Usage
```bash
docker build -f docker/Dockerfile -t ibs-ncgs-ml .
# Data S3: download msphere.00856-25-s0003.xlsx from the paper and place it in data/s3.xlsx
docker run --rm -v $PWD:/app -w /app --entrypoint snakemake ibs-ncgs-ml -c1 --snakefile workflow/Snakefile
# quick test: docker run --rm --entrypoint bash ibs-ncgs-ml tests/test_pipeline.sh
```
The `bracken_out/*.bracken` profiles are produced with `processing/` (the Kraken2 database `k2_standard_08gb_20240605` and GRCh38 must be provided).

## Evaluation audit (`src/audit_evaluation.py`, `results/audit_evaluation.json`)
- **No leakage detected:** the prevalence filter, standardisation and choice of C are learned on training rows only; CLR is per sample (a guard test is included).
- **AUC definition:** mean per repeat = 0.79 (SD 0.03); AUC of predictions averaged over 20 repeats = 0.81. Quote it as "AUC ≈ 0.8".
- **95% CI:** DeLong 0.64–0.98; bootstrap 0.61–0.96. Both treat the out-of-fold scores as fixed, so the true uncertainty is somewhat larger.
- **Permutations:** the null is centred on 0.49 (SD 0.14, so noise is substantial at n = 28). p = 0.015 (single run) and 0.030 (mean of 5 repeats, matched estimand).
- **Test calibration:** on 40 label-shuffled datasets, p ≤ 0.10 in 10% of cases and p ≤ 0.05 in 2.5% (mean p 0.48): no excess false positives.

## External validation
No independent public IBS vs NCGS cohort with shotgun metagenomes was found (see `DATASET.md`). The frozen model (`src/freeze_model.py`, `models/frozen_model_v1.joblib`) and the `src/validate_external.py` script (pre-specified decision rules) are ready for a future cohort. Applying the frozen model to the discovery data gives AUC = 1.0: that is a mechanics check, not a validation.

## Data
See [`DATASET.md`](DATASET.md) and `data/dataset_table.csv` (42 runs in the BioProject, 28 used = all IBS and NCGS; NCGS_PG excluded because it is the same patients after the diet).

## Klebsiella follow-up (`src/klebsiella_followup.py`, exploratory, uncorrected p-values)
- The signal comes mostly from *K. pneumoniae* (median 0.66% vs 0.04%, p = 0.010), with the same direction for *K. quasipneumoniae* and *K. variicola*. *K. oxytoca*, *grimontii* and *michiganensis* do not differ.
- No difference for *Escherichia*, *Enterobacter* or *Citrobacter*: this is not a general bloom of Enterobacteriaceae.
- A wider panel of butyrate producers does not differ (p = 0.21): only *Anaerostipes* stands out. The *Anaerostipes*–*Klebsiella* correlation is negative overall (rho = -0.36) but not within groups: no evidence of a direct link.
