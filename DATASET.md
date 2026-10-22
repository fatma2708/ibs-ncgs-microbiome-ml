# Dataset documentation

## Source
- **BioProject [PRJNA981497](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA981497)**: stool shotgun metagenomes, Illumina HiSeq X, 2 × 150 bp, about 10 GB per sample (Medgenome Labs).
- **Paper:** Dixit et al. 2026, *mSphere*, "Multi-meta-omics reveal distinct microbial genomic profiles and metabolic dysregulation in non-celiac gluten sensitivity". Supplementary Data S3 holds the LC-MS table.
- **Companion paper (same clinic):** Dixit et al. 2024, *Gut Microbes Reports*, 16S profiling and a shotgun subset (BioProject PRJNA973664 for 16S).
- **Recruitment:** outpatient clinic, Department of Gastroenterology and Human Nutrition, AIIMS New Delhi.

## What the BioProject contains and what this project uses

| Group in SRA metadata | Meaning (per the paper) | Runs in BioProject | Used here |
|---|---|---|---|
| IBS | IBS patients (Rome IV, Bristol chart) | 14 | **14** |
| NCGS | NCGS patients **before** a 6-week gluten-free diet (GFD) | 14 | **14** |
| NCGS_PG | the same NCGS patients **after** the GFD ("post gluten-free") | 14 | **0 (excluded by design)** |

Total in BioProject: 42 runs. **Used: 28 = all IBS + all NCGS runs. There was no sub-sampling or cherry-picking among eligible samples.**

### Why NCGS_PG is excluded
- The paper describes NCGS_PG as the *same patients* sampled again after the diet, analysed for within-patient comparisons. Adding them would duplicate patients across classes and break the independence needed for cross-validation.
- The question asked (IBS vs NCGS) is a between-group question; NCGS_PG answers a different one (diet effect).

### Other exclusions
- No individual sample was dropped for quality or outlier reasons (all 28 processed; depths in the table below).
- Human reads (GRCh38) are removed before classification; *Homo sapiens* is dropped from the profile and abundances are renormalised.
- The model uses only taxa present in at least 20% of training samples, learned inside each training fold.

## Clinical criteria as reported by the authors
*(taken from the papers; verify against the full text before quoting)*
- **IBS:** Rome IV criteria and Bristol stool chart. In the 2024 companion paper the IBS comparison group were anti-gliadin-antibody (AGA)-negative IBS patients.
- **NCGS:** more than 30% symptom reduction on a 6-week GFD, symptom recurrence on gluten re-introduction, and AGA IgA/IgG above 30 U/mL as part of the work-up.
- **Celiac disease ruled out:** negative anti-tTG IgA and no villous atrophy on small-bowel biopsy.
- **Exclusion criteria (patients):** regular alcohol (>20 g/week), smoking (>5/week), abdominal surgery or resection, active NSAID use, antibiotics or prebiotics within 8 weeks, malabsorption syndromes.

## Class counts and depth
- 28 samples, **14 IBS / 14 NCGS** (balanced). 4-fold stratified CV gives 7 test samples per fold (3 or 4 per class).
- Bracken-classified reads per sample (median): IBS 15.6 M, NCGS 14.5 M. Depth alone predicts class with AUC 0.56, so depth is not driving the result.
- Per-sample table: [`data/dataset_table.csv`](data/dataset_table.csv).

## Known limitations of the data
- **No healthy-control group** in this BioProject (the authors say so too). This project compares two patient groups only.
- **Single clinic, single country, single sequencing batch.** Batch and site effects cannot be separated from disease.
- **Clinical covariates** (age, sex, BMI, diet, medication, IBS subtype) were not analysed, so they could confound the result.
- **Metabolomics link is inferred.** The paper does not publish a Mg_ ↔ IBS_k/NCGS_k key. This repo matches them by sorted order within group. It affects only the fusion/metabolomics views, not the taxonomy result, which uses no mapping.
- **Independence of patients** is assumed (one stool sample per patient in the IBS and NCGS groups).
- **Sample size:** with n = 28 the AUC confidence interval is wide (see `results/audit_evaluation.json`).
