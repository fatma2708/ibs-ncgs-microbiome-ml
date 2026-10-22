"""Recoupement entre les variables du modèle et les taxons différentiels publiés (Data S2, IBS vs NCGS, DESeq2/Kraken2)."""
import argparse, json, re, numpy as np, pandas as pd
from scipy.stats import hypergeom, spearmanr, mannwhitneyu
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3", "s2"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); ap.add_argument("--top", type=int, default=30)
a = ap.parse_args()
norm = lambda s: re.sub(r"\s+", " ", str(s).replace("_", " ")).strip().lower()
S = pd.read_excel(a.s2, sheet_name="Kraken2_Differential_Bacteria", header=1).iloc[:, :8].dropna(subset=["Taxonomy"])
S["n"] = S.Taxonomy.map(norm); S["padj"] = pd.to_numeric(S.padj, errors="coerce"); S["lfc"] = pd.to_numeric(S.log2FoldChange, errors="coerce")
sig = S[S.padj < 0.05]; print(f"S2: {len(S)} taxons (tous significatifs, padj<0.05)")
T, M, y, _ = load(a.bracken, a.map, a.s3)
shap = pd.read_csv(f"{a.res}/shap_taxonomy.csv"); shap["n"] = shap.feature.map(norm)
univ = pd.DataFrame({"n": [norm(c) for c in T.columns], "p": [mannwhitneyu(T[c][y == 1], T[c][y == 0]).pvalue if T[c].nunique() > 1 else 1 for c in T.columns],
                     "d": [np.log2((T[c][y == 1].mean() + 1e-6) / (T[c][y == 0].mean() + 1e-6)) for c in T.columns]})
# S2 ne liste que les taxons significatifs : univers = taxons testables par le modèle (prévalence >= 20 %)
tested = set(univ.n[(T > 0).mean(0).values >= 0.2]); universe = tested; sigset = set(sig.n) & universe
top = set(shap.head(a.top).n) & universe
k = len(top & sigset); N, K, n = len(universe), len(sigset), len(top)
pv = float(hypergeom.sf(k - 1, N, K, n)) if n else float("nan")
print(f"taxons testés: {N}; dans S2 (significatifs): {K}; top-{a.top} du modèle: {len(top)}; recoupement: {k}; p(hypergéo)={pv:.3g}")
print("recoupement:", sorted(top & sigset))
mm = univ[univ.n.isin(tested)].merge(S[["n", "lfc"]].drop_duplicates("n"), on="n")
rho = spearmanr(mm.d, mm.lfc)[0]
print(f"corrélation de Spearman entre log2FC observé (IBS/NCGS) et log2FC publié (n={len(mm)}): {rho:.2f}  (signe dépendant de la convention de contraste ; mêmes échantillons : contrôle croisé des pipelines, pas une réplication)")
json.dump(dict(n_common=N, n_sig_s2=K, n_top=n, overlap=k, overlap_taxa=sorted(top & sigset), hypergeom_p=pv, spearman_lfc=float(rho)), open(f"{a.res}/deseq2_comparison.json", "w"), indent=1)
