"""Analyse différentielle avec correction FDR (Benjamini-Hochberg) sur nos profils retraités + taille d'effet de Klebsiella avec IC bootstrap."""
import argparse, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); a = ap.parse_args()
T, M, y, _ = load(a.bracken, a.map, a.s3); rng = np.random.RandomState(0)
def bh(p):
    p = np.asarray(p); o = np.argsort(p); r = np.empty(len(p)); n = len(p); q = p[o] * n / (np.arange(n) + 1)
    r[o] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1); return r
def table(X, prev=0.2):
    cols = X.columns[(X > 0).mean(axis=0) >= prev]; rows = []
    for c in cols:
        u = mannwhitneyu(X[c][y == 1], X[c][y == 0]); a_ = u.statistic / ((y == 1).sum() * (y == 0).sum())
        rows.append((c, float(X[c][y == 1].median()), float(X[c][y == 0].median()), float(a_), float(u.pvalue)))
    d = pd.DataFrame(rows, columns=["taxon", "median_IBS", "median_NCGS", "auc_IBS_gt_NCGS", "p"]); d["q_BH"] = bh(d.p.values); return d.sort_values("p")
sp = table(T); G = T.T.groupby(T.columns.str.split().str[0]).sum().T; ge = table(G)
sp.to_csv(f"{a.res}/differential_species_FDR.csv", index=False); ge.to_csv(f"{a.res}/differential_genus_FDR.csv", index=False)
print(f"espèces testées {len(sp)} ; q<0.10 : {(sp.q_BH<.1).sum()} ; q<0.05 : {(sp.q_BH<.05).sum()}"); print(sp.head(10).round(4).to_string(index=False))
print(f"\ngenres testés {len(ge)} ; q<0.10 : {(ge.q_BH<.1).sum()}"); print(ge.head(8).round(4).to_string(index=False))
k = G["Klebsiella"].values; pres = {"IBS": float((k[y == 1] > 0).mean()), "NCGS": float((k[y == 0] > 0).mean())}
q = lambda v: np.log10(v + 1e-5)
b = []
for _ in range(2000):
    i = rng.randint(0, len(y), len(y))
    if len(set(y[i])) == 2: b.append(np.median(q(k[i][y[i] == 1])) - np.median(q(k[i][y[i] == 0])))
out = dict(klebsiella_prevalence=pres, klebsiella_log10_median_diff=float(np.median(q(k[y == 1])) - np.median(q(k[y == 0]))), ci95=[float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))])
print("\n", json.dumps(out)); json.dump(out, open(f"{a.res}/klebsiella_effect.json", "w"), indent=1)
