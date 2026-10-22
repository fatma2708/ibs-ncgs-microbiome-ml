"""Test direct de l'hypothèse : IBS enrichi en Enterobacteriaceae vs NCGS (abondance relative cumulée, Mann-Whitney)."""
import argparse, json, numpy as np
from scipy.stats import mannwhitneyu
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); a = ap.parse_args()
T, M, y, _ = load(a.bracken, a.map, a.s3)
G = {"Escherichia", "Shigella", "Salmonella", "Klebsiella", "Enterobacter", "Citrobacter", "Cronobacter", "Kluyvera", "Raoultella",
     "Leclercia", "Lelliottia", "Pantoea", "Kosakonia", "Atlantibacter", "Buttiauxella", "Franconibacter", "Phytobacter", "Superficieibacter"}
cols = [c for c in T.columns if c.split()[0] in G]; ab = T[cols].sum(axis=1).values
u = mannwhitneyu(ab[y == 1], ab[y == 0], alternative="greater")
auc = u.statistic / ((y == 1).sum() * (y == 0).sum())
out = dict(n_species=len(cols), median_IBS=float(np.median(ab[y == 1])), median_NCGS=float(np.median(ab[y == 0])),
           mean_IBS=float(ab[y == 1].mean()), mean_NCGS=float(ab[y == 0].mean()), p_one_sided=float(u.pvalue), auc=float(auc))
print(json.dumps(out, indent=1)); json.dump(out, open(f"{a.res}/hypothesis_enterobacteriaceae.json", "w"), indent=1)
