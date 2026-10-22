"""Follow-up on the Klebsiella hypothesis (exploratory, uncorrected p-values): which species, specificity vs other Enterobacteriaceae,
butyrate-producer panel, and Anaerostipes-Klebsiella correlation. Usage: python klebsiella_followup.py --bracken B --map M --s3 S --res results"""
import argparse, json, numpy as np, pandas as pd, warnings
from scipy.stats import spearmanr, mannwhitneyu
from load_data import load
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); a = ap.parse_args()
T, M, y, _ = load(a.bracken, a.map, a.s3); G = T.T.groupby(T.columns.str.split().str[0]).sum().T * 100
mw = lambda v: float(mannwhitneyu(v[y == 1], v[y == 0]).pvalue)
out = {"klebsiella_species": {}, "other_enterobacteriaceae": {}}
for c in [c for c in T.columns if c.startswith("Klebsiella")]:
    v = T[c] * 100
    if (v > 0).mean() >= .2: out["klebsiella_species"][c] = dict(median_IBS=float(v[y == 1].median()), median_NCGS=float(v[y == 0].median()), prevalence=float((v > 0).mean()), p_uncorrected=mw(v))
for g in ("Escherichia", "Enterobacter", "Citrobacter"):
    if g in G: out["other_enterobacteriaceae"][g] = dict(median_IBS=float(G[g][y == 1].median()), median_NCGS=float(G[g][y == 0].median()), p_uncorrected=mw(G[g]))
panel = [g for g in ["Anaerostipes", "Roseburia", "Faecalibacterium", "Coprococcus", "Butyricicoccus", "Eubacterium", "Agathobacter", "Subdoligranulum"] if g in G]
bp = G[panel].sum(1); out["butyrate_panel"] = dict(genera=panel, median_IBS=float(bp[y == 1].median()), median_NCGS=float(bp[y == 0].median()), p_uncorrected=mw(bp))
kl = G["Klebsiella"]; r, p = spearmanr(G["Anaerostipes"], kl); out["anaerostipes_vs_klebsiella"] = dict(rho_all=float(r), p_all=float(p),
    rho_IBS=float(spearmanr(G["Anaerostipes"][y == 1], kl[y == 1])[0]), rho_NCGS=float(spearmanr(G["Anaerostipes"][y == 0], kl[y == 0])[0]),
    note="negative overall but not within groups: no support for a direct link")
json.dump(out, open(f"{a.res}/klebsiella_followup.json", "w"), indent=1); print(json.dumps(out, indent=1))
