"""Courbe d'apprentissage, IC bootstrap, confusion par profondeur de séquençage."""
import argparse, os, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score
from common import nested
from load_data import load
ap = argparse.ArgumentParser()
ap.add_argument("--bracken", required=True); ap.add_argument("--map", required=True); ap.add_argument("--s3", required=True)
ap.add_argument("--out", default="results"); ap.add_argument("--repeats", type=int, default=10)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.RandomState(0)
T, M, y, depth = load(a.bracken, a.map, a.s3); n = len(y)

# 1) IC bootstrap sur prédictions hors-pli moyennées sur les répétitions
P = {k: np.zeros(n) for k in ("tax", "met", "fusion")}
for s in range(a.repeats):
    _, p = nested(y, T, M, s)
    for k in P: P[k] += p[k] / a.repeats
ci = {}
for k, v in P.items():
    bs = []
    for _ in range(1000):
        i = rng.randint(0, n, n)
        if len(set(y[i])) == 2: bs.append(roc_auc_score(y[i], v[i]))
    ci[k] = dict(auc=float(roc_auc_score(y, v)), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])

# 2) courbe d'apprentissage (sous-échantillonnage stratifié)
rows = []
for size in sorted({max(8, int(f * n)) for f in (0.4, 0.55, 0.7, 0.85)} | {n}):
    for rep in range(a.repeats):
        idx = []
        for c in (0, 1):
            ic = np.where(y == c)[0]; idx += list(rng.choice(ic, max(2, int(round(size * len(ic) / n))), replace=False))
        idx = np.array(sorted(idx))
        if min(np.bincount(y[idx])) < 4: continue
        r, _ = nested(y[idx], T.iloc[idx], M.iloc[idx], rep)
        rows.append(dict(n=len(idx), rep=rep, **r))
lc = pd.DataFrame(rows).groupby("n").agg(["mean", "std"]); lc.to_csv(f"{a.out}/learning_curve.csv"); print(lc.round(3))

# 3) confusion par profondeur
conf = {}
if depth is not None:
    ld = np.log10(depth)
    conf = dict(auc_depth_alone=float(roc_auc_score(y, ld)), mannwhitney_p=float(mannwhitneyu(ld[y == 1], ld[y == 0]).pvalue),
                spearman_depth_vs_tax_score=float(spearmanr(ld, P["tax"])[0]))
pd.DataFrame({"y":y,"tax":P["tax"],"met":P["met"],"fusion":P["fusion"]}).to_csv(f"{a.out}/oof_predictions.csv",index=False)
json.dump(dict(bootstrap_ci=ci, depth_confounding=conf), open(f"{a.out}/extra.json", "w"), indent=1)
print(json.dumps(dict(bootstrap_ci=ci, depth_confounding=conf), indent=1))
