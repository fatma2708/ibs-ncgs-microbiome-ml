"""Évaluation honnête de la procédure 'trouver un rapport de deux taxons' : la sélection est refaite dans chaque pli d'entraînement."""
import argparse, json, collections, numpy as np
from scipy.stats import mannwhitneyu
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); ap.add_argument("--level", choices=["species", "genus"], default="species"); ap.add_argument("--repeats", type=int, default=20); ap.add_argument("--perms", type=int, default=200)
a = ap.parse_args(); T, M, y, _ = load(a.bracken, a.map, a.s3); rng = np.random.RandomState(0)
if a.level == "genus": T = T.T.groupby(T.columns.str.split().str[0]).sum().T
X = T.values; names = np.array(T.columns); EPS = 1e-5
def select(Xtr, ytr, prev=0.2):
    keep = np.where((Xtr > 0).mean(0) >= prev)[0]; up, dn = [], []
    for j in keep:
        p = mannwhitneyu(Xtr[ytr == 1, j], Xtr[ytr == 0, j]).pvalue
        d = np.median(Xtr[ytr == 1, j]) - np.median(Xtr[ytr == 0, j]) if np.median(Xtr[ytr == 1, j]) != np.median(Xtr[ytr == 0, j]) else Xtr[ytr == 1, j].mean() - Xtr[ytr == 0, j].mean()
        (up if d > 0 else dn).append((p, j))
    return min(up)[1], min(dn)[1]
def run(yv, seed, track=None):
    P = np.zeros(len(yv))
    for tr, te in RepeatedStratifiedKFold(n_splits=4, n_repeats=1, random_state=seed).split(X, yv):
        i, j = select(X[tr], yv[tr]); P[te] = np.log10((X[te, i] + EPS) / (X[te, j] + EPS))
        if track is not None: track[(names[i], names[j])] += 1
    return roc_auc_score(yv, P)
track = collections.Counter(); aucs = [run(y, s, track) for s in range(a.repeats)]
null = [run(rng.permutation(y), 1000 + p) for p in range(a.perms)]
pv = (1 + np.sum(np.array(null) >= np.mean(aucs))) / (1 + a.perms); tot = sum(track.values())
top = [dict(taxon_up=u, taxon_down=d, share=c / tot) for (u, d), c in track.most_common(6)]
out = dict(auc_mean=float(np.mean(aucs)), auc_sd=float(np.std(aucs)), perm_p=float(pv), most_selected_pairs=top)
print(json.dumps(out, indent=1, ensure_ascii=False)); json.dump(out, open(f"{a.res}/nested_ratio_{a.level}.json", "w"), indent=1, ensure_ascii=False)
