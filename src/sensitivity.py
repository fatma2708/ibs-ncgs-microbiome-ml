"""Analyse de sensibilité PRÉ-DÉCLARÉE (4 configurations, pas de recherche d'optimum) — taxonomie seule."""
import argparse, json, numpy as np, pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from common import prep_tax, fit_view
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); ap.add_argument("--repeats", type=int, default=10)
a = ap.parse_args()
T, M, y, _ = load(a.bracken, a.map, a.s3)
G = T.T.groupby(T.columns.str.split().str[0]).sum().T  # niveau genre
configs = {"espèce, prévalence ≥10 %": (T, .1), "espèce, prévalence ≥20 % (principal)": (T, .2), "espèce, prévalence ≥40 %": (T, .4), "genre, prévalence ≥20 %": (G, .2)}
out = {}
for name, (X, prev) in configs.items():
    aucs = []
    for s in range(a.repeats):
        p = np.zeros(len(y))
        for tr, te in RepeatedStratifiedKFold(n_splits=4, n_repeats=1, random_state=s).split(X, y):
            A, B, _ = prep_tax(X.iloc[tr], X.iloc[te], prev); p[te] = fit_view(A, y[tr], B)
        aucs.append(roc_auc_score(y, p))
    out[name] = dict(auc_mean=float(np.mean(aucs)), auc_sd=float(np.std(aucs)), n_features=int(((X > 0).mean(0) >= prev).sum()))
    print(f"{name:42s} AUC {np.mean(aucs):.3f} ± {np.std(aucs):.3f}  ({out[name]['n_features']} variables)")
json.dump(out, open(f"{a.res}/sensitivity.json", "w"), indent=1, ensure_ascii=False)
