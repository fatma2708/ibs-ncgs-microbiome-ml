"""Comparaison de modèles et de références (taxonomie). Inclut une signature publiée À BIAIS (variables choisies sur les mêmes échantillons)."""
import argparse, json, re, numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from common import prep_tax, fit_view
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--s2"); ap.add_argument("--res", default="results"); ap.add_argument("--repeats", type=int, default=10)
a = ap.parse_args(); T, M, y, _ = load(a.bracken, a.map, a.s3)
norm = lambda s: re.sub(r"\s+", " ", str(s).replace("_", " ")).strip().lower()
shannon = -(T.replace(0, np.nan).pipe(lambda d: d * np.log(d))).sum(axis=1).values.reshape(-1, 1)
def cv(fn):
    r = []
    for s in range(a.repeats):
        p = np.zeros(len(y))
        for tr, te in RepeatedStratifiedKFold(n_splits=4, n_repeats=1, random_state=s).split(T, y): p[te] = fn(tr, te, s)
        r.append(roc_auc_score(y, p))
    return float(np.mean(r)), float(np.std(r))
def lr(tr, te, s):
    A, B, _ = prep_tax(T.iloc[tr], T.iloc[te]); return fit_view(A, y[tr], B)
def rf(tr, te, s):
    A, B, _ = prep_tax(T.iloc[tr], T.iloc[te]); return RandomForestClassifier(300, random_state=s).fit(A, y[tr]).predict_proba(B)[:, 1]
def shan(tr, te, s):
    sc = StandardScaler().fit(shannon[tr]); return LogisticRegression().fit(sc.transform(shannon[tr]), y[tr]).predict_proba(sc.transform(shannon[te]))[:, 1]
res = {"Régression logistique L2 (modèle principal)": cv(lr), "Forêt aléatoire": cv(rf), "Diversité de Shannon seule (référence)": cv(shan)}
if a.s2:
    S = pd.read_excel(a.s2, sheet_name="Kraken2_Differential_Bacteria", header=1).iloc[:, :8].dropna(subset=["Taxonomy"])
    S["n"] = S.Taxonomy.map(norm); S["padj"] = pd.to_numeric(S.padj, errors="coerce"); S = S.sort_values("padj")
    cols = {norm(c): c for c in T.columns}; sel = [cols[n] for n in S.n if n in cols][:30]
    def leaky(tr, te, s):
        X = np.log(T[sel].values + 1e-4); X = X - X.mean(1, keepdims=True); return fit_view(X[tr], y[tr], X[te])
    res[f"Signature publiée (top {len(sel)} de S2) — OPTIMISTE, variables choisies sur tous les échantillons"] = cv(leaky)
for k, (m, s) in res.items(): print(f"{k:95s} AUC {m:.3f} ± {s:.3f}")
json.dump({k: dict(auc_mean=m, auc_sd=s) for k, (m, s) in res.items()}, open(f"{a.res}/model_comparison.json", "w"), indent=1, ensure_ascii=False)
