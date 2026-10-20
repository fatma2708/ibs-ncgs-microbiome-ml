"""IBS vs NCGS — taxonomie (Kraken2/Bracken) + métabolomique LC-MS (Data S3). Fusion tardive, CV imbriquée, SHAP.
Usage: python run_ml.py --bracken bracken_out --map sample_id_map.csv --s3 msphere...s0003.xlsx --out results
"""
import argparse, glob, os, re, warnings, json
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
import shap
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--bracken", required=True); ap.add_argument("--map", required=True)
ap.add_argument("--s3", required=True); ap.add_argument("--out", default="results")
ap.add_argument("--repeats", type=int, default=20); ap.add_argument("--perms", type=int, default=200)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rng = np.random.RandomState(0)

# ---------- données ----------
m = pd.read_csv(a.map); m = m[m.Group.isin(["IBS", "NCGS"])].copy()
m["mg"] = m.Mg_name.str.extract(r"(\d+)").astype(int)
m = m.sort_values(["Group", "mg"]); m["k"] = m.groupby("Group").cumcount() + 1
m["met_id"] = m.Group + "_" + m.k.astype(str)   # HYPOTHÈSE: k-ième Mg_ trié == IBS_k / NCGS_k
ngroup = m.groupby("Group").size().to_dict()
runs = {os.path.basename(f)[:-8]: f for f in glob.glob(f"{a.bracken}/*.bracken")}
m = m[m.Run.isin(runs)].copy()
tax = {}
for r in m.Run:
    d = pd.read_csv(runs[r], sep="\t"); tax[r] = d.set_index("name")["fraction_total_reads"]
T = pd.DataFrame(tax).fillna(0).T.drop(columns=["Homo sapiens"], errors="ignore")
T = T.div(T.sum(axis=1), axis=0)
L = pd.read_excel(a.s3, sheet_name="LCMS"); L["v"] = pd.to_numeric(L["Response Ratio"], errors="coerce")
M = L.pivot_table(index="Sample_ID", columns="Metabolite", values="v").fillna(0)
m = m[m.met_id.isin(M.index)]
T = T.loc[m.Run]; M = M.loc[m.met_id]; M.index = m.Run.values
y = (m.Group == "IBS").astype(int).values  # 1=IBS
print("n =", len(y), "IBS", y.sum(), "NCGS", (1 - y).sum())

# ---------- transformations (apprises dans le pli d'entraînement) ----------
def prep_tax(Xtr, Xte, prev=0.2):
    keep = (Xtr > 0).mean(0).values >= prev
    def clr(X):
        X = X.loc[:, keep].values + 1e-4; X = X / X.sum(1, keepdims=True); l = np.log(X); return l - l.mean(1, keepdims=True)
    return clr(Xtr), clr(Xte), Xtr.columns[keep]
def prep_met(Xtr, Xte):
    return np.log1p(Xtr.values), np.log1p(Xte.values), Xtr.columns
def fit_view(Xtr, ytr, Xte):
    sc = StandardScaler().fit(Xtr); A, B = sc.transform(Xtr), sc.transform(Xte)
    best, bs = None, -1
    inner = StratifiedKFold(3, shuffle=True, random_state=1)
    for C in [1e-3, 1e-2, 1e-1, 1]:
        s = []
        for i, j in inner.split(A, ytr):
            if len(set(ytr[i])) < 2 or len(set(ytr[j])) < 2: continue
            s.append(roc_auc_score(ytr[j], LogisticRegression(C=C, max_iter=2000).fit(A[i], ytr[i]).decision_function(A[j])))
        sc_ = np.mean(s) if s else 0
        if sc_ > bs: bs, best = sc_, C
    clf = LogisticRegression(C=best, max_iter=2000).fit(A, ytr)
    return clf.predict_proba(B)[:, 1]
def rank01(p): return pd.Series(p).rank().values / len(p)

def nested(y, T, M, seed, do_perm=False):
    cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=1, random_state=seed)
    P = {k: np.zeros(len(y)) for k in ("tax", "met", "fusion")}
    for tr, te in cv.split(T, y):
        ytr = y[tr]
        A, B, _ = prep_tax(T.iloc[tr], T.iloc[te]); pt = fit_view(A, ytr, B)
        A, B, _ = prep_met(M.iloc[tr], M.iloc[te]); pm = fit_view(A, ytr, B)
        P["tax"][te] = pt; P["met"][te] = pm; P["fusion"][te] = (pt + pm) / 2
    return {k: roc_auc_score(y, v) for k, v in P.items()}, P

res = []
for s in range(a.repeats):
    r, _ = nested(y, T, M, s); res.append(r)
R = pd.DataFrame(res); print(R.agg(["mean", "std"]).round(3))
R.to_csv(f"{a.out}/auc_repeats.csv", index=False)

# permutation (1 répétition par permutation)
null = {k: [] for k in R.columns}
for p in range(a.perms):
    yp = rng.permutation(y); r, _ = nested(yp, T, M, 1000 + p)
    for k in null: null[k].append(r[k])
summ = {}
for k in R.columns:
    obs = R[k].mean(); pv = (1 + np.sum(np.array(null[k]) >= obs)) / (1 + a.perms)
    summ[k] = dict(auc_mean=float(obs), auc_sd=float(R[k].std()), perm_p=float(pv))
json.dump(summ, open(f"{a.out}/summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))

# ---------- SHAP (ajustement sur tous les échantillons — descriptif seulement) ----------
for name, prep, X in (("taxonomy", prep_tax, T), ("metabolomics", prep_met, M)):
    A, _, cols = prep(X, X); sc = StandardScaler().fit(A); A = sc.transform(A)
    clf = LogisticRegression(C=0.01, max_iter=2000).fit(A, y)
    sv = shap.LinearExplainer(clf, A).shap_values(A)
    imp = pd.DataFrame({"feature": cols, "mean_abs_shap": np.abs(sv).mean(0), "mean_shap_IBS": sv[y == 1].mean(0)}).sort_values("mean_abs_shap", ascending=False)
    imp.to_csv(f"{a.out}/shap_{name}.csv", index=False)
    print(name, "top10:", imp.feature.head(10).tolist())
