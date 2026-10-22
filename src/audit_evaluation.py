"""Audit of the evaluation procedure: leakage guards, estimand consistency, CI check, permutation calibration.
Usage: python audit_evaluation.py --bracken B --map M --s3 S --res results --out results [--quick]"""
import argparse, json, os, warnings, numpy as np, pandas as pd
from multiprocessing import Pool
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
warnings.filterwarnings("ignore")
from common import nested, prep_tax
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); ap.add_argument("--out", default="results"); ap.add_argument("--quick", action="store_true")
a = ap.parse_args(); T, M, y, _ = load(a.bracken, a.map, a.s3); n = len(y); out = {}
Q = a.quick
# ---- 1. leakage guards -------------------------------------------------------
q_ = n // 4; tr, te1, te2 = np.arange(0, n - 2 * q_), np.arange(n - 2 * q_, n - q_), np.arange(n - q_, n)
k1 = prep_tax(T.iloc[tr], T.iloc[te1])[2]; k2 = prep_tax(T.iloc[tr], T.iloc[te2])[2]
assert list(k1) == list(k2), "feature filter depends on test data"
A, B, _ = prep_tax(T.iloc[tr], T.iloc[te1]); assert np.allclose(A.sum(1), 0) and np.allclose(B.sum(1), 0)
out["guards"] = "prevalence filter uses training rows only; CLR is per-sample; C is chosen by an inner 3-fold CV on training rows only"
# ---- 2. estimands ------------------------------------------------------------
R = 5 if Q else 20; per, P = [], np.zeros(n); fold_auc = []
for s in range(R):
    r, p = nested(y, T, M, s); per.append(r["tax"]); P += p["tax"] / R
    for _, te in StratifiedKFold(4, shuffle=True, random_state=s).split(T, y):
        if len(set(y[te])) == 2: fold_auc.append(roc_auc_score(y[te], p["tax"][te]))
per = np.array(per); ens = roc_auc_score(y, P)
out["estimands"] = dict(per_repeat_pooled_auc_mean=float(per.mean()), per_repeat_sd=float(per.std()), repeat_averaged_oof_auc=float(ens), mean_single_fold_auc=float(np.mean(fold_auc)), sd_single_fold_auc=float(np.std(fold_auc)), n_repeats=R)
# ---- 3. DeLong CI vs bootstrap CI -------------------------------------------
def delong(yv, s):
    pos, neg = s[yv == 1], s[yv == 0]; m_, n_ = len(pos), len(neg)
    psi = lambda p, q: (p > q) + 0.5 * (p == q)
    v10 = np.array([psi(p, neg).mean() for p in pos]); v01 = np.array([psi(pos, q).mean() for q in neg])
    auc = v10.mean(); se = np.sqrt(v10.var(ddof=1) / m_ + v01.var(ddof=1) / n_); return auc, auc - 1.96 * se, min(1, auc + 1.96 * se)
d = delong(y, P); rng = np.random.RandomState(0); bs = []
for _ in range(2000):
    i = rng.randint(0, n, n)
    if len(set(y[i])) == 2: bs.append(roc_auc_score(y[i], P[i]))
out["ci"] = dict(delong=[float(d[1]), float(d[2])], bootstrap_percentile=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                 note="both treat OOF scores as fixed; neither captures model-fitting variance, so true uncertainty is somewhat larger")
# ---- 4. permutation null: matched estimands ---------------------------------
NP = 60 if Q else 200
def single(seed): return nested(np.random.RandomState(seed).permutation(y), T, M, 1000 + seed)[0]["tax"]
def ensem(seed, k=5):
    yp = np.random.RandomState(seed).permutation(y); Pp = np.zeros(n)
    for j in range(k): Pp += nested(yp, T, M, 5000 + seed * 10 + j)[1]["tax"] / k
    return roc_auc_score(yp, Pp)
def nest_p(seed, nperm=40):
    rs = np.random.RandomState(777 + seed); ys = rs.permutation(y)
    obs = nested(ys, T, M, 0)[0]["tax"]; null = [nested(rs.permutation(ys), T, M, 100 + j)[0]["tax"] for j in range(nperm)]
    return (1 + sum(v >= obs for v in null)) / (1 + nperm)
with Pool(max(1, os.cpu_count())) as pool:
    ns = np.array(pool.map(single, range(NP))); ne = np.array(pool.map(ensem, range(NP // 2 if Q else 100)))
    out["permutation"] = dict(single_run_null_mean=float(ns.mean()), single_run_null_sd=float(ns.std()), p_single_run_obs_mean=float((1 + (ns >= per.mean()).sum()) / (1 + len(ns))),
        ensemble5_null_mean=float(ne.mean()), ensemble5_null_sd=float(ne.std()), p_ensemble_obs=float((1 + (ne >= ens).sum()) / (1 + len(ne))), n_perm_single=NP, n_perm_ensemble=len(ne))
    NS = 12 if Q else 40
    pv = np.array(pool.map(nest_p, range(NS)))
out["type1_calibration"] = dict(n_label_shuffled_datasets=NS, perms_each=40, fraction_p_below_010=float((pv <= .10).mean()), fraction_p_below_005=float((pv <= .05).mean()), mean_p=float(pv.mean()), note="uniform p expected: mean ~0.5, fractions ~ threshold")
json.dump(out, open(f"{a.out}/audit_evaluation.json", "w"), indent=1); print(json.dumps(out, indent=1))
