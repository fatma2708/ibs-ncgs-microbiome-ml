"""Fonctions partagées."""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

# ---------- transformations (apprises dans le pli d'entraînement) ----------
def prep_tax(Xtr, Xte, prev=0.2):
    keep = (Xtr > 0).mean(0).values >= prev
    def clr(X):
        X = X.loc[:, keep].values + 1e-4; X = X / X.sum(1, keepdims=True); l = np.log(X); return l - l.mean(1, keepdims=True)
    return clr(Xtr), clr(Xte), Xtr.columns[keep]
def prep_met(Xtr, Xte):
    return np.log1p(Xtr.values), np.log1p(Xte.values), Xtr.columns
def fit_model(Xtr, ytr):
    """Standardise, choose C by inner 3-fold CV (training rows only), refit. Returns (scaler, classifier, C)."""
    sc = StandardScaler().fit(Xtr); A = sc.transform(Xtr)
    best, bs = None, -1
    inner = StratifiedKFold(3, shuffle=True, random_state=1)
    for C in [1e-3, 1e-2, 1e-1, 1]:
        s = []
        for i, j in inner.split(A, ytr):
            if len(set(ytr[i])) < 2 or len(set(ytr[j])) < 2: continue
            s.append(roc_auc_score(ytr[j], LogisticRegression(C=C, max_iter=2000).fit(A[i], ytr[i]).decision_function(A[j])))
        sc_ = np.mean(s) if s else 0
        if sc_ > bs: bs, best = sc_, C
    return sc, LogisticRegression(C=best, max_iter=2000).fit(A, ytr), best
def fit_view(Xtr, ytr, Xte):
    sc, clf, _ = fit_model(Xtr, ytr)
    return clf.predict_proba(sc.transform(Xte))[:, 1]

def nested(y, T, M, seed, do_perm=False):
    cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=1, random_state=seed)
    P = {k: np.zeros(len(y)) for k in ("tax", "met", "fusion")}
    for tr, te in cv.split(T, y):
        ytr = y[tr]
        A, B, _ = prep_tax(T.iloc[tr], T.iloc[te]); pt = fit_view(A, ytr, B)
        A, B, _ = prep_met(M.iloc[tr], M.iloc[te]); pm = fit_view(A, ytr, B)
        P["tax"][te] = pt; P["met"][te] = pm; P["fusion"][te] = (pt + pm) / 2
    return {k: roc_auc_score(y, v) for k, v in P.items()}, P

