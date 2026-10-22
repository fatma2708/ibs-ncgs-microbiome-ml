"""Apply the FROZEN model to an independent cohort of Bracken species profiles. No refitting, no tuning.
Pre-specified outputs: (1) AUC of the frozen IBS-likeness score with DeLong CI and permutation p; (2) directional tests
for the two hypotheses (Anaerostipes lower in IBS, Klebsiella higher in IBS), one-sided Mann-Whitney.
labels.csv columns: Run,Group  (Group in {IBS, <comparator>}; comparator may be NCGS or healthy controls)
Usage: python validate_external.py --model models/frozen_model_v1.joblib --bracken DIR --labels labels.csv --positive IBS --out results"""
import argparse, glob, json, os, joblib, numpy as np, pandas as pd, warnings
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True); ap.add_argument("--bracken", required=True); ap.add_argument("--labels", required=True)
ap.add_argument("--positive", default="IBS"); ap.add_argument("--out", default="results"); ap.add_argument("--perms", type=int, default=10000)
a = ap.parse_args(); md = joblib.load(a.model); lab = pd.read_csv(a.labels)
prof = {}
for r in lab.Run:
    f = f"{a.bracken}/{r}.bracken"
    if os.path.exists(f): d = pd.read_csv(f, sep="\t"); prof[r] = d.set_index("name")["fraction_total_reads"]
X = pd.DataFrame(prof).fillna(0).T.drop(columns=["Homo sapiens"], errors="ignore"); X = X.div(X.sum(1), axis=0)
lab = lab[lab.Run.isin(X.index)].reset_index(drop=True); X = X.loc[lab.Run]; y = (lab.Group == a.positive).astype(int).values
cov = len(set(md["features"]) & set(X.columns)) / len(md["features"])
Z = X.reindex(columns=md["features"]).fillna(0).values + md["pseudocount"]; Z = Z / Z.sum(1, keepdims=True); L = np.log(Z); L = L - L.mean(1, keepdims=True)
score = md["clf"].decision_function(md["scaler"].transform(L))
def delong(yv, s):
    pos, neg = s[yv == 1], s[yv == 0]; psi = lambda p, q: (p > q) + 0.5 * (p == q)
    v10 = np.array([psi(p, neg).mean() for p in pos]); v01 = np.array([psi(pos, q).mean() for q in neg]); auc = v10.mean()
    se = np.sqrt(v10.var(ddof=1) / len(pos) + v01.var(ddof=1) / len(neg)); return auc, max(0, auc - 1.96 * se), min(1, auc + 1.96 * se)
auc, lo, hi = delong(y, score); rng = np.random.RandomState(0)
null = np.array([roc_auc_score(rng.permutation(y), score) for _ in range(a.perms)]); p = (1 + (null >= auc).sum()) / (1 + a.perms)
G = X.T.groupby(X.columns.str.split().str[0]).sum().T; hyp = {}
for g, alt in (("Anaerostipes", "less"), ("Klebsiella", "greater")):
    if g in G: hyp[g] = dict(direction_predicted="lower in IBS" if alt == "less" else "higher in IBS", one_sided_p=float(mannwhitneyu(G[g][y == 1], G[g][y == 0], alternative=alt).pvalue),
                             median_IBS=float(G[g][y == 1].median()), median_comparator=float(G[g][y == 0].median()))
out = dict(n=int(len(y)), n_positive=int(y.sum()), frozen_model=os.path.basename(a.model), feature_coverage=float(cov), auc=float(auc), ci95_delong=[float(lo), float(hi)], perm_p=float(p), hypotheses=hyp,
           decision_rule="Model transfers if AUC CI lower bound > 0.5 and perm p < 0.05. Hypotheses supported if one-sided p < 0.05 in the predicted direction; Bonferroni 0.025 for two tests.")
os.makedirs(a.out, exist_ok=True); json.dump(out, open(f"{a.out}/external_validation.json", "w"), indent=1); print(json.dumps(out, indent=1))
