"""Figures du README à partir de results/ (palette daltonien-compatible)."""
import argparse, json, os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.decomposition import PCA
from common import prep_tax
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--res", default="results"); ap.add_argument("--out", default="figures")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
C = {"IBS": "#0072B2", "NCGS": "#D55E00", "tax": "#0072B2", "met": "#009E73", "fusion": "#CC79A7"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
NAME = {"tax": "Taxonomie", "met": "Métabolomique (LC-MS)", "fusion": "Fusion"}
ex = json.load(open(f"{a.res}/extra.json")); oof = pd.read_csv(f"{a.res}/oof_predictions.csv")

# 1 ROC
fig, ax = plt.subplots(figsize=(4.6, 4.4))
for k in ("tax", "met", "fusion"):
    f, t, _ = roc_curve(oof.y, oof[k]); lo, hi = ex["bootstrap_ci"][k]["ci95"]
    ax.plot(f, t, color=C[k], lw=2, label=f'{NAME[k]}  AUC {ex["bootstrap_ci"][k]["auc"]:.2f} [{lo:.2f}–{hi:.2f}]')
ax.plot([0, 1], [0, 1], "--", color="grey", lw=1); ax.set_xlabel("1 − spécificité"); ax.set_ylabel("Sensibilité")
ax.set_title(f"IBS vs NCGS (n = {len(oof)}, prédictions hors-pli)", fontsize=10); ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(f"{a.out}/roc.png", dpi=200); plt.close()

# 2 courbe d'apprentissage
lc = pd.read_csv(f"{a.res}/learning_curve.csv", header=[0, 1], index_col=0)
fig, ax = plt.subplots(figsize=(4.6, 3.8))
for k in ("tax", "met", "fusion"):
    m, s = lc[(k, "mean")], lc[(k, "std")]; ax.plot(lc.index, m, "-o", color=C[k], label=NAME[k]); ax.fill_between(lc.index, m - s, m + s, color=C[k], alpha=.15)
ax.axhline(.5, ls="--", color="grey", lw=1); ax.set_xlabel("Nombre d'échantillons"); ax.set_ylabel("AUC (validation croisée imbriquée)")
ax.set_title("Courbe d'apprentissage", fontsize=10); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(f"{a.out}/learning_curve.png", dpi=200); plt.close()

# 3 SHAP (descriptif)
fig, axs = plt.subplots(1, 2, figsize=(9, 4.2))
for ax, nm, ttl in zip(axs, ("taxonomy", "metabolomics"), ("Taxonomie (top 12)", "Métabolites LC-MS (top 12)")):
    d = pd.read_csv(f"{a.res}/shap_{nm}.csv").head(12).iloc[::-1]
    ax.barh(d.feature, d.mean_abs_shap, color=[C["IBS"] if v > 0 else C["NCGS"] for v in d.coef_IBS])
    ax.set_title(ttl, fontsize=10); ax.set_xlabel("|SHAP| moyen"); ax.tick_params(axis="y", labelsize=8)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.text(0.5, 0.01, "bleu : associé à IBS   orange : associé à NCGS (signe du coefficient) — ajustement sur tous les échantillons, descriptif", ha="center", fontsize=8)
none=None; fig.savefig(f"{a.out}/shap.png", dpi=200); plt.close()

# 4 PCA des profils CLR
T, M, y, _ = load(a.bracken, a.map, a.s3); X, _, _ = prep_tax(T, T)
Z = PCA(2).fit(X); S = Z.transform(X); v = Z.explained_variance_ratio_
fig, ax = plt.subplots(figsize=(4.6, 4))
for g, lab in ((1, "IBS"), (0, "NCGS")): ax.scatter(S[y == g, 0], S[y == g, 1], s=45, color=C[lab], label=lab, edgecolor="white")
ax.set_xlabel(f"PC1 ({v[0]:.0%})"); ax.set_ylabel(f"PC2 ({v[1]:.0%})"); ax.set_title("PCA des profils taxonomiques (CLR)", fontsize=10)
ax.legend(frameon=False); fig.tight_layout(); fig.savefig(f"{a.out}/pca.png", dpi=200); plt.close(); print("ok")
