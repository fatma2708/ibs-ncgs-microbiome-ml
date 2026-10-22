"""Freeze the taxonomy model on all samples of the discovery cohort (n=28). Nothing is tuned afterwards.
Usage: python freeze_model.py --bracken B --map M --s3 S --out models/frozen_model_v1.joblib"""
import argparse, json, hashlib, os, joblib, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from common import prep_tax, fit_model
from load_data import load
ap = argparse.ArgumentParser()
for k in ("bracken", "map", "s3"): ap.add_argument(f"--{k}", required=True)
ap.add_argument("--out", default="models/frozen_model_v1.joblib"); a = ap.parse_args()
T, M, y, _ = load(a.bracken, a.map, a.s3)
A, _, cols = prep_tax(T, T)
sc, clf, C = fit_model(A, y)
model = dict(features=list(cols), scaler=sc, clf=clf, C=C, pseudocount=1e-4, prevalence=0.2, positive_class="IBS", n_train=int(len(y)),
             pipeline="fastp > bowtie2(GRCh38) > kraken2 k2_standard_08gb_20240605 > bracken (species, read length 150)")
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True); joblib.dump(model, a.out)
h = hashlib.sha256(open(a.out, "rb").read()).hexdigest()[:16]
json.dump(dict(sha256_16=h, n_features=len(cols), C=C, n_train=int(len(y))), open(a.out.replace(".joblib", ".json"), "w"), indent=1)
print(f"frozen: {len(cols)} species, C={C}, sha256[:16]={h}")
