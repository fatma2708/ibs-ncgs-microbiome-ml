"""Jeu de données synthétique minimal (CI) : mêmes formats que les vraies données."""
import os, sys, numpy as np, pandas as pd
out = sys.argv[1]; rng = np.random.RandomState(0)
os.makedirs(f"{out}/bracken_out", exist_ok=True)
rows, taxa = [], [f"Species {i}" for i in range(60)]
for g in ("NCGS", "IBS"):
    for k in range(1, 9):
        run = f"SRRTOY{g}{k}"; mg = 10 * k + (5 if g == "IBS" else 0)
        rows.append((len(rows) + 1, run, f"Mg_{mg}", g))
        p = rng.dirichlet(np.ones(60) * 0.3)
        if g == "IBS": p[:5] *= 3
        p /= p.sum()
        pd.DataFrame({"name": taxa + ["Homo sapiens"], "fraction_total_reads": list(p * .99) + [.01]}).to_csv(f"{out}/bracken_out/{run}.bracken", sep="\t", index=False)
pd.DataFrame(rows, columns=["SampleN", "Run", "Mg_name", "Group"]).to_csv(f"{out}/sample_id_map.csv", index=False)
L = [(f"Met{m}", f"{g}_{k}", float(rng.gamma(2) * (1.5 if (g == "IBS" and m < 4) else 1)))
     for g in ("IBS", "NCGS") for k in range(1, 9) for m in range(30)]
with pd.ExcelWriter(f"{out}/s3.xlsx") as w:
    pd.DataFrame(L, columns=["Metabolite", "Sample_ID", "Response Ratio"]).to_excel(w, sheet_name="LCMS", index=False)
