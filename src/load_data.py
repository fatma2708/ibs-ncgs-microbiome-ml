"""Chargement commun (mapping hypothétique k-ième Mg_ trié -> IBS_k/NCGS_k)."""
import glob, os, numpy as np, pandas as pd
def load(bracken, mapfile, s3):
    m = pd.read_csv(mapfile); m = m[m.Group.isin(["IBS", "NCGS"])].copy()
    m["mg"] = m.Mg_name.str.extract(r"(\d+)").astype(int)
    m = m.sort_values(["Group", "mg"]); m["k"] = m.groupby("Group").cumcount() + 1
    m["met_id"] = m.Group + "_" + m.k.astype(str)
    runs = {os.path.basename(f)[:-8]: f for f in glob.glob(f"{bracken}/*.bracken")}
    m = m[m.Run.isin(runs)].copy()
    tax, depth = {}, {}
    for r in m.Run:
        d = pd.read_csv(runs[r], sep="\t"); tax[r] = d.set_index("name")["fraction_total_reads"]
        if "new_est_reads" in d: depth[r] = float(d.new_est_reads.sum())
    T = pd.DataFrame(tax).fillna(0).T.drop(columns=["Homo sapiens"], errors="ignore"); T = T.div(T.sum(axis=1), axis=0)
    L = pd.read_excel(s3, sheet_name="LCMS"); L["v"] = pd.to_numeric(L["Response Ratio"], errors="coerce")
    M = L.pivot_table(index="Sample_ID", columns="Metabolite", values="v").fillna(0)
    m = m[m.met_id.isin(M.index)]; T = T.loc[m.Run]; M = M.loc[m.met_id]; M.index = m.Run.values
    y = (m.Group == "IBS").astype(int).values
    dp = pd.Series(depth).reindex(m.Run).values if depth else None
    return T, M, y, dp
