"""
4.2.3.5 Enriched/Depleted by Pathology — Pathology-Depleted Cargo (strict criteria)
================================================================================
Reuses the exact ES/ΔES computation from the Constitutive Cargo script
(../GO enrichment figures/Constitutive cargo/01_build_constitutive_cargo_list.py)
— same RA formula, same epsilon, same filtered dataset — applying the
pathology-depleted gate instead of the constitutive gate:

  ES_Normal > 1  AND  ΔES_<condition> < -1

("Actively enriched in Normal EV cargo, whose packaging then dropped
substantially under HG or HL" — the strict criteria; a looser version using
ES_Normal > -1 instead of > 1 was floated but not yet run.)

HG and HL evaluated SEPARATELY (Revital's call, 2026-08-23) — not requiring
both simultaneously as a *criterion*. Note the ES/dES computation itself
still requires each protein to be detected in Normal, HG, AND HL (3-way
inner join) before either the HG or HL gate is evaluated — this script
draws on all three conditions at once, so per Revital's filing rule
(2026-08-25) it stays at the top level of the folder rather than being
split into HG/HL subfolders; only its outputs are filed per-condition.

Input:
  ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
Output: Depleted/HG/01_pathology_depleted_cargo_HG.xlsx
        Depleted/HL/01_pathology_depleted_cargo_HL.xlsx
"""

import os
import numpy as np
import pandas as pd

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_BY_COND = {"HG": os.path.join(HERE, "Depleted", "HG"),
               "HL": os.path.join(HERE, "Depleted", "HL")}

EV_COLS_BY_COND   = {"Normal": ["Normal-1", "Normal-2", "Normal-3"],
                     "HG": ["HG-1", "HG-2"], "HL": ["HL-1", "HL-2", "HL-3"]}
CELL_COLS_BY_COND = {"Normal": ["Normal-1", "Normal-2", "Normal-3"],
                     "HG": ["HG-1", "HG-2", "HG-3"], "HL": ["HL-1", "HL-2", "HL-3"]}
EPSILON = 1e-6

def log2_to_raw(series):
    s = series.copy()
    s[s <= 13] = np.nan
    return np.power(2, s)

def compute_ra(df, cols):
    df = df.copy()
    df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")
    raw = pd.DataFrame(index=df.index)
    for col in cols:
        raw[col] = log2_to_raw(df[col])
    ra = pd.DataFrame(index=df.index)
    for col in cols:
        col_sum = raw[col].sum(skipna=True)
        ra[col + "_RA"] = raw[col] / col_sum
    df["mean_RA"] = ra[[c + "_RA" for c in cols]].mean(axis=1)
    return df[["Protein.Group", "Genes", "mean_RA"]].rename(columns={"mean_RA": "RA"})

evs_raw   = pd.read_csv(os.path.join(DATA_DIR, "EVs_df_filtered.csv"))
cells_raw = pd.read_csv(os.path.join(DATA_DIR, "Cells_df_filtered.csv"))

ev_ra   = {c: compute_ra(evs_raw, EV_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}
cell_ra = {c: compute_ra(cells_raw, CELL_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}

es = {}
for c in ["Normal", "HG", "HL"]:
    m = ev_ra[c].merge(cell_ra[c][["Protein.Group", "RA"]], on="Protein.Group",
                        how="inner", suffixes=("_ev", "_cell"))
    m = m[m["RA_ev"].notna() & m["RA_cell"].notna() & (m["RA_cell"] > 0)].copy()
    m["ES"] = np.log2((m["RA_ev"] + EPSILON) / (m["RA_cell"] + EPSILON))
    es[c] = m[["Protein.Group", "Genes", "ES"]].rename(columns={"ES": f"ES_{c}"})
    print(f"ES_{c}: {len(es[c])} proteins (detected in both compartments)")

df = es["Normal"].merge(es["HG"][["Protein.Group", "ES_HG"]], on="Protein.Group", how="inner")
df = df.merge(es["HL"][["Protein.Group", "ES_HL"]], on="Protein.Group", how="inner")
print(f"\nProteins with ES in all 3 conditions: {len(df)}")

df["dES_HG"] = df["ES_HG"] - df["ES_Normal"]
df["dES_HL"] = df["ES_HL"] - df["ES_Normal"]
df["Gene"] = (
    df["Genes"].astype(str).str.split(";").str[0]
    .str.replace(r"^cRAP-", "", regex=True).str.strip()
)

for cond in ["HG", "HL"]:
    mask = (df["ES_Normal"] > 1) & (df[f"dES_{cond}"] < -1)
    depleted = df[mask].copy()
    print(f"\nPathology-depleted cargo, {cond} (ES_Normal>1, dES_{cond}<-1): {len(depleted)} proteins")

    out_df = depleted[["Protein.Group", "Gene", "ES_Normal", f"ES_{cond}", f"dES_{cond}"]]
    out_df = out_df.sort_values(f"dES_{cond}")
    out_path = os.path.join(OUT_BY_COND[cond], f"01_pathology_depleted_cargo_{cond}.xlsx")
    out_df.to_excel(out_path, index=False)
    print("Saved:", out_path)
