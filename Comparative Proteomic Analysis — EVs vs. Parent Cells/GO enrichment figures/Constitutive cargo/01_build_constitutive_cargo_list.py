"""
Build constitutive cargo protein list — filtered dataset, finalized threshold
================================================================================
Computes ES for Normal/HG/HL and ΔES (HG, HL vs Normal) from the single-
peptide-filtered dataset, then applies the finalized 2-gate definition:

  ES_Normal > 1  AND  |ΔES_HG| < 1  AND  |ΔES_HL| < 1

(Finalized 2026-08-13 — supersedes both older note-06 versions: 560 proteins
under ES>1-all-conditions, and 199 proteins under the old |ΔES|<0.5 threshold.)

Input:
  ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
Output: 01_constitutive_cargo.xlsx (Protein.Group, Gene, ES_Normal, dES_HG, dES_HL)
"""

import os
import numpy as np
import pandas as pd

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
OUT = os.path.join(BASE, "Comparative Proteomic Analysis — EVs vs. Parent Cells",
                    "GO enrichment figures", "Constitutive cargo")

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

# -- RA per condition, per compartment ---------------------------------------
ev_ra   = {c: compute_ra(evs_raw, EV_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}
cell_ra = {c: compute_ra(cells_raw, CELL_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}

# -- ES per condition (merge EV and Cell RA on Protein.Group) -----------------
es = {}
for c in ["Normal", "HG", "HL"]:
    m = ev_ra[c].merge(cell_ra[c][["Protein.Group", "RA"]], on="Protein.Group",
                        how="inner", suffixes=("_ev", "_cell"))
    m = m[m["RA_ev"].notna() & m["RA_cell"].notna() & (m["RA_cell"] > 0)].copy()
    m["ES"] = np.log2((m["RA_ev"] + EPSILON) / (m["RA_cell"] + EPSILON))
    es[c] = m[["Protein.Group", "Genes", "ES"]].rename(columns={"ES": f"ES_{c}"})
    print(f"ES_{c}: {len(es[c])} proteins (detected in both compartments)")

# -- Merge all three conditions on Protein.Group (inner — need all 3) --------
df = es["Normal"].merge(es["HG"][["Protein.Group", "ES_HG"]], on="Protein.Group", how="inner")
df = df.merge(es["HL"][["Protein.Group", "ES_HL"]], on="Protein.Group", how="inner")
print(f"\nProteins with ES in all 3 conditions: {len(df)}")

df["dES_HG"] = df["ES_HG"] - df["ES_Normal"]
df["dES_HL"] = df["ES_HL"] - df["ES_Normal"]

# -- Apply constitutive cargo gate --------------------------------------------
mask = (df["ES_Normal"] > 1) & (df["dES_HG"].abs() < 1) & (df["dES_HL"].abs() < 1)
constitutive = df[mask].copy()
print(f"Constitutive cargo proteins (ES_Normal>1, |dES_HG|<1, |dES_HL|<1): {len(constitutive)}")

constitutive["Gene"] = (
    constitutive["Genes"].astype(str).str.split(";").str[0]
    .str.replace(r"^cRAP-", "", regex=True).str.strip()
)

out_df = constitutive[["Protein.Group", "Gene", "ES_Normal", "dES_HG", "dES_HL"]]
out_path = os.path.join(OUT, "01_constitutive_cargo.xlsx")
out_df.to_excel(out_path, index=False)
print("\nSaved:", out_path)
