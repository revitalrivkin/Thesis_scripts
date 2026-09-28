"""
ES scatter — Normal condition, filtered dataset, no annotations
==================================================================
Recreated from 27_ES_scatter_Normal_annotated.py, using the single-peptide-
filtered dataset (see note 17 - MS Proteomics Data - Filtering & Provenance)
and with the explanatory annotation panel / in-plot threshold labels removed
per Revital's request — scatter only.

ES = log2((EV_RA + eps) / (Cell_RA + eps)), eps = 1e-6
Zone reference lines (ES = 0, ±1) and background shading kept, since these
are interpretive gridlines rather than descriptive annotations.

Input:
  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
Output: ES_scatter_Normal.png (+ .svg)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT = "Arial"
plt.rcParams["font.family"] = FONT

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
OUT  = os.path.join(BASE, "Comparative Proteomic Analysis — EVs vs. Parent Cells", "ES figures")

NORMAL_EV_COLS   = ["Normal-1", "Normal-2", "Normal-3"]
NORMAL_CELL_COLS = ["Normal-1", "Normal-2", "Normal-3"]

EPSILON = 1e-6
BLUE    = (0.15, 0.35, 0.85)
GREY    = "#444444"

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
    return df[["Protein.Group", "mean_RA"]]

# -- Load filtered data, compute Normal RA per compartment -------------------
evs   = compute_ra(pd.read_csv(os.path.join(DATA_DIR, "EVs_df_filtered.csv")), NORMAL_EV_COLS)
cells = compute_ra(pd.read_csv(os.path.join(DATA_DIR, "Cells_df_filtered.csv")), NORMAL_CELL_COLS)

# -- ES for proteins with real RA in both compartments ------------------------
evs = evs[evs["mean_RA"].notna()].copy()
evs["ev_pct"] = evs["mean_RA"] * 100

cell_idx = cells.set_index("Protein.Group")["mean_RA"]
cell_ra = evs["Protein.Group"].map(cell_idx)
mask = cell_ra.notna() & (cell_ra > 0)

df = evs[mask].copy()
cell_ra = cell_ra[mask]
df["ES"] = np.log2((df["mean_RA"] + EPSILON) / (cell_ra + EPSILON))
df["log10_ev"] = np.log10(df["ev_pct"].clip(lower=1e-5))
df = df[["log10_ev", "ES"]].dropna()
n = len(df)
print(f"Proteins plotted (detected in both compartments, Normal): {n}")

# -- Axis limits (match original approach) ------------------------------------
x_min = np.percentile(df["log10_ev"], 0.5)
x_max = np.percentile(df["log10_ev"], 99.5)
y_min = np.percentile(df["ES"],       0.2)
y_max = np.percentile(df["ES"],       99.8)
x_pad = (x_max - x_min) * 0.05
y_pad = (y_max - y_min) * 0.05
xlim = (x_min - x_pad, x_max + x_pad)
ylim = (y_min - y_pad, y_max + y_pad)

# -- Figure: scatter only, no annotation panel --------------------------------
fig, ax = plt.subplots(figsize=(8, 6.5))

ax.scatter(df["log10_ev"], df["ES"], color=BLUE, s=10, alpha=0.50,
           linewidths=0, rasterized=True)

ax.axhline(0,  color="#666666", linewidth=1.0, zorder=3)
ax.axhline(1,  color="#444444", linewidth=0.9, linestyle="--", zorder=3)
ax.axhline(-1, color="#444444", linewidth=0.9, linestyle="--", zorder=3)

ax.axhspan(1,       ylim[1], color="#2356A022", zorder=0)
ax.axhspan(-1,       1,      color="#F0F0F022", zorder=0)
ax.axhspan(ylim[0], -1,      color="#A62B1F18", zorder=0)

ax.set_xlim(*xlim)
ax.set_ylim(*ylim)
ax.set_xlabel("Relative Abundance in EVs (log$_{10}$%)", fontsize=11,
              color=GREY, fontweight="bold", fontfamily=FONT)
ax.set_ylabel("Enrichment Score (ES)", fontsize=11, color=GREY,
              fontweight="bold", fontfamily=FONT)
ax.tick_params(labelsize=10)
for lbl in ax.get_xticklabels() + ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color("#bbbbbb")

fig.suptitle(f"EV Enrichment Score vs. Relative Abundance — Normal ({n:,} proteins)",
             fontsize=15, fontweight="bold", color=GREY, y=0.98, fontfamily=FONT)

plt.tight_layout()
for ext in ("png", "svg"):
    path = os.path.join(OUT, f"ES_scatter_Normal.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("Saved:", path)
