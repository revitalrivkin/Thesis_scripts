"""
Delta Enrichment Score (dES) scatter - HG and HL vs Normal, rebuilt on the canonical filtered dataset
==============================================================================
REBUILD (2026-09-27) of dES_scatter_HG_HL.py, which read a stale pipeline
(Data analysis projects\05-2026_Comparative_Analysis\02-normalization\
04_EVs_relative_abundance.xlsx / 03_cells_relative_abundance.xlsx, dated
2026-04-30 - a leftover from before this whole workstream moved onto the
single-peptide-filtered dataset the same day as the sibling RA/ES scripts,
2026-08-15). That script also never applied its own "reliable" filter to the
cells file. Rebuilt here to read the same filtered CSVs as
plot_top10_RA_Normal_EVs.py / plot_ES_scatter_Normal.py, with the same
RA/ES computation (log2_to_raw, compute_ra, pseudocount c = 1e-6), so all of
4.2.3.1/4.2.3.2/4.2.3.3 are built from one consistent dataset.

dES = ES_condition - ES_Normal, where ES = log2[(EV_RA + c) / (Cell_RA + c)].
A protein is included only if it has real (non-floor) RA in EVs and in cells,
for BOTH the condition and Normal (4 quantities per protein).

Input:
  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
Output: dES_scatter_HG_HL_rebuilt.png (+ .svg)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FONT = "Arial"
plt.rcParams["font.family"] = FONT

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
OUT = os.path.dirname(os.path.abspath(__file__))

EV_COLS = {"Normal": ["Normal-1", "Normal-2", "Normal-3"], "HG": ["HG-1", "HG-2"], "HL": ["HL-1", "HL-2", "HL-3"]}
CELL_COLS = {"Normal": ["Normal-1", "Normal-2", "Normal-3"], "HG": ["HG-1", "HG-2", "HG-3"], "HL": ["HL-1", "HL-2", "HL-3"]}

EPSILON = 1e-6
GREY = "#444444"
CONDITIONS = [("HG", (0.99, 0.57, 0.00)), ("HL", (0.95, 0.15, 0.05))]
COL_GAINED, COL_STABLE, COL_LOST = "#DAE8F5", "#F2F2F2", "#FADADD"


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


evs_raw = pd.read_csv(os.path.join(DATA_DIR, "EVs_df_filtered.csv"))
cells_raw = pd.read_csv(os.path.join(DATA_DIR, "Cells_df_filtered.csv"))

ev_ra = {c: compute_ra(evs_raw, EV_COLS[c]).set_index("Protein.Group")["mean_RA"] for c in EV_COLS}
cell_ra = {c: compute_ra(cells_raw, CELL_COLS[c]).set_index("Protein.Group")["mean_RA"] for c in CELL_COLS}


def get_dES(cond):
    idx = ev_ra[cond].dropna().index
    df = pd.DataFrame(index=idx)
    df["ev_cond"] = ev_ra[cond]
    df["cell_cond"] = df.index.map(cell_ra[cond])
    df["ev_norm"] = df.index.map(ev_ra["Normal"])
    df["cell_norm"] = df.index.map(cell_ra["Normal"])
    df = df.dropna()  # real RA required in all 4 quantities (EV/cell x condition/Normal)
    es_cond = np.log2((df["ev_cond"] + EPSILON) / (df["cell_cond"] + EPSILON))
    es_norm = np.log2((df["ev_norm"] + EPSILON) / (df["cell_norm"] + EPSILON))
    out = pd.DataFrame({"log10_ev": np.log10((df["ev_cond"] * 100).clip(lower=1e-5)), "dES": es_cond - es_norm})
    n = len(out)
    gained, lost = int((out["dES"] > 1).sum()), int((out["dES"] < -1).sum())
    stable = n - gained - lost
    print(f"{cond} vs Normal: n={n}  gained={gained} ({100*gained/n:.1f}%)  "
          f"stable={stable} ({100*stable/n:.1f}%)  lost={lost} ({100*lost/n:.1f}%)")
    return out


all_data = {cond: get_dES(cond) for cond, _ in CONDITIONS}

all_x = np.concatenate([d["log10_ev"].values for d in all_data.values()])
all_y = np.concatenate([d["dES"].values for d in all_data.values()])
x_min, x_max = np.percentile(all_x, 0.5), np.percentile(all_x, 99.5)
y_min, y_max = np.percentile(all_y, 0.2), np.percentile(all_y, 99.8)
x_pad, y_pad = (x_max - x_min) * 0.05, (y_max - y_min) * 0.05
xlim, ylim = (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad)

fig = plt.figure(figsize=(11, 5.5))
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.10, left=0.09, right=0.97, top=0.86, bottom=0.14)

for col_i, (cond, color) in enumerate(CONDITIONS):
    ax = fig.add_subplot(gs[col_i])
    df = all_data[cond]
    n = len(df)

    ax.axhspan(1, ylim[1], color=COL_GAINED, zorder=0)
    ax.axhspan(-1, 1, color=COL_STABLE, zorder=0)
    ax.axhspan(ylim[0], -1, color=COL_LOST, zorder=0)

    ax.scatter(df["log10_ev"], df["dES"], color=color, s=10, alpha=0.50, linewidths=0, rasterized=True, zorder=2)

    ax.axhline(0, color="#666666", linewidth=1.0, zorder=3)
    ax.axhline(1, color="#888888", linewidth=0.9, linestyle="--", zorder=3)
    ax.axhline(-1, color="#888888", linewidth=0.9, linestyle="--", zorder=3)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_title(f"{cond} vs Normal\n({n:,} proteins)", fontsize=13, fontweight="bold", color=color, pad=6,
                 fontfamily=FONT)
    ax.set_xlabel("Relative Abundance in EVs (log$_{10}$%)", fontsize=11, color=GREY, fontfamily=FONT)

    if col_i == 0:
        ax.set_ylabel("\u0394ES  (vs Normal)", fontsize=11, color=GREY, fontfamily=FONT)
    else:
        ax.set_yticklabels([])

    ax.tick_params(labelsize=10)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#bbbbbb")

fig.suptitle("Change in EV Packaging Preference - Diabetic vs Normal", fontsize=13, fontweight="bold", color=GREY,
             y=0.99, fontfamily=FONT)

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"dES_scatter_HG_HL_rebuilt.{ext}"), dpi=300, bbox_inches="tight")
plt.close()
print("Saved: dES_scatter_HG_HL_rebuilt.png / .svg")
