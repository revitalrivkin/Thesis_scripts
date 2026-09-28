"""
Top 10 most abundant EV proteins by RA — HG and HL, side by side
===================================================================
Two triptychs (Cells | protein names | EVs) side by side: HG (left half),
HL (right half). Same style/logic as the Normal figure — ranked by EV RA,
parent-cell RA shown for context, equal 0-15 scale both sides per condition.

Note: EVs has only 2 replicates for HG (HG-3 excluded, known outlier);
Cells has all 3 HG replicates.

Input:
  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
Output: Top10_RA_HG_HL_EVs_vs_Cells.png (+ .svg)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FONT = "Arial"
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
OUT  = os.path.join(BASE, "Comparative Proteomic Analysis — EVs vs. Parent Cells", "RA figures")

EV_COLS_BY_COND   = {"HG": ["HG-1", "HG-2"],            "HL": ["HL-1", "HL-2", "HL-3"]}
CELL_COLS_BY_COND = {"HG": ["HG-1", "HG-2", "HG-3"],    "HL": ["HL-1", "HL-2", "HL-3"]}

CELL_COLOR = (0, 0, 1)
EV_COLOR   = (0.53, 0.71, 1.0)
GREY       = "#444444"
FS_TITLE, FS_PANEL, FS_AXIS, FS_PROT, FS_VAL, FS_TICK = 15, 13, 12, 12, 10, 10
AXIS_MAX = None  # set after computing data, must fit the true max across both conditions

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
    return df

evs_raw   = pd.read_csv(os.path.join(DATA_DIR, "EVs_df_filtered.csv"))
cells_raw = pd.read_csv(os.path.join(DATA_DIR, "Cells_df_filtered.csv"))

def get_top10(cond):
    evs   = compute_ra(evs_raw, EV_COLS_BY_COND[cond])
    cells = compute_ra(cells_raw, CELL_COLS_BY_COND[cond])

    top10 = evs.nlargest(10, "mean_RA")[["Protein.Group", "Genes", "mean_RA"]].copy()
    top10["ev_pct"] = top10["mean_RA"] * 100
    cell_lkp = cells.set_index("Protein.Group")["mean_RA"]
    top10["cell_pct"] = top10["Protein.Group"].map(cell_lkp).fillna(0) * 100
    top10["label"] = (
        top10["Genes"].astype(str).str.split(";").str[0]
        .str.replace(r"^cRAP-", "", regex=True).str.strip()
    )
    top10 = top10.sort_values("ev_pct", ascending=True).reset_index(drop=True)
    print(f"\nTop 10 EV proteins by RA ({cond}):")
    print(top10[["label", "ev_pct", "cell_pct"]].to_string(index=False))
    return top10

data = {cond: get_top10(cond) for cond in ["HG", "HL"]}

# Axis must fit the true max across both conditions (avoid silently clipping bars)
max_val = max(top10["ev_pct"].max() for top10 in data.values())
AXIS_MAX = int(np.ceil((max_val * 1.1) / 2) * 2)  # round up to next even number, +10% padding
print(f"\nMax EV RA value: {max_val:.2f}% -> axis max set to {AXIS_MAX}")

# -- Figure: 2 triptychs side by side -----------------------------------------
fig = plt.figure(figsize=(20, 6))
outer = gridspec.GridSpec(1, 2, wspace=0.12, left=0.02, right=0.98, top=0.82, bottom=0.13)

for panel_idx, cond in enumerate(["HG", "HL"]):
    top10 = data[cond]
    labels, ev_vals, cell_vals = top10["label"].tolist(), top10["ev_pct"].tolist(), top10["cell_pct"].tolist()
    y = np.arange(len(labels))

    inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[panel_idx],
                                              width_ratios=[4.5, 1.6, 4.5], wspace=0.0)
    ax_c = fig.add_subplot(inner[0]); ax_l = fig.add_subplot(inner[1]); ax_e = fig.add_subplot(inner[2])
    BAR_H = 0.58

    ax_l.set_xlim(0, 1); ax_l.set_ylim(-0.5, len(y) - 0.5)
    ax_l.set_yticks([]); ax_l.set_xticks([])
    for spine in ax_l.spines.values():
        spine.set_visible(False)
    for i, name in enumerate(labels):
        ax_l.text(0.5, i, name, va="center", ha="center", fontsize=FS_PROT,
                  color=GREY, fontweight="bold", fontfamily=FONT)

    ax_c.barh(y, cell_vals, height=BAR_H, color=CELL_COLOR, align="center")
    ax_c.invert_xaxis()
    ax_c.set_ylim(-0.5, len(y) - 0.5); ax_c.set_yticks([])
    ax_c.set_title(f"Cells ({cond})", fontsize=FS_PANEL, fontweight="bold", color=CELL_COLOR, pad=6, fontfamily=FONT)
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.tick_params(axis="x", labelsize=FS_TICK)
    for i, v in enumerate(cell_vals):
        if v >= 0.001:
            ax_c.text(v, i, f" {v:.2f}", va="center", ha="right", fontsize=FS_VAL,
                      color=GREY, fontfamily=FONT, fontweight="bold")

    ax_e.barh(y, ev_vals, height=BAR_H, color=EV_COLOR, align="center")
    ax_e.set_ylim(-0.5, len(y) - 0.5); ax_e.set_yticks([])
    ax_e.set_title(f"EVs ({cond})", fontsize=FS_PANEL, fontweight="bold", color=EV_COLOR, pad=6, fontfamily=FONT)
    ax_e.spines[["top", "left", "right"]].set_visible(False)
    ax_e.tick_params(axis="x", labelsize=FS_TICK)
    for i, v in enumerate(ev_vals):
        ax_e.text(v, i, f"  {v:.1f}", va="center", ha="left", fontsize=FS_VAL,
                  color=GREY, fontfamily=FONT, fontweight="bold")

    ticks = np.arange(0, AXIS_MAX + 1, 2)
    ax_e.set_xlim(0, AXIS_MAX); ax_c.set_xlim(AXIS_MAX, 0)
    ax_e.set_xticks(ticks); ax_c.set_xticks(ticks)
    for ax in (ax_c, ax_e):
        ax.axvline(0, color="#bbbbbb", linewidth=0.8, zorder=0)
        for lbl in ax.get_xticklabels():
            lbl.set_fontfamily(FONT)

fig.suptitle(
    "Top 10 Most Abundant EV Proteins — HG and HL, with Parent-Cell Abundance",
    fontsize=FS_TITLE, fontweight="bold", y=0.95, color=GREY, fontfamily=FONT
)
fig.text(0.5, 0.03, "Relative abundance (%)", ha="center", va="center",
          fontsize=FS_AXIS, fontweight="bold", color=GREY, fontfamily=FONT)

for ext in ("png", "svg"):
    path = os.path.join(OUT, f"Top10_RA_HG_HL_EVs_vs_Cells.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("\nSaved:", path)
