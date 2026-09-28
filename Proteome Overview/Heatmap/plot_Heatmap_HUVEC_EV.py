"""
Whole-Proteome Heatmap — HUVEC-EVs, Normal / HG / HL (filtered dataset)
==============================================================================
Recreated from HUVEC-EVs 01-2026/04-proteome_heatmap/Proteome_heatmap.py, adapted to:
  - Read from the filtered dataset (single-peptide proteins excluded),
    same file used for 4.2.1's PCA/Pearson and 4.2.3's RA/ES/ΔES analysis —
    see [[17 - MS Proteomics Data - Filtering & Provenance]]
  - Row-clustered only (col_cluster=False), no dendrogram shown (matches
    the reference script's show_dendrogram=False configuration)
  - Colorbar merged into one combined figure instead of a separate window,
    manually positioned (same approach as the Pearson figure), with no
    outline
  - Interactive: opens a matplotlib window (plt.show()) instead of saving
    headless — drag/resize the window to whatever aspect ratio (square or
    rectangular) you want, then use the window's save icon (floppy disk in
    the toolbar) to export the final PNG/SVG. RUN THIS LOCALLY, not headless.

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: (manual, via the toolbar's save button) Heatmap_HUVEC_EV.png / .svg
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
import matplotlib.pyplot as plt
# NOTE: no matplotlib.use("Agg") here — this script opens an INTERACTIVE window.
# Run it locally to resize the figure to your preferred aspect ratio, then use
# the window's save icon (floppy disk in the toolbar) to export.

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
OUT  = os.path.join(BASE, "Proteome Overview", "Heatmap")

GENE_COL   = "Genes"
LEVEL_COLS = ["Normal-1", "Normal-2", "Normal-3", "HG-1", "HG-2", "HL-1", "HL-2", "HL-3"]
TITLE      = "Whole-Proteome Heatmap of HUVEC-EVs"
VMIN, VMAX = -3, 3

# -- Load filtered EVs data -----------------------------------------------------
df = pd.read_csv(DATA)
df[LEVEL_COLS] = df[LEVEL_COLS].apply(pd.to_numeric, errors="coerce")
df_data = df.set_index(GENE_COL)[LEVEL_COLS].dropna()

print(f"Proteins: {len(df_data)}, Samples: {len(LEVEL_COLS)}")

# -- Row z-score + hierarchical row clustering (no column clustering) ----------
df_zscore = df_data.apply(zscore, axis=1)
df_zscore = df_zscore.fillna(0)

row_linkage = linkage(df_zscore.values, method="average", metric="euclidean")
row_order = leaves_list(row_linkage)
df_plot = df_zscore.iloc[row_order]

# -- Figure: heatmap + manually placed colorbar (no outline) -------------------
fig, ax = plt.subplots(figsize=(14, 12))
sns.heatmap(df_plot, ax=ax, cmap="coolwarm", xticklabels=True, yticklabels=False,
            cbar=False, vmin=VMIN, vmax=VMAX)

ax.set_ylabel("")
plt.setp(ax.get_xticklabels(), rotation=0, fontsize=14, fontweight="bold", fontfamily=FONT)
ax.set_title(TITLE, fontsize=22, fontweight="bold", fontfamily=FONT, pad=22)

fig.subplots_adjust(right=0.90)
fig.canvas.draw()
pos = ax.get_position()
cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.02, pos.height])
sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=VMIN, vmax=VMAX))
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Z-Score", fontsize=14, fontweight="bold", fontfamily=FONT, labelpad=10)
cbar.ax.tick_params(labelsize=11)
for lbl in cbar.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cbar.outline.set_visible(False)

print("Resize the window to your preferred aspect ratio, then use the toolbar's save icon to export.")
plt.show()
