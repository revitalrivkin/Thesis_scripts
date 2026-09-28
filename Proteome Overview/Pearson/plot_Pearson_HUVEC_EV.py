"""
Pairwise Pearson Correlation — HUVEC-EVs, Normal / HG / HL (filtered dataset)
==============================================================================
Recreated from HUVEC-EVs 01-2026/06-Pearson/pearson_HUVEC-EVs.py, adapted to:
  - Read from the filtered dataset (single-peptide proteins excluded),
    same file used for 4.2.1's PCA and 4.2.3's RA/ES/ΔES analysis — see
    [[17 - MS Proteomics Data - Filtering & Provenance]]
  - House style: Arial throughout, bold title/ticks
  - Interactive: opens a matplotlib window (plt.show()) instead of saving
    headless — resize the window to your preferred dimensions (heatmap
    cells stay square since this is a correlation matrix, but the overall
    figure/margins adjust), then use the window's save icon (floppy disk
    in the toolbar) to export the final PNG/SVG. RUN THIS LOCALLY, not
    headless.

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: (manual, via the toolbar's save button) Pearson_HUVEC_EV.png / .svg
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
# NOTE: no matplotlib.use("Agg") here — this script opens an INTERACTIVE window.
# Run it locally to resize the figure, then use the window's save icon
# (floppy disk in the toolbar) to export.

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
OUT  = os.path.join(BASE, "Proteome Overview", "Pearson")

SAMPLE_COLS = ["Normal-1", "Normal-2", "Normal-3", "HG-1", "HG-2", "HL-1", "HL-2", "HL-3"]

# -- Load filtered EVs data, keep only sample columns --------------------------
df = pd.read_csv(DATA)
df_filtered = df[SAMPLE_COLS]

print(f"Samples: {len(SAMPLE_COLS)}, Proteins (rows): {len(df_filtered)}")

# -- Pairwise Pearson correlation ----------------------------------------------
corr_matrix = df_filtered.corr()
mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

fig, ax = plt.subplots(figsize=(8, 8))
cmap = plt.cm.Greens

# No auto colorbar here — placed manually below so its height matches the
# heatmap's own vertical extent (top under the title, bottom at the x-tick labels)
sns.heatmap(corr_matrix, mask=mask, annot=True, cmap=cmap, ax=ax,
            square=True, linewidths=.5, cbar=False, annot_kws={"size": 12})

plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=12, fontweight="bold", fontfamily=FONT)
plt.setp(ax.get_yticklabels(), rotation=0, fontsize=12, fontweight="bold", fontfamily=FONT)

ax.set_title("Pairwise Correlation Map", fontsize=17, fontweight="bold", fontfamily=FONT, pad=22)
fig.subplots_adjust(right=0.86)

# Manually place the colorbar so it spans exactly the heatmap's y-extent
fig.canvas.draw()
pos = ax.get_position()
cax = fig.add_axes([pos.x1 + 0.03, pos.y0, 0.03, pos.height])
cbar = fig.colorbar(ax.collections[0], cax=cax)

# Round colorbar ticks to clean 0.05 increments spanning the data range
vmin, vmax = corr_matrix.values[~mask].min(), corr_matrix.values[~mask].max()
tick_lo = np.floor(vmin / 0.05) * 0.05
tick_hi = np.ceil(vmax / 0.05) * 0.05
cbar_ticks = np.round(np.arange(tick_lo, tick_hi + 1e-9, 0.05), 2)
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels([f"{t:.2f}" for t in cbar_ticks])
cbar.ax.tick_params(labelsize=10)
for lbl in cbar.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cbar.outline.set_visible(False)

print("Resize the window to your preferred dimensions, then use the toolbar's save icon to export.")
plt.show()
