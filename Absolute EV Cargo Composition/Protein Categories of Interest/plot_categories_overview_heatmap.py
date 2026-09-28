"""
Categories of Interest — Log2FC Heatmap (all detected proteins, non-DEP scope)
==============================================================================
Overview heatmap of every detected protein belonging to a category of
interest (no DEP filtering — that visualization question is deferred).
Rows grouped into category blocks (fixed order); within each block, rows
are hierarchically clustered (average-linkage, Euclidean distance — same
method as the 4.2.1.3 whole-proteome heatmap) on the two displayed columns.
No row labels (too many proteins to label readably at this scope).

Columns: log2 fold change vs. Normal, for HG and HL (Student's t-test
Difference columns — group-mean-level log2FC, same values used to define
DEPs in the volcano, but shown here unfiltered by significance).

Category assignment: same GO-derived gene sets and priority-resolution
order as the scatter/PPI figures (Absolute EV Cargo Composition/Volcano/
Volcano with categories/category_genesets.csv) — a gene belonging to
multiple categories is assigned to the most specific one.

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
        ../Volcano/Volcano with categories/category_genesets.csv
Output: Categories_overview_heatmap.png (+ .svg)
"""

import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "category_genesets.csv")

# Fixed category order + colors (same as the scatter/PPI figures)
CATEGORY_ORDER = [
    "Adhesion / Docking / Uptake",
    "EV Biogenesis",
    "Immune-evasion",
    "ECM Organization",
    "Cardioprotective",
]
CATEGORY_COLORS = {
    "Adhesion / Docking / Uptake": "#1f77b4",
    "EV Biogenesis": "#9467bd",
    "Immune-evasion": "#2ca02c",
    "ECM Organization": "#8c564b",
    "Cardioprotective": "#d62728",
}
# Priority order for resolving multi-category genes (most specific first)
CATEGORY_PRIORITY = ["EV Biogenesis", "ECM Organization", "Immune-evasion", "Cardioprotective", "Adhesion / Docking / Uptake"]

# -- Load data --------------------------------------------------------------------
df = pd.read_csv(DATA)
d = df[["Genes", "Student's T-test Difference HG_Normal", "Student's T-test Difference HL_Normal"]].copy()
d = d.rename(columns={"Student's T-test Difference HG_Normal": "HG",
                       "Student's T-test Difference HL_Normal": "HL"}).dropna()
d = d.set_index("Genes")
print(f"Detected proteins with both log2FC values: {len(d)}")

genesets = pd.read_csv(GENESETS)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()

def assign_category(gene):
    cats = gene_to_cats.get(gene, set())
    for cat in CATEGORY_PRIORITY:
        if cat in cats:
            return cat
    return None

d["Category"] = [assign_category(g) for g in d.index]
d = d[d["Category"].notna()]
print("Rows per category:")
print(d["Category"].value_counts())

# -- Cluster rows within each category block (Euclidean, average-linkage) ------
ordered_blocks = []
for cat in CATEGORY_ORDER:
    sub = d[d["Category"] == cat][["HG", "HL"]]
    if len(sub) < 3:
        ordered_blocks.append(sub)
        continue
    Z = linkage(sub.values, method="average", metric="euclidean")
    order = leaves_list(Z)
    ordered_blocks.append(sub.iloc[order])

plot_df = pd.concat(ordered_blocks)
block_sizes = [len(b) for b in ordered_blocks]
print("Final block sizes (plot order):", dict(zip(CATEGORY_ORDER, block_sizes)))

# -- Plot ---------------------------------------------------------------------------
n = len(plot_df)
fig_h = max(8, n * 0.018 + 2)
fig, ax = plt.subplots(figsize=(6, fig_h))

vmax = np.percentile(np.abs(plot_df[["HG", "HL"]].values), 99)
vmin = -vmax
im = ax.imshow(plot_df[["HG", "HL"]].values, aspect="auto", cmap="coolwarm", vmin=vmin, vmax=vmax,
                interpolation="none")

ax.set_xticks([0, 1])
ax.set_xticklabels(["HG", "HL"], fontsize=13, fontweight="bold", fontfamily=FONT)
ax.set_yticks([])
ax.set_xlabel("")

# Category sidebar (colored bar to the left of the heatmap) + block labels
sidebar_x = -0.9
y0 = 0
for cat, size in zip(CATEGORY_ORDER, block_sizes):
    ax.add_patch(plt.Rectangle((sidebar_x, y0 - 0.5), 0.5, size, color=CATEGORY_COLORS[cat],
                                clip_on=False, transform=ax.transData))
    ax.text(sidebar_x - 0.3, y0 + size / 2 - 0.5, f"{cat} (n={size})", ha="right", va="center",
            fontsize=10, fontfamily=FONT, color=CATEGORY_COLORS[cat], fontweight="bold")
    if y0 > 0:
        ax.axhline(y0 - 0.5, color="white", linewidth=1.5)
    y0 += size

ax.set_xlim(-0.5, 1.5)
ax.set_ylim(n - 0.5, -0.5)

cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.04, pad=0.04)
cbar.set_label("Log2 Fold Change (vs. Normal)", fontsize=12, fontweight="bold", fontfamily=FONT, labelpad=8)
cbar.ax.tick_params(labelsize=10)
for lbl in cbar.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cbar.outline.set_visible(False)

ax.set_title("Categories of Interest — Log2FC Heatmap\n(all detected proteins, HG/HL vs. Normal)",
              fontsize=15, fontweight="bold", fontfamily=FONT, pad=14)

plt.tight_layout()
for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Categories_overview_heatmap.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("Saved:", path)
