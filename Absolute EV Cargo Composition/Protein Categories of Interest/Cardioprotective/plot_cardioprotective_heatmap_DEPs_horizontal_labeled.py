"""
Cardioprotective Category - Log2FC Heatmap by Subcategory, DEPs only, gene-labeled (thesis version)
==============================================================================
Finalized main-text figure for the Cardioprotective category. Combines:
  - DEP-only scope (plot_cardioprotective_heatmap_DEPs.py) - proteins that
    are DEPs in HG or HL (p<0.05 & |log2FC|>0.85), not all detected. The
    all-detected version showed mostly-flat signal with sparse standouts;
    DEP-only isolates the real story and was chosen as the main-text figure
    (all-detected version held in reserve as supplementary if needed).
  - Horizontal layout (plot_cardioprotective_heatmap_horizontal.py) -
    proteins along the x-axis (grouped into subcategory blocks, no labels),
    HG/HL as the 2 rows. Presentation-friendly aspect ratio.

Rows (now the x-axis) hierarchically clustered within each subcategory
block (average-linkage, Euclidean - same method as the 4.2.1.3 whole-
proteome heatmap).

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
        ../category_genesets.csv
Thesis version (2026-09-26): subcategory strips/names moved ABOVE the heatmap and gene labels added
below it. Each protein is shown once, in the subcategory chosen by SUBCATEGORY_PRIORITY
(subcategories overlap, so block n's differ from the non-exclusive counts).
Output: Cardioprotective_heatmap_DEPs_horizontal_labeled.png (+ .svg)
"""

import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")

SIG_THRESHOLD, LOG2FC_THRESHOLD = 0.05, 0.85

SUBCATEGORY_ORDER = ["Ca2+ homeostasis", "AMPK signaling", "Antioxidants/Redox", "Glycolysis", "UPR", "HSPs"]
SUBCATEGORY_COLORS = {
    "Ca2+ homeostasis": "#ff7f0e",
    "AMPK signaling": "#bcbd22",
    "Antioxidants/Redox": "#17becf",
    "Glycolysis": "#e377c2",
    "UPR": "#7f7f7f",
    "HSPs": "#b8860b",
}
SUBCATEGORY_PRIORITY = ["Glycolysis", "Ca2+ homeostasis", "AMPK signaling", "UPR", "HSPs", "Antioxidants/Redox"]

# -- Load data, compute log2FC + DEP flags --------------------------------------
df = pd.read_csv(DATA)
d = df[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
        "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d = d.rename(columns={"Student's T-test Difference HG_Normal": "HG",
                       "Student's T-test Difference HL_Normal": "HL",
                       "Student's T-test p-value HG_Normal": "p_HG",
                       "Student's T-test p-value HL_Normal": "p_HL"})
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")

is_dep_hg = (d["p_HG"] < SIG_THRESHOLD) & (d["HG"].abs() > LOG2FC_THRESHOLD)
is_dep_hl = (d["p_HL"] < SIG_THRESHOLD) & (d["HL"].abs() > LOG2FC_THRESHOLD)
d = d[(is_dep_hg.fillna(False)) | (is_dep_hl.fillna(False))]

genesets = pd.read_csv(GENESETS)
cardio = genesets[genesets["Category"] == "Cardioprotective"]
gene_to_subcats = cardio.groupby("Gene")["Subcategory"].apply(set).to_dict()

def assign_subcategory(gene):
    subs = gene_to_subcats.get(gene, set())
    for sub in SUBCATEGORY_PRIORITY:
        if sub in subs:
            return sub
    return None

d["Subcategory"] = [assign_subcategory(g) for g in d.index]
d = d[d["Subcategory"].notna()]
print(f"Cardioprotective DEPs: {len(d)}")
print("Rows per subcategory:")
print(d["Subcategory"].value_counts())

# -- Cluster within each subcategory block ---------------------------------------
ordered_blocks = []
for sub in SUBCATEGORY_ORDER:
    block = d[d["Subcategory"] == sub][["HG", "HL"]]
    if len(block) < 3:
        ordered_blocks.append(block)
        continue
    Z = linkage(block.values, method="average", metric="euclidean")
    order = leaves_list(Z)
    ordered_blocks.append(block.iloc[order])

plot_df = pd.concat(ordered_blocks)
block_sizes = [len(b) for b in ordered_blocks]
print("Block sizes (plot order):", dict(zip(SUBCATEGORY_ORDER, block_sizes)))

# -- Plot (transposed: proteins on x, HG/HL on y) --------------------------------
n = len(plot_df)
fig_w = max(9, n * 0.22 + 2.6)
fig, ax = plt.subplots(figsize=(fig_w, 5.6))

data_T = plot_df[["HG", "HL"]].values.T
vmax = np.percentile(np.abs(data_T), 99)
vmin = -vmax
im = ax.imshow(data_T, aspect="auto", cmap="coolwarm", vmin=vmin, vmax=vmax, interpolation="none")

ax.set_yticks([0, 1])
ax.set_yticklabels(["HG", "HL"], fontsize=17, fontweight="bold", fontfamily=FONT)
ax.set_xticks(range(n))
ax.set_xticklabels(plot_df.index, rotation=45, ha="right", rotation_mode="anchor",
                    fontsize=10, fontweight="bold", fontfamily=FONT, color="#222222")
ax.tick_params(axis="x", length=0, pad=3)

# subcategory strip + name ABOVE the heatmap
strip_y, strip_h = -0.95, 0.32
x0 = 0
for sub, size in zip(SUBCATEGORY_ORDER, block_sizes):
    if size == 0:
        continue
    ax.add_patch(plt.Rectangle((x0 - 0.5, strip_y), size, strip_h, color=SUBCATEGORY_COLORS[sub],
                                clip_on=False, transform=ax.transData))
    ax.text(x0 + size / 2 - 0.5, strip_y - 0.12, f"{sub}\n(n={size})", ha="center", va="bottom",
            fontsize=12.5, fontfamily=FONT, color=SUBCATEGORY_COLORS[sub], fontweight="bold")
    if x0 > 0:
        ax.axvline(x0 - 0.5, color="white", linewidth=1.5)
    x0 += size

ax.set_ylim(1.5, -0.5)
ax.set_xlim(-0.5, n - 0.5)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False)

cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.03, pad=0.02)
cbar.set_label("Log2FC\n(vs. Normal)", fontsize=14, fontweight="bold", fontfamily=FONT, labelpad=8)
cbar.ax.tick_params(labelsize=13)
for lbl in cbar.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cbar.outline.set_visible(False)

ax.set_title("Cardioprotective DEPs - Grouped by Subcategory",
              fontsize=20, fontweight="bold", fontfamily=FONT, color="#1f3864", pad=105)

plt.tight_layout()
for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Cardioprotective_heatmap_DEPs_horizontal_labeled.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("Saved:", path)
