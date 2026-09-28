"""
Immune-evasion and EV Biogenesis Categories - Log2FC Heatmap, DEPs only, gene-labeled (thesis version)
==============================================================================
Combined two-block figure for the two smallest categories (Immune-evasion: 15 DEPs,
EV Biogenesis: 2 DEPs). DEP = Student's t-test p < 0.05 and |log2FC| > 0.85 vs Normal in
HG and/or HL. Each protein is shown once (EV Biogenesis takes priority, as in the category
priority order). Within each block, proteins are ordered by hierarchical clustering of their
HG and HL fold changes (average linkage, Euclidean distance).

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
        ../category_genesets.csv
Output: Immune_evasion_EV_Biogenesis_heatmap_DEPs_horizontal_labeled.png (+ .svg)
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
BLOCKS = [("Immune-evasion", "#2ca02c"), ("EV Biogenesis", "#9467bd")]   # plot order
PRIORITY = ["EV Biogenesis", "Immune-evasion"]                          # one protein, one block

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
d = d[(is_dep_hg.fillna(False)) | (is_dep_hl.fillna(False))][["HG", "HL"]]

gs = pd.read_csv(GENESETS)
sets = {c: set(gs.loc[gs["Category"] == c, "Gene"]) for c in PRIORITY}
assigned = {}
for g in d.index:
    for c in PRIORITY:
        if g in sets[c]:
            assigned[g] = c
            break
d["Block"] = [assigned.get(g) for g in d.index]
d = d[d["Block"].notna()]

blocks = []
for name, _ in BLOCKS:
    sub = d[d["Block"] == name][["HG", "HL"]]
    if len(sub) >= 3:
        sub = sub.iloc[leaves_list(linkage(sub.values, method="average", metric="euclidean"))]
    blocks.append(sub)
    print(name, len(sub))
plot_df = pd.concat(blocks)
block_sizes = [len(b) for b in blocks]

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

strip_y, strip_h = -0.95, 0.32
x0 = 0
for (name, color), size in zip(BLOCKS, block_sizes):
    ax.add_patch(plt.Rectangle((x0 - 0.5, strip_y), size, strip_h, color=color, clip_on=False, transform=ax.transData))
    ax.text(x0 + size / 2 - 0.5, strip_y - 0.12, f"{name}\n(n={size})", ha="center", va="bottom",
            fontsize=12.5, fontfamily=FONT, color=color, fontweight="bold")
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

ax.set_title("Immune-evasion and EV Biogenesis DEPs",
              fontsize=20, fontweight="bold", fontfamily=FONT, color="#1f3864", pad=105)

plt.tight_layout()
for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Immune_evasion_EV_Biogenesis_heatmap_DEPs_horizontal_labeled.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("Saved:", path)
