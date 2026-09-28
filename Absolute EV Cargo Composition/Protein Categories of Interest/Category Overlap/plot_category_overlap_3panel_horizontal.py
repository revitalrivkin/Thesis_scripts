"""
Category overlap - three stacked, centered horizontal heatmaps (layout of Adhesion_Docking_Uptake_clusters_detected_bg.png)
==============================================================================
A: Adhesion/Docking/Uptake x Cardioprotective
B: Adhesion/Docking/Uptake x ECM Organization
C: Adhesion/Docking/Uptake x EV Biogenesis

Same layout as the Adhesion cluster figure: three panels stacked vertically, each centered with its width
proportional to its number of proteins, a two-colour strip (the two categories) above each panel with the
panel letter to its left, proteins along the x axis (gene labels below), HG and HL as the two rows, one
shared vertical colour bar on the right, centered title.

All detected proteins that belong to both categories (not only DEPs). Proteins are hierarchically clustered
within each panel by their HG and HL log2 fold changes (average linkage, Euclidean distance). One shared colour
scale (99th percentile of |log2FC| over the three sets).

Designed at 10.4 in wide (scaled to about 87% on a 9 in landscape Word page); gene labels are tilted 45 degrees like the reference.
Output: Category_overlap_3panel_horizontal.png (+ .svg)
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")

ADH = "Adhesion / Docking / Uptake"
COLORS = {ADH: "#1f77b4", "Cardioprotective": "#d62728", "ECM Organization": "#8c564b", "EV Biogenesis": "#9467bd"}
PANELS = [("A", "Cardioprotective"), ("B", "ECM Organization"), ("C", "EV Biogenesis")]

FIG_W, FIG_H = 10.4, 6.3
LEFT, RIGHT = 0.7, 1.15          # inches: left margin (panel letters), right margin (colour bar)
TITLE_H = 0.5
BLOCK_H = 1.85                    # height of one panel block (header, strip, heatmap, labels)
STRIP_H, HM_H = 0.1, 0.66
LABEL_PT = 6.0

df = pd.read_csv(DATA)
d = df[["Genes", "Student's T-test Difference HG_Normal", "Student's T-test Difference HL_Normal"]].copy()
d.columns = ["Genes", "HG", "HL"]
d = d.dropna().set_index("Genes")
gs = pd.read_csv(GENESETS)
sets = {c: set(gs.loc[gs["Category"] == c, "Gene"]) & set(d.index)
        for c in [ADH, "Cardioprotective", "ECM Organization", "EV Biogenesis"]}

data = {}
for letter, other in PANELS:
    genes = sets[ADH] & sets[other]
    sub = d.loc[sorted(genes), ["HG", "HL"]]
    order = leaves_list(linkage(sub.values, method="average", metric="euclidean"))
    data[letter] = sub.iloc[order]
    print(f"{letter}: {ADH} x {other}: {len(sub)} proteins")

vmax = np.percentile(np.abs(np.concatenate([v.values.ravel() for v in data.values()])), 99)
max_n = max(len(v) for v in data.values())
pitch = (FIG_W - LEFT - RIGHT) / max_n           # inches per protein column
center_x = LEFT + (FIG_W - LEFT - RIGHT) / 2

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)


def rect(x, y_top, w, h):
    return [x / FIG_W, 1 - (y_top + h) / FIG_H, w / FIG_W, h / FIG_H]


fig.text(center_x / FIG_W, 1 - 0.24 / FIG_H, "Proteins Shared Between Categories", fontsize=15, fontweight="bold",
         fontfamily=FONT, color="#1f3864", ha="center", va="center")

im = None
for i, (letter, other) in enumerate(PANELS):
    sub = data[letter]
    n = len(sub)
    w = n * pitch
    x0 = center_x - w / 2
    y0 = TITLE_H + i * BLOCK_H + 0.12

    fig.text((x0 - 0.12) / FIG_W, 1 - (y0 + 0.02) / FIG_H, letter, fontsize=18, fontweight="bold", fontfamily=FONT,
             color="#0b0b0b", ha="right", va="top")
    fig.text(x0 / FIG_W, 1 - (y0 + 0.02) / FIG_H,
             f"{ADH} \u00d7 {other} (n={n})", fontsize=10, fontweight="bold", fontfamily=FONT,
             color="#1f3864", ha="left", va="top")
    sa = fig.add_axes(rect(x0, y0 + 0.26, w, STRIP_H))
    sa.add_patch(plt.Rectangle((0, 0), 0.5, 1, color=COLORS[ADH], transform=sa.transAxes))
    sa.add_patch(plt.Rectangle((0.5, 0), 0.5, 1, color=COLORS[other], transform=sa.transAxes))
    sa.axis("off")

    ax = fig.add_axes(rect(x0, y0 + 0.26 + STRIP_H + 0.07, w, HM_H))
    im = ax.imshow(sub.values.T, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax, interpolation="none")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["HG", "HL"], fontsize=10, fontweight="bold", fontfamily=FONT)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub.index, rotation=45, ha="right", rotation_mode="anchor", fontsize=LABEL_PT, fontweight="bold",
                       fontfamily=FONT, color="#222222")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

cax = fig.add_axes(rect(FIG_W - RIGHT + 0.4, (FIG_H - 2.8) / 2, 0.16, 2.8))
cb = fig.colorbar(im, cax=cax, orientation="vertical")
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=8, length=2)
for lbl in cb.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cb.set_label("Log2 Fold Change\n(vs. Normal)", fontsize=9, fontweight="bold", fontfamily=FONT, labelpad=6)

for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Category_overlap_3panel_horizontal.{ext}")
    fig.savefig(path, dpi=300)
    print("Saved:", path)
