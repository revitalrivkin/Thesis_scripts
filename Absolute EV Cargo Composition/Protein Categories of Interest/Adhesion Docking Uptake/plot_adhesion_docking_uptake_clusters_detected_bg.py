"""
Adhesion / Docking / Uptake DEPs - Clustered by Direction (Option A version)
==============================================================================
Style matches the "Pathology-Enriched Cargo" reference figure (multi-panel
heatmap, category-colored strip above each panel, shared colorbar, rotated
readable gene labels) - see [[19 - Faculty Retreat Presentation]].

Clusters are DATA-DRIVEN (direction of change across HG and HL), not visual/
arbitrary cuts of the dendrogram (tested via fcluster - didn't yield clean
thematic breaks, dominated by one large near-zero cluster). Each cluster's
label is a genuine GO BP enrichment result (Fisher's exact, BH FDR), not
an eyeballed guess:

  - Consistently UP (both HG & HL > 0): Golgi / secretory vesicle transport
  - Consistently DOWN (both HG & HL < 0): cell-matrix adhesion +
    receptor-mediated endocytosis - matches the PI's own network claim
  - Opposite direction (HG and HL move in opposite directions): extracellular
    matrix / basement membrane organization

Input:  ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
        ../category_genesets.csv
Output: Adhesion_Docking_Uptake_clusters_detected_bg.png (+ .svg)
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
CATEGORY = "Adhesion / Docking / Uptake"

CLUSTER_COLORS = {
    "up": "#2a9d5c",
    "down": "#1f77b4",
    "mixed": "#8c564b",
}

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
cat_genes = set(genesets.loc[genesets["Category"] == CATEGORY, "Gene"])
d = d[d.index.isin(cat_genes)]

# -- Direction-based clusters, each internally sorted by hierarchical clustering --
def cluster_order(sub):
    if len(sub) < 2:
        return sub
    Z = linkage(sub[["HG", "HL"]].values, method="average", metric="euclidean")
    return sub.iloc[leaves_list(Z)]

both_up = cluster_order(d[(d["HG"] > 0) & (d["HL"] > 0)])
both_down = cluster_order(d[(d["HG"] < 0) & (d["HL"] < 0)])
mixed = cluster_order(d[~(((d["HG"] > 0) & (d["HL"] > 0)) | ((d["HG"] < 0) & (d["HL"] < 0)))])

PANELS = [
    ("up", f"Consistently Upregulated (n={len(both_up)})\nGolgi / Secretory Vesicle Transport", both_up),
    ("down", f"Consistently Downregulated (n={len(both_down)})\nCell-Matrix Adhesion & Receptor-Mediated Endocytosis", both_down),
    ("mixed", f"Opposite Direction in HG vs. HL (n={len(mixed)})\nECM / Basement Membrane Organization", mixed),
]
# Theme labels come from the cluster GO enrichment against the detected EV proteome
# (04_cluster_GO_enrichment_detected_background.py); FDRs are reported in the text, not the figure.

# Shared color scale across all panels
all_vals = d[["HG", "HL"]].values
vmax = np.percentile(np.abs(all_vals), 99)
vmin = -vmax

# -- Figure: stacked panels, each centered, width proportional to gene count ---
# Axes are placed manually (inches -> figure fractions) so that every panel is
# centered in the plotting area and its width is proportional to its gene count.
max_n = max(len(p[2]) for p in PANELS)
FIG_W, LEFT, RIGHT = 15.3, 0.9, 2.2           # inches; RIGHT reserved for the colorbar
TITLE_H, BLOCK_H = 0.9, 3.4
fig_h = TITLE_H + BLOCK_H * len(PANELS)
gw = (FIG_W - LEFT - RIGHT) / max_n           # inches per gene column
center_x = LEFT + (FIG_W - LEFT - RIGHT) / 2
STRIP_H, HEAT_H = 0.22, 1.45

fig = plt.figure(figsize=(FIG_W, fig_h))


def rect(x0, y_top, w, h):
    """Inches from the top-left -> [left, bottom, width, height] in figure fractions."""
    return [x0 / FIG_W, 1 - (y_top + h) / fig_h, w / FIG_W, h / fig_h]


for i, (key, header, sub) in enumerate(PANELS):
    n = len(sub)
    w = n * gw
    x0 = center_x - w / 2
    y0 = TITLE_H + i * BLOCK_H + 0.65           # room above the strip for the 2-line header

    strip_ax = fig.add_axes(rect(x0, y0, w, STRIP_H))
    strip_ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=CLUSTER_COLORS[key], transform=strip_ax.transAxes))
    strip_ax.axis("off")
    fig.text((x0 - 0.15) / FIG_W, 1 - (y0 - 0.55) / fig_h, "ABC"[i], fontsize=22, fontweight="bold",
             fontfamily=FONT, color="#0b0b0b", va="top", ha="right")
    fig.text(x0 / FIG_W, 1 - (y0 - 0.08) / fig_h, header, fontsize=13, fontweight="bold", fontfamily=FONT,
             color=CLUSTER_COLORS[key], va="bottom", ha="left", linespacing=1.15)

    ax = fig.add_axes(rect(x0, y0 + STRIP_H + 0.25, w, HEAT_H))
    data_T = sub[["HG", "HL"]].values.T
    im = ax.imshow(data_T, aspect="auto", cmap="coolwarm", vmin=vmin, vmax=vmax, interpolation="none")

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["HG", "HL"], fontsize=13, fontweight="bold", fontfamily=FONT)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub.index, rotation=45, ha="right", rotation_mode="anchor",
                        fontsize=10, fontweight="bold", fontfamily=FONT, color="#222222")
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

# -- Shared colorbar on the right ------------------------------------------------
sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(vmin=vmin, vmax=vmax))
cax = fig.add_axes(rect(FIG_W - RIGHT + 0.55, TITLE_H + 0.6, 0.3, BLOCK_H * len(PANELS) * 0.7))
cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
cbar.set_label("Log2 Fold Change\n(vs. Normal)", fontsize=13, fontweight="bold", fontfamily=FONT, labelpad=10)
cbar.ax.tick_params(labelsize=11)
for lbl in cbar.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
cbar.outline.set_visible(False)

fig.suptitle(f"{CATEGORY} DEPs - Clustered by Direction of Change",
             fontsize=20, fontweight="bold", fontfamily=FONT, color="#1f3864", x=center_x / FIG_W, y=0.985, ha="center")

for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Adhesion_Docking_Uptake_clusters_detected_bg.{ext}")
    fig.savefig(path, dpi=300, transparent=True)
    print("Saved:", path)
