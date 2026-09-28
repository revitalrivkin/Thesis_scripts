"""
Categories of Interest — Simple Explainer Graphic (for PI communication)
==============================================================================
Minimal version (2026-08-21): title + one box per category only. Big
category name, big protein count. No descriptions, no overlap callout, no
methodology caption — those were cut to minimize text/whitespace per
Revital's request. (Overlap finding still available separately if wanted
back in — see category_overlap analysis in note 12.)

Input:  category_genesets.csv
        ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: Categories_explainer_graphic.png (+ .svg)
"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "category_genesets.csv")

# (display name with manual line breaks, lookup key for counts, color)
CATEGORIES = [
    ("Adhesion / Docking\n/ Uptake", "Adhesion / Docking / Uptake", "#1f77b4"),
    ("EV Biogenesis", "EV Biogenesis", "#9467bd"),
    ("Immune-evasion", "Immune-evasion", "#2ca02c"),
    ("ECM\nOrganization", "ECM Organization", "#8c564b"),
    ("Cardioprotective", "Cardioprotective", "#d62728"),
]

df = pd.read_csv(DATA)
detected = set(df["Genes"])
genesets = pd.read_csv(GENESETS)
gs_detected = genesets[genesets["Gene"].isin(detected)]
counts = gs_detected.groupby("Category")["Gene"].nunique().to_dict()

# -- Layout: 5 tall boxes filling the frame, minimal margins --------------------
fig, ax = plt.subplots(figsize=(16, 5))
ax.set_xlim(0, 16)
ax.set_ylim(0, 5)
ax.axis("off")

card_w, gap = 2.95, 0.2
n = len(CATEGORIES)
total_w = n * card_w + (n - 1) * gap
x0 = (16 - total_w) / 2
y0, card_h = 0.15, 4.0

for i, (display_name, key, color) in enumerate(CATEGORIES):
    x = x0 + i * (card_w + gap)
    ax.add_patch(FancyBboxPatch((x, y0), card_w, card_h, boxstyle="round,pad=0.02,rounding_size=0.1",
                                 linewidth=2.5, edgecolor=color, facecolor="white", zorder=2))
    ax.text(x + card_w / 2, y0 + card_h * 0.62, display_name, ha="center", va="center",
            fontsize=15.5, fontweight="bold", color=color, fontfamily=FONT, linespacing=1.3,
            zorder=3)
    n_prot = counts.get(key, 0)
    ax.text(x + card_w / 2, y0 + card_h * 0.28, str(n_prot), ha="center", va="center",
            fontsize=42, fontweight="bold", color=color, fontfamily=FONT, zorder=3)
    ax.text(x + card_w / 2, y0 + card_h * 0.13, "proteins", ha="center", va="center",
            fontsize=13, color=color, fontfamily=FONT, zorder=3)

ax.set_title("Protein Categories of Interest in HUVEC-EV Cargo",
              fontsize=21, fontweight="bold", fontfamily=FONT, pad=8)

plt.tight_layout(pad=0.6)
for ext in ("png", "svg"):
    path = os.path.join(HERE, f"Categories_explainer_graphic.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print("Saved:", path)
