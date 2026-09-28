"""
Volcano Plots — HUVEC-EVs, HG vs Normal & HL vs Normal (filtered dataset)
==============================================================================
Recreated from HUVEC-EVs 01-2026/02-volcano/volcano_HG_v_normal.py and
volcano_HL_v_normal.py, combined into one side-by-side figure, adapted to:
  - Read from the filtered dataset (single-peptide proteins excluded),
    same file used for 4.2.1's PCA/Pearson/heatmap and 4.2.3's RA/ES/ΔES
    analysis — see [[17 - MS Proteomics Data - Filtering & Provenance]]
  - Both comparisons (HG vs Normal, HL vs Normal) share equal x/y axis
    scales for direct visual comparison (house style)
  - Single centered y-axis label instead of duplicated per panel (house
    style); x-axis labels stay separate (different group per panel)
  - Condition colors match the rest of the thesis: HG orange (0.99,0.57,0),
    HL red (0.99,0.14,0)
  - Interactive: opens a matplotlib window (plt.show()) — hover over any
    point to see its protein name (as in the reference scripts), resize
    the window to your preferred dimensions, then use the window's save
    icon (floppy disk in the toolbar) to export the final PNG/SVG.
    RUN THIS LOCALLY, not headless.

Thresholds (unchanged from reference): p < 0.05, |log2FC| > 0.85

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: (manual, via the toolbar's save button) Volcano_HUVEC_EV.png / .svg
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
# NOTE: no matplotlib.use("Agg") here — this script opens an INTERACTIVE window.
# Run it locally to hover over points (protein names) and resize the figure,
# then use the window's save icon (floppy disk in the toolbar) to export.

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
OUT  = os.path.join(BASE, "Proteome Overview", "Volcano")

PROTEIN_COL = "Genes"
SIG_THRESHOLD    = 0.05
LOG2FC_THRESHOLD = 0.85

COMPARISONS = [
    {"group": "HG", "pval_col": "Student's T-test p-value HG_Normal",
     "log2fc_col": "Student's T-test Difference HG_Normal", "color": (0.99, 0.57, 0)},
    {"group": "HL", "pval_col": "Student's T-test p-value HL_Normal",
     "log2fc_col": "Student's T-test Difference HL_Normal", "color": (0.99, 0.14, 0)},
]
NOT_SIG_COLOR = (0.7, 0.7, 0.7)

# -- Load filtered EVs data -----------------------------------------------------
df_all = pd.read_csv(DATA)

# -- Prepare each comparison: drop rows with missing p-value/diff, classify ----
prepared = []
for comp in COMPARISONS:
    d = df_all[[PROTEIN_COL, comp["pval_col"], comp["log2fc_col"]]].dropna().copy()
    d["-log10(p-value)"] = -np.log10(d[comp["pval_col"]])
    d["Significance"] = "Not Significant"
    d.loc[(d[comp["pval_col"]] < SIG_THRESHOLD) & (d[comp["log2fc_col"]] > LOG2FC_THRESHOLD), "Significance"] = "Upregulated"
    d.loc[(d[comp["pval_col"]] < SIG_THRESHOLD) & (d[comp["log2fc_col"]] < -LOG2FC_THRESHOLD), "Significance"] = "Downregulated"
    n_up = (d["Significance"] == "Upregulated").sum()
    n_down = (d["Significance"] == "Downregulated").sum()
    print(f"{comp['group']} vs Normal — proteins: {len(d)} (dropped {len(df_all) - len(d)} missing), "
          f"Upregulated: {n_up}, Downregulated: {n_down}")
    prepared.append(d)

# -- Shared axis limits, so both panels are directly comparable ----------------
x_max = max(d[comp["log2fc_col"]].abs().max() for d, comp in zip(prepared, COMPARISONS)) * 1.1
y_max = max(d["-log10(p-value)"].max() for d in prepared) * 1.08

# -- Figure: two volcano panels side by side ------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 7))
all_texts = []  # for hover lookups, one dict per panel

for ax, comp, d in zip(axes, COMPARISONS, prepared):
    rgb_colors = {
        "Upregulated": comp["color"],
        "Downregulated": comp["color"],
        "Not Significant": NOT_SIG_COLOR,
    }
    for category, color in rgb_colors.items():
        subset = d[d["Significance"] == category]
        ax.scatter(subset[comp["log2fc_col"]], subset["-log10(p-value)"],
                   color=[color], alpha=0.7, s=25, edgecolors="none", picker=True)

    ax.axhline(y=-np.log10(SIG_THRESHOLD), color="gray", linestyle="--", linewidth=1)
    ax.axvline(x=LOG2FC_THRESHOLD, color="gray", linestyle="--", linewidth=1)
    ax.axvline(x=-LOG2FC_THRESHOLD, color="gray", linestyle="--", linewidth=1)

    ax.set_xlim(-x_max, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel(f"Log2 Fold Change ({comp['group']} / Normal)", fontsize=13, fontweight="bold", fontfamily=FONT)
    ax.set_title(f"{comp['group']} vs Normal", fontsize=17, fontweight="bold", fontfamily=FONT)
    ax.tick_params(axis="both", which="major", labelsize=12)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight("bold")
        lbl.set_fontfamily(FONT)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#bbbbbb")

# Bold panel letters (A/B), top-left of each panel - standing house-style
# rule for multi-panel figures (feedback_panel_lettering)
for ax, letter in zip(axes, ["A", "B"]):
    ax.text(-0.10, 1.08, letter, transform=ax.transAxes, fontsize=18,
            fontweight="bold", color="#0b0b0b", fontfamily=FONT,
            va="bottom", ha="left")

# Single centered y-axis label instead of duplicating it on both panels
axes[0].set_ylabel("-Log10(p-value)", fontsize=13, fontweight="bold", fontfamily=FONT)
axes[1].set_ylabel("")

fig.suptitle("Volcano Plots of HUVEC-EVs", fontsize=19, fontweight="bold", fontfamily=FONT)

# -- Hover annotation (per-panel protein lookup) --------------------------------
annots = []
point_lookups = []
for ax, comp, d in zip(axes, COMPARISONS, prepared):
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                         bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"),
                         fontfamily=FONT)
    annot.set_visible(False)
    annots.append(annot)
    lookup = {(round(row[comp["log2fc_col"]], 4), round(row["-log10(p-value)"], 4)): row[PROTEIN_COL]
              for _, row in d.iterrows()}
    point_lookups.append(lookup)

def hover(event):
    for ax, annot, lookup in zip(axes, annots, point_lookups):
        if event.inaxes == ax:
            for scatter_obj in ax.collections:
                cont, ind = scatter_obj.contains(event)
                if cont:
                    pos = scatter_obj.get_offsets()[ind["ind"][0]]
                    x, y = pos
                    closest = min(lookup.keys(), key=lambda p: abs(p[0] - x) + abs(p[1] - y))
                    annot.xy = pos
                    annot.set_text(lookup[closest])
                    annot.get_bbox_patch().set_alpha(0.8)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
            if annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()
        else:
            if annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)

plt.tight_layout()
fig.subplots_adjust(wspace=0.35)  # extra horizontal gap between the two panels
print("Hover over points to see protein names. Resize the window, then use the toolbar's save icon to export.")
plt.show()
