"""
PCA — HUVEC-EVs, Normal / HG / HL (filtered dataset, interactive labels)
======================================================================
Recreated from HUVEC-EVs 01-2026/05-PCA/PCA__HUVEC-EV.py, adapted to:
  - Read from the filtered dataset (single-peptide proteins excluded),
    same file used for 4.2.3's RA/ES/ΔES analysis — see
    [[17 - MS Proteomics Data - Filtering & Provenance]]
  - Option (a): raw log2 intensities used as-is, no imputation — matches
    the original script's apparent assumption of already-clean input
  - Interactive: opens a matplotlib window (plt.show()); drag any sample
    label with the mouse to reposition it, then use the window's save
    icon (floppy disk in the toolbar) to export the final PNG/SVG once
    you're happy with the placement. RUN THIS LOCALLY, not headless.

Input:  ../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: (manual, via the toolbar's save button) PCA_HUVEC_EV.png / .svg
"""

import os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
import matplotlib.pyplot as plt
# NOTE: no matplotlib.use("Agg") here — this script opens an INTERACTIVE window.
# Run it locally (not headless) to drag labels by hand, then use the window's
# save icon (floppy disk in the toolbar) to export once labels are placed.

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
OUT  = os.path.join(BASE, "Proteome Overview", "PCA")

SAMPLE_COLS_BY_GROUP = {
    "Normal": ["Normal-1", "Normal-2", "Normal-3"],
    "HG":     ["HG-1", "HG-2"],
    "HL":     ["HL-1", "HL-2", "HL-3"],
}
GROUP_COLORS = {"Normal": (0, 0, 1), "HG": (0.99, 0.57, 0), "HL": (0.99, 0.14, 0)}

# -- Load filtered EVs data, build wide (samples x proteins) matrix ----------
evs = pd.read_csv(DATA)
all_sample_cols = [c for cols in SAMPLE_COLS_BY_GROUP.values() for c in cols]
evs[all_sample_cols] = evs[all_sample_cols].apply(pd.to_numeric, errors="coerce")

# Transpose: rows become samples, columns become proteins (Protein.Group as feature ID)
X_wide = evs.set_index("Protein.Group")[all_sample_cols].T
X_wide.index.name = "Sample"
X_wide = X_wide.reset_index()

sample_to_group = {s: g for g, cols in SAMPLE_COLS_BY_GROUP.items() for s in cols}
X_wide["Group"] = X_wide["Sample"].map(sample_to_group)

# Reorder columns to match original script's expectation: [Sample, Group, protein...]
protein_cols = [c for c in X_wide.columns if c not in ("Sample", "Group")]
df = X_wide[["Sample", "Group"] + protein_cols]

print(f"Samples: {len(df)}, Proteins (features): {len(protein_cols)}")
print(df[["Sample", "Group"]].to_string(index=False))

# -- PCA (option a: raw log2 intensities as-is, no imputation) --------------
X = df[protein_cols].values
y = df["Group"]

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
principalDf = pd.DataFrame(data=principalComponents, columns=["PC1", "PC2"])
explained_variance = pca.explained_variance_ratio_ * 100

# -- Figure -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 5.5))

texts = []
for group, color in GROUP_COLORS.items():
    mask = (y == group).values
    ax.scatter(principalDf.loc[mask, "PC1"], principalDf.loc[mask, "PC2"],
               color=[color], label=group, s=60, edgecolors="white", linewidths=0.6, zorder=3)

    indices = list(principalDf[mask].index)
    for i in indices:
        texts.append(ax.text(principalDf.loc[i, "PC1"], principalDf.loc[i, "PC2"],
                     str(df.loc[i, "Sample"]), fontsize=9, fontfamily=FONT, color="#444444",
                     picker=5))

ax.set_xlabel(f"PC 1 ({explained_variance[0]:.2f}%)", fontsize=13, fontweight="bold", fontfamily=FONT)
ax.set_ylabel(f"PC 2 ({explained_variance[1]:.2f}%)", fontsize=13, fontweight="bold", fontfamily=FONT)
ax.set_title("PCA of HUVEC-EVs", fontsize=17, fontweight="bold", fontfamily=FONT)
legend = ax.legend(prop={"family": FONT})
for text in legend.get_texts():
    text.set_fontfamily(FONT)
ax.tick_params(labelsize=10)
for lbl in ax.get_xticklabels() + ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color("#bbbbbb")

plt.tight_layout()

# -- Interactive label dragging -----------------------------------------------
# Click and drag any sample-name label to reposition it. Release to drop it.
dragging = {"text": None, "dx": 0.0, "dy": 0.0}

def on_press(event):
    if event.inaxes != ax:
        return
    for t in texts:
        contains, _ = t.contains(event)
        if contains:
            tx, ty = t.get_position()
            dragging["text"] = t
            dragging["dx"] = tx - event.xdata
            dragging["dy"] = ty - event.ydata
            break

def on_motion(event):
    if dragging["text"] is None or event.inaxes != ax:
        return
    dragging["text"].set_position((event.xdata + dragging["dx"], event.ydata + dragging["dy"]))
    fig.canvas.draw_idle()

def on_release(event):
    dragging["text"] = None

fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("motion_notify_event", on_motion)
fig.canvas.mpl_connect("button_release_event", on_release)

print("Drag sample labels to reposition. Use the toolbar's save icon to export once done.")
plt.show()
