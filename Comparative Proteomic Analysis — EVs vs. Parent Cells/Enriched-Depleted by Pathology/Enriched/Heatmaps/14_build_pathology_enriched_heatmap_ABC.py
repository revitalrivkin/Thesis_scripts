"""
4.2.3.5 Pathology-Enriched Cargo - 3-panel dES heatmap, house style (matches the
Adhesion/Docking/Uptake cluster figure: plot_adhesion_docking_uptake_clusters_detected_bg.py)
==============================================================================
Rebuild of 13_build_pathology_enriched_heatmap_dES_horizontal_3panel.py (data/
categories unchanged) with the current thesis house style for multi-panel
figures (2026-09-27, Revital: "go with A/B/C letters, match Adhesion style"):
  - panels stacked vertically, each CENTERED, width proportional to n genes
  - bold A/B/C panel letters top-left of each panel's category strip
  - centered navy suptitle, hyphen not em dash
  - larger gene labels (10pt, up from 8.8pt)
Panel order top to bottom: A = HL-enriched only, B = HG-enriched only,
C = enriched in both.

Input:
  ../../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
  ../../../../MS proteomics data/Filtered Datasets/Cells_df_filtered.csv
  ../../../../Absolute EV Cargo Composition/Protein Categories of Interest/category_genesets.csv
Output: Pathology_Enriched_ABC.png/.svg, Pathology_Enriched_ABC_data.xlsx
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

FONT = "Arial"
plt.rcParams["font.family"] = FONT
TITLE_COLOR = "#1f3864"
GREY = "#444444"

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
GENESETS = os.path.join(BASE, "Absolute EV Cargo Composition",
                         "Protein Categories of Interest", "category_genesets.csv")

EV_COLS_BY_COND = {"Normal": ["Normal-1", "Normal-2", "Normal-3"],
                   "HG": ["HG-1", "HG-2"], "HL": ["HL-1", "HL-2", "HL-3"]}
CELL_COLS_BY_COND = {"Normal": ["Normal-1", "Normal-2", "Normal-3"],
                     "HG": ["HG-1", "HG-2", "HG-3"], "HL": ["HL-1", "HL-2", "HL-3"]}
EPSILON = 1e-6

# ---------------------------------------------------------------- category definitions (unchanged from 13)
TRAFFICKING_GENES = {
    "ABCA3", "ACP2", "ARL8B", "ATP6V0A2", "B4GALT1", "CHMP2A", "CLN3", "DPP4",
    "GLMP", "GNB1", "LAMP1", "LAMP2", "LAPTM4A", "MFSD8", "NPC1", "P2RX4",
    "PI4K2A", "PIP4P1", "SCARB2", "SLC15A4", "SLC17A5", "SLC49A4", "SPNS1",
    "SPPL2A", "TMEM179B", "TMEM59",
    "GNAS", "RALA", "RHOA", "RHOJ", "RIT1",
}
ECM_ADD = {"TINAGL1", "SLIT2"}
IMMUNE_ADD = {"TNFRSF1A", "LTBR", "CFI", "BST1"}
OLD_CATEGORY_PRIORITY = ["ECM Organization", "Immune-evasion", "Cardioprotective", "Adhesion / Docking / Uptake"]
CATEGORY_ORDER = ["Endosomal-Lysosomal Trafficking", "ECM Organization", "Immune-evasion",
                  "Cardioprotective", "Adhesion / Docking / Uptake", "Other"]
CATEGORY_COLORS = {
    "Endosomal-Lysosomal Trafficking": "#1b9e77", "ECM Organization": "#d95f02",
    "Immune-evasion": "#7570b3", "Cardioprotective": "#e7298a",
    "Adhesion / Docking / Uptake": "#66a61e", "Other": "#bbbbbb",
}


def log2_to_raw(series):
    s = series.copy()
    s[s <= 13] = np.nan
    return np.power(2, s)


def compute_ra(df, cols):
    df = df.copy()
    df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")
    raw = pd.DataFrame(index=df.index)
    for col in cols:
        raw[col] = log2_to_raw(df[col])
    ra = pd.DataFrame(index=df.index)
    for col in cols:
        col_sum = raw[col].sum(skipna=True)
        ra[col + "_RA"] = raw[col] / col_sum
    df["mean_RA"] = ra[[c + "_RA" for c in cols]].mean(axis=1)
    return df[["Protein.Group", "Genes", "mean_RA"]].rename(columns={"mean_RA": "RA"})


evs_raw = pd.read_csv(os.path.join(DATA_DIR, "EVs_df_filtered.csv"))
cells_raw = pd.read_csv(os.path.join(DATA_DIR, "Cells_df_filtered.csv"))
ev_ra = {c: compute_ra(evs_raw, EV_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}
cell_ra = {c: compute_ra(cells_raw, CELL_COLS_BY_COND[c]) for c in ["Normal", "HG", "HL"]}

es = {}
for c in ["Normal", "HG", "HL"]:
    m = ev_ra[c].merge(cell_ra[c][["Protein.Group", "RA"]], on="Protein.Group", how="inner",
                        suffixes=("_ev", "_cell"))
    m = m[m["RA_ev"].notna() & m["RA_cell"].notna() & (m["RA_cell"] > 0)].copy()
    m["ES"] = np.log2((m["RA_ev"] + EPSILON) / (m["RA_cell"] + EPSILON))
    es[c] = m[["Protein.Group", "Genes", "ES"]].rename(columns={"ES": f"ES_{c}"})

df = es["Normal"].merge(es["HG"][["Protein.Group", "ES_HG"]], on="Protein.Group", how="inner")
df = df.merge(es["HL"][["Protein.Group", "ES_HL"]], on="Protein.Group", how="inner")
df["dES_HG"] = df["ES_HG"] - df["ES_Normal"]
df["dES_HL"] = df["ES_HL"] - df["ES_Normal"]
df["Gene"] = df["Genes"].astype(str).str.split(";").str[0].str.replace(r"^cRAP-", "", regex=True).str.strip()

mask_hg = (df["ES_HG"] > 1) & (df["dES_HG"] > 1)
mask_hl = (df["ES_HL"] > 1) & (df["dES_HL"] > 1)
union = df[mask_hg | mask_hl].copy()
union["Enriched_in"] = np.select(
    [mask_hg.loc[union.index] & mask_hl.loc[union.index], mask_hg.loc[union.index]],
    ["Both", "HG only"], default="HL only")
print(f"HG-enriched: {mask_hg.sum()}  HL-enriched: {mask_hl.sum()}  "
      f"Union: {len(union)}  Overlap: {(mask_hg & mask_hl).sum()}")

genesets = pd.read_csv(GENESETS)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()


def categorize(gene):
    if gene in TRAFFICKING_GENES:
        return "Endosomal-Lysosomal Trafficking"
    if gene in ECM_ADD:
        return "ECM Organization"
    if gene in IMMUNE_ADD:
        return "Immune-evasion"
    cats = gene_to_cats.get(gene, set())
    for cat in OLD_CATEGORY_PRIORITY:
        if cat in cats:
            return cat
    return "Other"


union["Category"] = union["Gene"].apply(categorize)
union["cat_rank"] = union["Category"].map({c: i for i, c in enumerate(CATEGORY_ORDER)})

out_cols = ["Gene", "Category", "Enriched_in", "ES_Normal", "ES_HG", "ES_HL", "dES_HG", "dES_HL"]
union.sort_values(["Enriched_in", "cat_rank"])[out_cols].to_excel(
    os.path.join(HERE, "Pathology_Enriched_ABC_data.xlsx"), index=False)

vmax = max(2.0, np.nanpercentile(np.abs(union[["dES_HG", "dES_HL"]].to_numpy()), 95))
norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
cmap = plt.cm.RdBu_r

# panel order top to bottom: A = HL only, B = HG only, C = Both
GROUPS = [("A", "HL only", "HL-enriched only"), ("B", "HG only", "HG-enriched only"),
          ("C", "Both", "Enriched in both HG and HL")]
panels = []
for letter, key, label in GROUPS:
    sub = union[union["Enriched_in"] == key].sort_values(["cat_rank", "dES_HG"], ascending=[True, False])
    panels.append((letter, label, sub.reset_index(drop=True)))
max_n = max(len(p[2]) for p in panels)

# ---------------------------------------------------------------- layout (Adhesion-cluster house style)
FIG_W, LEFT, RIGHT = 15.3, 0.9, 3.0     # RIGHT wider: needs room for colorbar + category legend
TITLE_H, BLOCK_H = 0.9, 3.55
STRIP_H, HEAT_H = 0.22, 1.45
fig_h = TITLE_H + BLOCK_H * len(panels)
gw = (FIG_W - LEFT - RIGHT) / max_n
center_x = LEFT + (FIG_W - LEFT - RIGHT) / 2

fig = plt.figure(figsize=(FIG_W, fig_h))


def rect(x0, y_top, w, h):
    return [x0 / FIG_W, 1 - (y_top + h) / fig_h, w / FIG_W, h / fig_h]


im_ref = None
for i, (letter, header, sub) in enumerate(panels):
    n = len(sub)
    w = n * gw
    x0 = center_x - w / 2
    y0 = TITLE_H + i * BLOCK_H + 0.80

    strip_ax = fig.add_axes(rect(x0, y0, w, STRIP_H))
    cat_colors = np.array([mcolors.to_rgb(CATEGORY_COLORS[c]) for c in sub["Category"]]).reshape(1, n, 3)
    strip_ax.imshow(cat_colors, aspect="auto")
    strip_ax.set_xticks([]); strip_ax.set_yticks([])
    strip_ax.set_xlim(-0.5, n - 0.5)
    boundary = 0
    for cat in CATEGORY_ORDER:
        size = int((sub["Category"] == cat).sum())
        if size == 0:
            continue
        if boundary > 0:
            strip_ax.axvline(boundary - 0.5, color="white", linewidth=1.6, zorder=4)
        boundary += size

    fig.text((x0 - 0.15) / FIG_W, 1 - (y0 - 0.55) / fig_h, letter, fontsize=24, fontweight="bold",
             fontfamily=FONT, color="#0b0b0b", va="bottom", ha="right")
    fig.text(x0 / FIG_W, 1 - (y0 - 0.08) / fig_h, f"{header}  (n={n})", fontsize=16, fontweight="bold",
             fontfamily=FONT, color=TITLE_COLOR, va="bottom", ha="left")

    ax = fig.add_axes(rect(x0, y0 + STRIP_H + 0.25, w, HEAT_H))
    mat = sub[["dES_HG", "dES_HL"]].to_numpy().T
    im = ax.imshow(mat, aspect="auto", cmap=cmap, norm=norm)
    im_ref = im
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["\u0394ES HG", "\u0394ES HL"], fontsize=15, fontweight="bold", fontfamily=FONT, color=GREY)
    ax.set_xticks(range(n))
    ax.set_xticklabels(sub["Gene"], rotation=45, ha="right", rotation_mode="anchor",
                        fontsize=11, fontweight="bold", fontfamily=FONT, color="#222222")
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.5, n - 0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    boundary = 0
    for cat in CATEGORY_ORDER:
        size = int((sub["Category"] == cat).sum())
        if size == 0:
            continue
        if boundary > 0:
            ax.axvline(boundary - 0.5, color="white", linewidth=1.6, zorder=4)
        boundary += size

# shared vertical colorbar
cax = fig.add_axes(rect(FIG_W - RIGHT + 0.55, TITLE_H + 0.35, 0.4, BLOCK_H * len(panels) * 0.62))
cb = fig.colorbar(im_ref, cax=cax, orientation="vertical")
cb.set_label("\u0394ES", fontsize=17, fontweight="bold", fontfamily=FONT, color=GREY, labelpad=10)
cb.ax.tick_params(labelsize=13, colors=GREY)
cb.outline.set_visible(False)
for lbl in cb.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)

# shared category legend, below the colorbar
cats_present = [c for c in CATEGORY_ORDER if c in set(union["Category"])]
legend_handles = [Patch(facecolor=CATEGORY_COLORS[c], label=c, linewidth=0) for c in cats_present]
leg = fig.legend(handles=legend_handles, loc="upper left",
                  bbox_to_anchor=((FIG_W - RIGHT + 0.5) / FIG_W, 1 - (TITLE_H + 0.35 + BLOCK_H * len(panels) * 0.62 + 0.45) / fig_h),
                  bbox_transform=fig.transFigure, fontsize=15.5, frameon=False,
                  handlelength=1.4, labelspacing=1.2, handletextpad=0.7)
for t in leg.get_texts():
    t.set_fontfamily(FONT); t.set_color(GREY)

fig.suptitle("Pathology-Enriched Cargo - \u0394ES by Group", fontsize=22, fontweight="bold",
             fontfamily=FONT, color=TITLE_COLOR, x=center_x / FIG_W, y=0.985, ha="center")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(HERE, f"Pathology_Enriched_ABC.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved: Pathology_Enriched_ABC.png/.svg")
