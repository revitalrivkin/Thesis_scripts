"""
4.2.3.6 Pathology-Depleted Cargo - 3-panel dES heatmap, house style
(mirrors 4.2.3.5's 14_build_pathology_enriched_heatmap_ABC.py)
==============================================================================
Criteria: ES_Normal > 1 AND dES_<condition> < -1, evaluated separately for HG
and HL (same as 01_build_pathology_depleted_cargo_list.py / the Constitutive
Cargo ES computation). Categories = the original 5 categories of interest
from 4.2.2.2 (NOT the Endosomal-Lysosomal Trafficking category added for
Enriched cargo - no equivalent bottom-up category was justified here).

Contaminant handling (Revital, 2026-09-28): a GO-based, reproducible rule
flags likely keratin / cornified-envelope / epidermal-differentiation
contaminants (common airborne/handling MS contamination, not real EV
biology). Flagged genes are forced into "Other" regardless of what GO
category they would otherwise match (this catches e.g. DSC1/DSG1/CDSN/FLG2,
which GO also tags "cell adhesion" but which are skin desmosome proteins,
not real vascular signal). No new category/color is created for them - they
stay grey, but are grouped together within the "Other" block with a thin
separator, so they're visible without being emphasized.

Input:
  ../01_pathology_depleted_cargo_HG.py output (Depleted/HG/01_..._HG.xlsx)
  ../01_pathology_depleted_cargo_HL.py output (Depleted/HL/01_..._HL.xlsx)
  ../../../../Absolute EV Cargo Composition/Protein Categories of Interest/category_genesets.csv
  GO_data (goa_human.gaf.gz, go-basic.obo)
Output: Pathology_Depleted_ABC.png/.svg, Pathology_Depleted_ABC_data.xlsx
"""

import gzip
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

FONT = "Arial"
plt.rcParams["font.family"] = FONT
TITLE_COLOR = "#1f3864"
GREY = "#444444"

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
GENESETS = os.path.join(BASE, "Absolute EV Cargo Composition",
                         "Protein Categories of Interest", "category_genesets.csv")
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"

CATEGORY_PRIORITY = ["EV Biogenesis", "ECM Organization", "Immune-evasion", "Cardioprotective",
                     "Adhesion / Docking / Uptake"]
CATEGORY_ORDER = CATEGORY_PRIORITY + ["Other"]
CATEGORY_COLORS = {
    "EV Biogenesis": "#9467bd", "ECM Organization": "#d95f02", "Immune-evasion": "#7570b3",
    "Cardioprotective": "#e7298a", "Adhesion / Docking / Uptake": "#66a61e", "Other": "#bbbbbb",
}

# ---------------------------------------------------------------- GO-based contaminant rule
CONTAM_TERMS = {
    "GO:0030216": "keratinocyte differentiation",
    "GO:0070268": "cornification",
    "GO:0031424": "keratinization",
}  # "epidermis development" (GO:0008544) deliberately excluded: too broad, catches real
   # signaling proteins that merely play a role in epidermis (e.g. CCN2/CTGF) rather than
   # structural skin contaminants
CONTAM_SUPPLEMENT = {  # classic cornified-envelope/exocrine-secretion MS contaminants not
    "FLG", "FLG2", "SBSN", "HRNR", "DSC1", "DSG1", "CALML5", "LGALS7", "LCN1", "PIP", "DCD", "KPRP",
}  # caught by the narrow GO terms above; same "GO core + hard-coded supplement" approach
   # already used for EV Biogenesis in 4.2.2.2
children = defaultdict(set)
cur = None
with open(os.path.join(GO_DATA, "go-basic.obo"), encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if line in ("[Term]", "[Typedef]"):
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif cur and line.startswith("is_a: "):
            children[line.split()[1]].add(cur)
        elif cur and line.startswith("relationship: part_of "):
            children[line.split()[2]].add(cur)


def descendants(t):
    out, stack = {t}, [t]
    while stack:
        x = stack.pop()
        for c in children.get(x, ()):
            if c not in out:
                out.add(c); stack.append(c)
    return out


contam_go = set().union(*[descendants(t) for t in CONTAM_TERMS])
gene2go = defaultdict(set)
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10 or c[8] != "P" or "NOT" in c[3]:
            continue
        gene2go[c[2]].add(c[4])


def is_contaminant(gene):
    if gene.startswith("KRT") and gene[3:4].isdigit():
        return True
    if gene in CONTAM_SUPPLEMENT:
        return True
    return bool(gene2go.get(gene, set()) & contam_go)


# ---------------------------------------------------------------- category + contaminant assignment
genesets = pd.read_csv(GENESETS)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()


def categorize(gene):
    if is_contaminant(gene):
        return "Other", True
    cats = gene_to_cats.get(gene, set())
    for cat in CATEGORY_PRIORITY:
        if cat in cats:
            return cat, False
    return "Other", False


DATA_DIR = os.path.join(BASE, "MS proteomics data", "Filtered Datasets")
EV_COLS_BY_COND = {"Normal": ["Normal-1", "Normal-2", "Normal-3"], "HG": ["HG-1", "HG-2"],
                   "HL": ["HL-1", "HL-2", "HL-3"]}
CELL_COLS_BY_COND = {"Normal": ["Normal-1", "Normal-2", "Normal-3"], "HG": ["HG-1", "HG-2", "HG-3"],
                     "HL": ["HL-1", "HL-2", "HL-3"]}
EPSILON = 1e-6


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

full = es["Normal"].merge(es["HG"][["Protein.Group", "ES_HG"]], on="Protein.Group", how="inner")
full = full.merge(es["HL"][["Protein.Group", "ES_HL"]], on="Protein.Group", how="inner")
full["dES_HG"] = full["ES_HG"] - full["ES_Normal"]
full["dES_HL"] = full["ES_HL"] - full["ES_Normal"]
full["Gene"] = full["Genes"].astype(str).str.split(";").str[0].str.replace(r"^cRAP-", "", regex=True).str.strip()
print(f"Proteins with ES in all 3 conditions: {len(full)}")

mask_hg = (full["ES_Normal"] > 1) & (full["dES_HG"] < -1)
mask_hl = (full["ES_Normal"] > 1) & (full["dES_HL"] < -1)
hg_genes = set(full.loc[mask_hg, "Gene"])
hl_genes = set(full.loc[mask_hl, "Gene"])
both_genes = hg_genes & hl_genes
print(f"HG-depleted: {len(hg_genes)}  HL-depleted: {len(hl_genes)}  Overlap (Both): {len(both_genes)}")

hg_dES_all = full.set_index("Gene")["dES_HG"]
hl_dES_all = full.set_index("Gene")["dES_HL"]


def build_group(genes, label):
    rows = []
    for g in sorted(genes):
        cat, contam = categorize(g)
        rows.append({"Gene": g, "Category": cat, "Contaminant": contam,
                     "dES_HG": hg_dES_all.get(g, np.nan), "dES_HL": hl_dES_all.get(g, np.nan)})
    df = pd.DataFrame(rows)
    df["cat_rank"] = df["Category"].map({c: i for i, c in enumerate(CATEGORY_ORDER)})
    print(f"\n{label} (n={len(df)}): categories = {df['Category'].value_counts().to_dict()}, "
          f"contaminant-suspected = {int(df['Contaminant'].sum())}")
    return df


hg_only = build_group(hg_genes - both_genes, "HG-depleted only")
hl_only = build_group(hl_genes - both_genes, "HL-depleted only")
both = build_group(both_genes, "Depleted in both HG and HL")

all_data = pd.concat([hg_only.assign(Group="HG only"), hl_only.assign(Group="HL only"),
                      both.assign(Group="Both")], ignore_index=True)
all_data.to_excel(os.path.join(HERE, "Pathology_Depleted_ABC_data.xlsx"), index=False)
print("\nSaved: Pathology_Depleted_ABC_data.xlsx")

# ---------------------------------------------------------------- ordering within each panel
def order_panel(df):
    # category, then contaminant-suspected genes grouped last within "Other", then by dES
    df = df.copy()
    df["contam_rank"] = df["Contaminant"].astype(int)
    sort_col = "dES_HG" if "dES_HG" in df and df["dES_HG"].notna().any() else "dES_HL"
    return df.sort_values(["cat_rank", "contam_rank", sort_col]).reset_index(drop=True)


panels_data = {"A": ("HL-depleted only", order_panel(hl_only), "dES_HL"),
              "B": ("HG-depleted only", order_panel(hg_only), "dES_HG"),
              "C": ("Depleted in both HG and HL", order_panel(both), "dES_HG")}
max_n = max(len(d) for _, d, _ in panels_data.values())

vmax = max(2.0, np.nanpercentile(np.abs(np.concatenate(
    [d[["dES_HG", "dES_HL"]].to_numpy(dtype=float).ravel() for _, d, _ in panels_data.values()])), 95))
norm = mcolors.TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
cmap = plt.cm.RdBu_r

FIG_W, LEFT, RIGHT = 15.3, 0.9, 3.0
TITLE_H, BLOCK_H = 0.9, 3.55
STRIP_H, HEAT_H = 0.22, 1.45
fig_h = TITLE_H + BLOCK_H * len(panels_data)
gw = (FIG_W - LEFT - RIGHT) / max_n
center_x = LEFT + (FIG_W - LEFT - RIGHT) / 2

fig = plt.figure(figsize=(FIG_W, fig_h))


def rect(x0, y_top, w, h):
    return [x0 / FIG_W, 1 - (y_top + h) / fig_h, w / FIG_W, h / fig_h]


im_ref = None
for i, letter in enumerate(["A", "B", "C"]):
    header, sub, own_col = panels_data[letter]
    n = len(sub)
    w = n * gw
    x0 = center_x - w / 2
    y0 = TITLE_H + i * BLOCK_H + 0.80

    strip_ax = fig.add_axes(rect(x0, y0, w, STRIP_H))
    cat_colors = np.array([mcolors.to_rgb(CATEGORY_COLORS[c]) for c in sub["Category"]]).reshape(1, n, 3)
    strip_ax.imshow(cat_colors, aspect="auto")
    strip_ax.set_xticks([]); strip_ax.set_yticks([]); strip_ax.set_xlim(-0.5, n - 0.5)
    boundary = 0
    for cat in CATEGORY_ORDER:
        size = int((sub["Category"] == cat).sum())
        if size == 0:
            continue
        if boundary > 0:
            strip_ax.axvline(boundary - 0.5, color="white", linewidth=1.6, zorder=4)
        boundary += size
    # thin separator inside "Other" between non-contaminant and contaminant-suspected genes
    other_mask = sub["Category"] == "Other"
    if other_mask.any() and sub.loc[other_mask, "Contaminant"].nunique() > 1:
        other_start = other_mask.idxmax()
        n_clean = int(((sub["Category"] == "Other") & (~sub["Contaminant"])).sum())
        sep_x = other_start + n_clean - 0.5
        strip_ax.axvline(sep_x, color="white", linewidth=1.0, linestyle=(0, (2, 1)), zorder=4)

    fig.text((x0 - 0.15) / FIG_W, 1 - (y0 - 0.55) / fig_h, letter, fontsize=24, fontweight="bold",
             fontfamily=FONT, color="#0b0b0b", va="bottom", ha="right")
    fig.text(x0 / FIG_W, 1 - (y0 - 0.08) / fig_h, f"{header}  (n={n})", fontsize=16, fontweight="bold",
             fontfamily=FONT, color=TITLE_COLOR, va="bottom", ha="left")

    ax = fig.add_axes(rect(x0, y0 + STRIP_H + 0.25, w, HEAT_H))
    mat = sub[["dES_HG", "dES_HL"]].to_numpy(dtype=float).T
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
    if other_mask.any() and sub.loc[other_mask, "Contaminant"].nunique() > 1:
        ax.axvline(sep_x, color="white", linewidth=1.0, linestyle=(0, (2, 1)), zorder=4)

cax = fig.add_axes(rect(FIG_W - RIGHT + 0.55, TITLE_H + 0.35, 0.4, BLOCK_H * len(panels_data) * 0.62))
cb = fig.colorbar(im_ref, cax=cax, orientation="vertical")
cb.set_label("\u0394ES", fontsize=17, fontweight="bold", fontfamily=FONT, color=GREY, labelpad=10)
cb.ax.tick_params(labelsize=13, colors=GREY)
cb.outline.set_visible(False)
for lbl in cb.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)

cats_present = [c for c in CATEGORY_ORDER if c in set(all_data["Category"])]
legend_handles = [Patch(facecolor=CATEGORY_COLORS[c], label=c, linewidth=0) for c in cats_present]
leg = fig.legend(handles=legend_handles, loc="upper left",
                  bbox_to_anchor=((FIG_W - RIGHT + 0.5) / FIG_W,
                                  1 - (TITLE_H + 0.35 + BLOCK_H * len(panels_data) * 0.62 + 0.45) / fig_h),
                  bbox_transform=fig.transFigure, fontsize=15.5, frameon=False,
                  handlelength=1.4, labelspacing=1.2, handletextpad=0.7)
for t in leg.get_texts():
    t.set_fontfamily(FONT); t.set_color(GREY)

fig.suptitle("Pathology-Depleted Cargo - \u0394ES by Group", fontsize=22, fontweight="bold",
             fontfamily=FONT, color=TITLE_COLOR, x=center_x / FIG_W, y=0.985, ha="center")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(HERE, f"Pathology_Depleted_ABC.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved: Pathology_Depleted_ABC.png/.svg")

# ---------------------------------------------------------------- summary: informative vs contaminant
print("\n=== Summary: contaminant-suspected vs. informative, per group ===")
for letter, (label, sub, _) in panels_data.items():
    n_contam = int(sub["Contaminant"].sum())
    n_cat = int((sub["Category"] != "Other").sum())
    n_other_clean = len(sub) - n_contam - n_cat
    print(f"{letter} ({label}, n={len(sub)}): {n_cat} in a named category, "
          f"{n_other_clean} Other/non-contaminant, {n_contam} contaminant-suspected")
