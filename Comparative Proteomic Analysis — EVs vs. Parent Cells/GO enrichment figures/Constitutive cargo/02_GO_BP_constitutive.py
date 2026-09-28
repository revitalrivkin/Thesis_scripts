"""
GO BP Enrichment — Constitutive EV cargo (filtered dataset, finalized ±1 threshold)
=====================================================================================
Recreated from 05-2026_Comparative_Analysis/13-Constitutive_cargo/GO_annotations/
02_GO_BP_constitutive.py — same methodology, updated to the filtered dataset and
the finalized constitutive cargo definition (425 proteins, see
01_build_constitutive_cargo_list.py).

Criteria (applied upstream):
  ES_Normal > 1
  |dES_HG| < 1  AND  |dES_HL| < 1

GO parameters (unchanged from reference):
  Ontology   : Biological Process (BP)
  Test       : Fisher's exact, one-sided, BH-corrected FDR
  Background : all human genes in UniProtKB-GOA
  Min genes per term: 3 | Min depth: 4 | Top terms: 15

Output: GO_BP_constitutive.png/.svg, GO_BP_constitutive_results.xlsx
"""

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from goatools.obo_parser import GODag

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["mathtext.default"] = "regular"  # in case any mathtext slips through

BASE    = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results\Comparative Proteomic Analysis — EVs vs. Parent Cells\GO enrichment figures\Constitutive cargo"
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"
CARGO   = os.path.join(BASE, "01_constitutive_cargo.xlsx")
OUT     = BASE

MIN_GENES = 3
MIN_DEPTH = 4
TOP_N     = 15
DOT_SCALE = 10
GREY      = "#444444"
BLUE_TITLE = (0.15, 0.35, 0.85)

# -- 1. Load constitutive cargo genes -----------------------------------------
df_cargo = pd.read_excel(CARGO)
query_genes = set(df_cargo["Gene"].dropna().str.strip())
print(f"Constitutive cargo genes: {len(query_genes):,}")

# -- 2. Parse GO BP annotations ------------------------------------------------
print("Parsing GO annotations ...")
gene2go_raw = {}
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), 'rt', encoding='utf-8') as fh:
    for line in fh:
        if line.startswith('!'): continue
        cols = line.strip().split('\t')
        if len(cols) < 10 or cols[8] != 'P' or 'NOT' in cols[3]: continue
        gene2go_raw.setdefault(cols[2], set()).add(cols[4])

background = frozenset(gene2go_raw.keys())
godag = GODag(os.path.join(GO_DATA, "go-basic.obo"))

parent_map = {}
for gid, term in godag.items():
    if term.namespace != 'biological_process': continue
    pids = set()
    for p in term.parents:
        pid = p.item_id if hasattr(p, 'item_id') else p.id
        if pid in godag and godag[pid].namespace == 'biological_process':
            pids.add(pid)
    parent_map[gid] = pids

def propagate(go_ids):
    visited = {gid for gid in go_ids if gid in parent_map}
    q = list(visited)
    while q:
        gid = q.pop()
        for pid in parent_map.get(gid, ()):
            if pid not in visited:
                visited.add(pid); q.append(pid)
    return visited

gene2go = {g: propagate(r) for g, r in gene2go_raw.items()}
term2bg = {}
for g, go_ids in gene2go.items():
    for gid in go_ids:
        term2bg.setdefault(gid, set()).add(g)

# -- 3. Fisher's exact + BH correction -----------------------------------------
print("Running GO enrichment ...")
query   = {g for g in query_genes if g in background}
n_query = len(query)
n_bg    = len(background)
print(f"Query genes found in background: {n_query} / {len(query_genes)}")

rows = []
for gid, bg_genes in term2bg.items():
    if gid not in godag or godag[gid].depth < MIN_DEPTH: continue
    overlap = query & bg_genes
    k = len(overlap)
    if k < MIN_GENES: continue
    K = len(bg_genes)
    table = [[k, n_query - k], [K - k, n_bg - K - (n_query - k)]]
    _, pval = fisher_exact(table, alternative='greater')
    rows.append({"GO_ID": gid, "Term": godag[gid].name,
                 "k": k, "depth": godag[gid].depth, "p_value": pval})

res_all = pd.DataFrame(rows)
_, res_all["p_adj"], _, _ = multipletests(res_all["p_value"], method="fdr_bh")
res_all = res_all.sort_values("p_adj").reset_index(drop=True)

n_sig = (res_all["p_adj"] < 0.05).sum()
print(f"  {n_sig} significant terms (FDR<0.05), showing top {TOP_N}")

res = res_all.head(TOP_N)
res_all.to_excel(os.path.join(OUT, "GO_BP_constitutive_results.xlsx"), index=False)

# -- 4. Figure ------------------------------------------------------------------
from matplotlib.lines import Line2D

def neg_log10(series):
    return -np.log10(pd.Series(series).clip(lower=1e-300).values)

nl_vals = neg_log10(res["p_adj"])
vmax    = nl_vals.max()
norm    = mcolors.Normalize(vmin=0, vmax=vmax)
cmap    = plt.cm.Blues

sub = res.iloc[::-1].reset_index(drop=True)
n   = len(sub)
y   = np.arange(n)
nl  = neg_log10(sub["p_adj"])

fig = plt.figure(figsize=(13, max(5.5, n * 0.62 + 2.8)))
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[12, 0.55], wspace=0.05,
                        left=0.03, right=0.97, top=0.88, bottom=0.12)

ax_dot = fig.add_subplot(gs[0])
ax_cb  = fig.add_subplot(gs[1])

sc = ax_dot.scatter(nl, y, s=sub["k"] * DOT_SCALE, c=nl, cmap=cmap, norm=norm,
                     edgecolors="#888888", linewidths=0.4, zorder=3)

FDR_LABEL = "-log10(FDR)"   # plain ASCII - unicode subscripts aren't in Arial's glyph set

ax_dot.set_yticks(y)
ax_dot.set_yticklabels(sub["Term"], fontsize=19.5, color=GREY, fontfamily=FONT)
ax_dot.set_ylim(-0.8, n - 0.2)
ax_dot.set_xlim(max(0, nl.min() - 2), vmax * 1.12)
ax_dot.set_xlabel(FDR_LABEL, fontsize=15, color=GREY, labelpad=6, fontfamily=FONT)
ax_dot.tick_params(axis="x", labelsize=14, colors=GREY)
ax_dot.tick_params(axis="y", length=0)
for lbl in ax_dot.get_xticklabels():
    lbl.set_fontfamily(FONT)
ax_dot.spines[["top", "right", "left"]].set_visible(False)
ax_dot.spines["bottom"].set_color("#bbbbbb")

fig.suptitle(
    "GO BP - Constitutive EV Cargo\n"
    "(|ΔES| < 1 in both conditions,  ES_Normal > 1)",
    fontsize=22, fontweight="bold", color=BLUE_TITLE, y=0.99, fontfamily=FONT
)

LEGEND_KS = [25, 50, 75, 100]
legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(0.60),
           markersize=np.sqrt(k * DOT_SCALE), markeredgecolor='#888888',
           markeredgewidth=0.4, label=str(k))
    for k in LEGEND_KS
]
legend = ax_dot.legend(legend_handles, [str(k) for k in LEGEND_KS], title="Gene count",
                        loc="lower right", frameon=False, title_fontsize=14, fontsize=14,
                        handletextpad=1.5, labelspacing=2.8, borderpad=1.0)
legend.get_title().set_color(GREY)
legend.get_title().set_fontfamily(FONT)
for txt in legend.get_texts():
    txt.set_fontfamily(FONT)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0, vmax=vmax))
sm.set_array([])
cb = fig.colorbar(sm, ax=ax_cb, orientation="vertical", fraction=0.9, aspect=18)
cb.set_label(FDR_LABEL, fontsize=14, color=GREY, fontfamily=FONT)
cb.ax.tick_params(labelsize=12, colors=GREY)
for lbl in cb.ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
ax_cb.axis("off")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"GO_BP_constitutive.{ext}"), dpi=300, bbox_inches="tight")
plt.close()

print("Saved: GO_BP_constitutive.png/.svg")
print("Saved: GO_BP_constitutive_results.xlsx")
