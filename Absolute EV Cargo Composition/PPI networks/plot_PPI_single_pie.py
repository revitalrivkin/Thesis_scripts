"""
Single PPI network, colored circle nodes (pie-split for multi-category proteins), designed at
TRUE PRINT SIZE for the Word document (one network per figure)
==============================================================================
FINAL design (2026-09-26, Revital): square 6.5 x 6.5 in figure, network laid out in a circle, legend at
the bottom, no explanatory text in the figure (goes in the caption), paler node fills, bold 6.5 pt black
gene labels with a thin white outline centered on each circle. CIRCLE SIZE = NUMBER OF INTERACTIONS of the
protein within the network (radius = 6 + 1.7 * sqrt(degree) pt); labels of small circles may extend slightly
beyond them. Proteins in 2+ categories are drawn as equal wedges (one per category).

The canvas is the printed size in inches, layout is in figure pixels (DPI = 200,
axes = [0, 0, 1, 1]), so a 8 pt label here is 8 pt on the page.

Usage: python plot_PPI_single_pie.py   (builds the networks listed in NETWORKS below)
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, Wedge
import matplotlib.patheffects as pe

from plot_PPI_networks_pie import CATEGORY_COLORS, LEGEND_ORDER, HERE, lighten, node_categories

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

DPI = 200
PT = DPI / 72.0
FILL = 0.80            # paler node fill (standalone figures used 0.65)
LABEL_PT = 6.5
TITLE_PT = 13
LEGEND_PT = 10
FIG_W, FIG_H = 6.5, 6.5      # square page figure
TITLE_H, LEGEND_H = 0.42, 0.72   # inches reserved at the top / bottom

NETWORKS = [
    ("HL downregulated", "string_interactions_HL_downreg.tsv", "PPI Network of HL-Downregulated DEPs by Category",
     "PPI_HL_downregulated_print", {"deg_coef": 1.0, "spacing": 1.06}),
    ("HG downregulated", "string_interactions_HG_downreg.tsv", "PPI Network of HG-Downregulated DEPs by Category",
     "PPI_HG_downregulated_print"),
    ("HG upregulated", "string_interactions_HG_upreg.tsv", "PPI Network of HG-Upregulated DEPs by Category",
     "PPI_HG_upregulated_print"),
    ("HL upregulated", "string_interactions_HL_upreg.tsv", "PPI Network of HL-Upregulated DEPs by Category",
     "PPI_HL_upregulated_print", {"fig_h": 7.6, "deg_coef": 0.6, "spacing": 1.08}),
]


def draw_node(ax, xy, r, cats):
    if len(cats) == 1:
        c = CATEGORY_COLORS[cats[0]]
        ax.add_patch(Circle(xy, r, facecolor=lighten(c, FILL), edgecolor=c, linewidth=1.4, zorder=3))
        return
    step = 360.0 / len(cats)
    for i, cat in enumerate(cats):
        c = CATEGORY_COLORS[cat]
        ax.add_patch(Wedge(xy, r, 90 + i * step, 90 + (i + 1) * step, facecolor=lighten(c, FILL),
                           edgecolor=c, linewidth=1.4, zorder=3))


def build(folder, tsv, title, out_name, fig_h=FIG_H, deg_coef=1.7, spacing=1.06, save=True):
    edges = pd.read_csv(os.path.join(HERE, folder, tsv), sep="\t")
    edges.columns = [c.lstrip("#") for c in edges.columns]
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["node1"], row["node2"], weight=row["combined_score"])
    nodes = list(G.nodes())
    deg = dict(G.degree())
    cats = {n: node_categories(n) for n in nodes}

    fig = plt.figure(figsize=(FIG_W, fig_h), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    W, H = FIG_W * DPI, fig_h * DPI
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # circle radius from the label width (label sits inside), plus a bump for connectivity
    texts = {n: ax.text(0, 0, n, fontsize=LABEL_PT, fontweight="bold", fontfamily=FONT) for n in nodes}
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    R = np.zeros(len(nodes))      # footprint radius used by the de-overlap (px)
    RD = np.zeros(len(nodes))     # drawn circle radius, set by the number of interactions (px)
    for i, n in enumerate(nodes):
        bb = texts[n].get_window_extent(renderer=rend)
        RD[i] = (6.0 + deg_coef * np.sqrt(deg[n])) * PT
        R[i] = max(RD[i], 0.5 * bb.width * 0.95 + 0.6 * PT)
        texts[n].remove()

    # plotting region between the title (top) and the legend (bottom): a circle for a square
    # figure, an ellipse when the figure is taller than wide
    region_h = H - (TITLE_H + LEGEND_H) * DPI
    cx, cy = W / 2, LEGEND_H * DPI + region_h / 2
    ax_a, ax_b = W / 2 - 0.04 * DPI, region_h / 2
    if ax_b <= ax_a:
        ax_a = ax_b                      # circle
    pos_u = nx.spring_layout(G, weight="weight", seed=42, k=2.2, iterations=200)
    U = np.array([pos_u[n] for n in nodes])
    U = U - U.mean(axis=0)
    if ax_a == ax_b:
        U = U / np.linalg.norm(U, axis=1).max()
    else:
        U = U / np.abs(U).max(axis=0)
    pts = np.column_stack([cx + U[:, 0] * (ax_a - R.max()), cy + U[:, 1] * (ax_b - R.max())])

    N = len(nodes)
    rng = np.random.default_rng(42)
    for it in range(1500):
        diff = pts[:, None, :] - pts[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(dist, np.inf)
        zero = dist < 1e-6
        if zero.any():
            jit = rng.uniform(-1, 1, size=(N, N, 2))
            diff = np.where(zero[:, :, None], jit, diff)
            dist = np.where(zero, np.linalg.norm(jit, axis=2), dist)
        min_dist = (R[:, None] + R[None, :] + 2) * spacing
        overlap = np.clip(min_dist - dist, 0, None)
        if overlap.max() < 0.5:
            break
        pts += ((overlap[:, :, None] / 2) * diff / dist[:, :, None]).sum(axis=1) * 0.5
        rel = pts - np.array([cx, cy])
        e = np.sqrt((rel[:, 0] / (ax_a - R)) ** 2 + (rel[:, 1] / (ax_b - R)) ** 2)
        over = e > 1
        pts[over] = np.array([cx, cy]) + rel[over] / e[over][:, None]
    n_overlap = int(((overlap > 0.5).sum()) // 2)
    print(f"{folder}: {N} nodes, {G.number_of_edges()} edges, de-overlap {it + 1} its, "
          f"{n_overlap} overlapping pairs left (max {overlap.max():.1f}px)")

    pos = {n: pts[i] for i, n in enumerate(nodes)}
    for u, v in G.edges():
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#c8c8c8", alpha=0.6,
                linewidth=0.2 + 0.6 * G[u][v]["weight"], zorder=1, solid_capstyle="round")
    for i, n in enumerate(nodes):
        draw_node(ax, pos[n], RD[i], cats[n])
    for n in nodes:
        ax.text(pos[n][0], pos[n][1], n, fontsize=LABEL_PT, fontweight="bold", fontfamily=FONT,
                ha="center", va="center", zorder=6, color="#000000",
                path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])

    ax.text(W / 2, H - 0.2 * DPI, title, fontsize=TITLE_PT, fontweight="bold", fontfamily=FONT,
            ha="center", va="center", color="#1f3864")
    per_row = 3
    col_w = W / per_row
    for i, c in enumerate(LEGEND_ORDER):
        r_, k = divmod(i, per_row)
        lx = 0.2 * DPI + k * col_w
        ly = (LEGEND_H - 0.24 - 0.27 * r_) * DPI
        col = CATEGORY_COLORS[c]
        ax.add_patch(Circle((lx, ly), 0.075 * DPI, facecolor=lighten(col, FILL), edgecolor=col, linewidth=1.4))
        ax.text(lx + 0.16 * DPI, ly, c, fontsize=LEGEND_PT, fontfamily=FONT, va="center", ha="left")

    for ext in (("png", "svg") if save else ()):
        path = os.path.join(HERE, folder, f"{out_name}.{ext}")
        fig.savefig(path, dpi=DPI)
        print("  Saved:", path)
    plt.close(fig)


if __name__ == "__main__":
    for cfg in NETWORKS:
        opts = cfg[4] if len(cfg) > 4 else {}
        build(*cfg[:4], **opts)
