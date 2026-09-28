"""
PPI networks of DEPs, nodes split like a pie chart when a protein belongs to 2+ categories
==============================================================================
Static (headless) rebuild of the interactive plot_PPI_*.py scripts. Same design
choices (STRING TSV edges with combined_score, cutoff 0.4 as exported; isolated
nodes excluded; node size = 200 + 60 * degree; light fill with a full-strength
outline; spring layout seeded by STRING confidence; automatic de-overlap pass),
with ONE change: a protein that belongs to several of the 5 categories is drawn
as a pie with one equal wedge per category (in CATEGORY_ORDER), so the
category priority rule is no longer needed for these networks. Proteins in no
category are grey ("Other DEP").

Everything is laid out in figure PIXELS (axes = [0, 0, 1, 1], 1 data unit = 1 px
at the base dpi), so circles stay circular and the de-overlap works in the same
units it draws in.

Inputs:  <condition folder>/string_interactions_<...>.tsv (STRING export)
         ../Protein Categories of Interest/category_genesets.csv
Outputs: <condition folder>/PPI_<...>_pie.png (+ .svg)
"""

import os

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Wedge, Circle

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
GENESETS_FILE = os.path.join(HERE, "..", "Protein Categories of Interest", "category_genesets.csv")

CATEGORY_ORDER = [
    "Adhesion / Docking / Uptake",
    "Cardioprotective",
    "ECM Organization",
    "Immune-evasion",
    "EV Biogenesis",
]
CATEGORY_COLORS = {
    "Adhesion / Docking / Uptake": "#1f77b4",
    "EV Biogenesis": "#9467bd",
    "Immune-evasion": "#2ca02c",
    "ECM Organization": "#8c564b",
    "Cardioprotective": "#d62728",
    "Other DEP": "#bbbbbb",
}
LEGEND_ORDER = ["Adhesion / Docking / Uptake", "EV Biogenesis", "Immune-evasion", "ECM Organization",
                "Cardioprotective", "Other DEP"]

NETWORKS = [
    ("HG downregulated", "string_interactions_HG_downreg.tsv", "PPI Network of HG-Downregulated DEPs by Category",
     "PPI_HG_downregulated_pie"),
    ("HL downregulated", "string_interactions_HL_downreg.tsv", "PPI Network of HL-Downregulated DEPs by Category",
     "PPI_HL_downregulated_pie"),
    ("HG upregulated", "string_interactions_HG_upreg.tsv", "PPI Network of HG-Upregulated DEPs by Category",
     "PPI_HG_upregulated_pie"),
    ("HL upregulated", "string_interactions_HL_upreg.tsv", "PPI Network of HL-Upregulated DEPs by Category",
     "PPI_HL_upregulated_pie"),
]

DPI = 100
FIG_IN = 16
PX = FIG_IN * DPI          # 1600 px square
LEGEND_W = 340             # reserved left column (px)
PT_TO_PX = DPI / 72.0


def lighten(hex_color, amount=0.65):
    r, g, b = mcolors.to_rgb(hex_color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


genesets = pd.read_csv(GENESETS_FILE)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()


def node_categories(gene):
    cats = [c for c in CATEGORY_ORDER if c in gene_to_cats.get(gene, set())]
    return cats or ["Other DEP"]


def draw_pie(ax, xy, radius, cats, z=3):
    """One equal wedge per category, light fill + full-strength outline in each wedge's color."""
    if len(cats) == 1:
        c = CATEGORY_COLORS[cats[0]]
        ax.add_patch(Circle(xy, radius, facecolor=lighten(c), edgecolor=c, linewidth=1.8, zorder=z))
        return
    step = 360.0 / len(cats)
    for i, cat in enumerate(cats):
        c = CATEGORY_COLORS[cat]
        ax.add_patch(Wedge(xy, radius, 90 + i * step, 90 + (i + 1) * step, facecolor=lighten(c),
                           edgecolor=c, linewidth=1.8, zorder=z))


def build(folder, tsv, title, out_name):
    edges = pd.read_csv(os.path.join(HERE, folder, tsv), sep="\t")
    edges.columns = [c.lstrip("#") for c in edges.columns]
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["node1"], row["node2"], weight=row["combined_score"])
    print(f"{folder}: nodes {G.number_of_nodes()}, edges {G.number_of_edges()}")

    cats = {n: node_categories(n) for n in G.nodes()}
    n_multi = sum(len(v) > 1 for v in cats.values())
    print(f"  proteins in 2+ categories: {n_multi}")
    for c in LEGEND_ORDER:
        print(f"  {c}: {sum(c in v for v in cats.values())} nodes")

    # layout in unit space, then mapped into the pixel box right of the legend column
    pos_u = nx.spring_layout(G, weight="weight", seed=42, k=2.2, iterations=200)
    nodes = list(G.nodes())
    U = np.array([pos_u[n] for n in nodes])
    lo, hi = U.min(axis=0), U.max(axis=0)
    box_x0, box_x1, box_y0, box_y1 = LEGEND_W + 30, PX - 40, 60, PX - 110
    scale = min((box_x1 - box_x0) / (hi[0] - lo[0]), (box_y1 - box_y0) / (hi[1] - lo[1]))
    P = (U - lo) * scale
    P[:, 0] += box_x0 + ((box_x1 - box_x0) - (hi[0] - lo[0]) * scale) / 2
    P[:, 1] += box_y0 + ((box_y1 - box_y0) - (hi[1] - lo[1]) * scale) / 2

    degrees = dict(G.degree())
    fig = plt.figure(figsize=(FIG_IN, FIG_IN), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, PX)
    ax.set_ylim(0, PX)
    ax.axis("off")

    # label extents (px) from a throwaway draw
    label_objs = {n: ax.text(0, 0, n, fontsize=11, fontweight="bold", fontfamily=FONT) for n in nodes}
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    radius = {}
    node_r = {}
    for n in nodes:
        node_size = 200 + degrees[n] * 60
        node_r[n] = 0.5 * (node_size ** 0.5) * PT_TO_PX
        bb = label_objs[n].get_window_extent(renderer=renderer)
        label_r = 0.5 * (bb.width ** 2 + bb.height ** 2) ** 0.5
        radius[n] = max(node_r[n], label_r) + 4
    for t in label_objs.values():
        t.remove()

    # automatic de-overlap pass (vectorized), same as the interactive scripts
    N = len(nodes)
    pts = P.copy()
    radii = np.array([radius[n] for n in nodes])
    rng = np.random.default_rng(42)
    for it in range(500):
        diff = pts[:, None, :] - pts[None, :, :]
        dist = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(dist, np.inf)
        zero = dist < 1e-6
        if zero.any():
            jit = rng.uniform(-1, 1, size=(N, N, 2))
            diff = np.where(zero[:, :, None], jit, diff)
            dist = np.where(zero, np.linalg.norm(jit, axis=2), dist)
        min_dist = (radii[:, None] + radii[None, :]) * 1.08
        overlap = np.clip(min_dist - dist, 0, None)
        if overlap.max() < 0.5:
            break
        direction = diff / dist[:, :, None]
        pts += (overlap[:, :, None] / 2 * direction).sum(axis=1) * 0.5
    print(f"  de-overlap: {it + 1} iterations, max remaining overlap {overlap.max():.1f}px")

    # if the pass pushed nodes outside the canvas, shrink the whole picture toward the center
    x_lo, x_hi = pts[:, 0].min() - 30, pts[:, 0].max() + 30
    y_lo, y_hi = pts[:, 1].min() - 30, pts[:, 1].max() + 30
    ax.set_xlim(min(0, x_lo), max(PX, x_hi))
    ax.set_ylim(min(0, y_lo), max(PX, y_hi))
    ax.set_aspect("equal", adjustable="box")

    pos = {n: pts[i] for i, n in enumerate(nodes)}
    for u, v in G.edges():
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#cccccc", alpha=0.7,
                linewidth=G[u][v]["weight"] * 3 * PT_TO_PX / PT_TO_PX, zorder=1, solid_capstyle="round")
    for n in nodes:
        draw_pie(ax, pos[n], node_r[n], cats[n])
    for n in nodes:
        ax.text(pos[n][0], pos[n][1], n, fontsize=11, fontweight="bold", fontfamily=FONT,
                ha="center", va="center", zorder=5)

    # title + legend (pixel coordinates of the base canvas)
    ax.text(PX / 2, PX - 45, title, fontsize=19, fontweight="bold", fontfamily=FONT, ha="center", va="center")
    ly = PX * 0.66
    for c in LEGEND_ORDER:
        if not any(c in v for v in cats.values()):
            continue
        col = CATEGORY_COLORS[c]
        ax.add_patch(Circle((40, ly), 9, facecolor=lighten(col), edgecolor=col, linewidth=1.8, zorder=6))
        ax.text(62, ly, c, fontsize=12, fontfamily=FONT, va="center", ha="left")
        ly -= 38
    if n_multi:
        ax.text(30, ly - 10, "Nodes split into wedges belong\nto more than one category", fontsize=11,
                fontfamily=FONT, va="top", ha="left", color="#444444")

    for ext in ("png", "svg"):
        path = os.path.join(HERE, folder, f"{out_name}.{ext}")
        fig.savefig(path, dpi=200 if ext == "png" else DPI)
        print("  Saved:", path)
    plt.close(fig)


if __name__ == "__main__":
    for cfg in NETWORKS:
        build(*cfg)
