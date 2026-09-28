"""
PPI Network — HG-Downregulated DEPs, colored by category of interest
==============================================================================
Built directly from STRING's TSV interaction export (networkx + matplotlib),
skipping Cytoscape, so the figure can match thesis house style directly.
See [[12 - Thesis Writing]] for the full list of non-default choices
(edge meaning = Confidence, cutoff kept at STRING default 0.4, isolated
nodes excluded).

Categories now use the finalized GO-derived gene sets (all 5 categories,
same source and priority-resolution order as the category scatter plot —
Absolute EV Cargo Composition/Volcano/Volcano with categories/
02_plot_volcano_categories_HG.py) instead of the earlier hand-picked
marker lists.

Input:  string_interactions_HG_downreg.tsv (STRING export, edges = Confidence,
        cutoff 0.4 — as configured by Revital)
        ../../Volcano/Volcano with categories/category_genesets.csv
Output: (manual, via the toolbar's save button) PPI_HG_downregulated.png / .svg

Interactive: opens a matplotlib window (plt.show()) — drag any node to
reposition it (edges and its label follow); the force-directed layout is
just a starting point, not a fixed result. Resize the window as needed,
then use the toolbar's save icon (floppy disk) to export. RUN THIS LOCALLY,
not headless.
"""

import os
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# NOTE: no matplotlib.use("Agg") here — this script opens an INTERACTIVE window.
# Run it locally to drag nodes, then use the window's save icon to export.


def lighten(hex_color, amount=0.65):
    """Blend a hex color toward white by `amount` (0=no change, 1=white)."""
    r, g, b = mcolors.to_rgb(hex_color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)

FONT = "Arial"
plt.rcParams["font.family"] = FONT
plt.rcParams["font.sans-serif"] = [FONT]

HERE = os.path.dirname(os.path.abspath(__file__))
EDGES_FILE = os.path.join(HERE, "string_interactions_HG_downreg.tsv")
GENESETS_FILE = os.path.join(HERE, "..", "..", "Protein Categories of Interest", "category_genesets.csv")

# Same priority order as the category scatter plot (most specific first, so
# rare categories aren't swamped by the much larger Adhesion/Docking/Uptake set)
CATEGORY_PRIORITY = [
    "EV Biogenesis",
    "ECM Organization",
    "Immune-evasion",
    "Cardioprotective",
    "Adhesion / Docking / Uptake",
]
CATEGORY_COLORS = {
    "Adhesion / Docking / Uptake": "#1f77b4",
    "EV Biogenesis": "#9467bd",
    "Immune-evasion": "#2ca02c",
    "ECM Organization": "#8c564b",
    "Cardioprotective": "#d62728",
    "Other DEP": "#bbbbbb",
}

# -- Load edges, build graph (isolated nodes excluded automatically) -----------
edges = pd.read_csv(EDGES_FILE, sep="\t")
edges.columns = [c.lstrip("#") for c in edges.columns]

G = nx.Graph()
for _, row in edges.iterrows():
    G.add_edge(row["node1"], row["node2"], weight=row["combined_score"])

print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

# -- Category assignment (same GO-derived sets + priority order as the scatter) --
genesets = pd.read_csv(GENESETS_FILE)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()

def categorize(gene):
    cats = gene_to_cats.get(gene, set())
    for cat in CATEGORY_PRIORITY:
        if cat in cats:
            return cat
    return "Other DEP"

categories = {n: categorize(n) for n in G.nodes()}
for cat in list(CATEGORY_PRIORITY) + ["Other DEP"]:
    n_cat = sum(1 for c in categories.values() if c == cat)
    print(f"  {cat}: {n_cat} nodes")

# -- Layout (force-directed, weighted by STRING confidence) --------------------
# k pushed further (1.3 -> 2.2) and iterations increased, combined with the
# larger square axes below, so nodes start with much less overlap and need
# only minor manual adjustment
pos = nx.spring_layout(G, weight="weight", seed=42, k=2.2, iterations=200)

# -- Draw -------------------------------------------------------------------------
# Large square axes filling most of the figure (per Revital's marked layout),
# with only a narrow reserved column on the left for the legend and a margin
# at top for the title — neither can overlap the network regardless of how
# nodes are dragged.
fig = plt.figure(figsize=(16, 16))
ax = fig.add_axes([0.24, 0.06, 0.74, 0.80])
ax.set_aspect("equal", adjustable="box")

edge_list = list(G.edges())  # fixed order, matches the LineCollection's segments
edge_widths = [G[u][v]["weight"] * 3 for u, v in G.edges()]
edge_artist = nx.draw_networkx_edges(G, pos, ax=ax, edgelist=edge_list, width=edge_widths,
                                      edge_color="#cccccc", alpha=0.7)

degrees = dict(G.degree())
node_collections = {}  # category -> (PathCollection, [nodes in draw order])
for cat, color in CATEGORY_COLORS.items():
    nodes_in_cat = [n for n in G.nodes() if categories[n] == cat]
    if not nodes_in_cat:
        continue
    sizes = [200 + degrees[n] * 60 for n in nodes_in_cat]
    coll = nx.draw_networkx_nodes(G, pos, nodelist=nodes_in_cat, ax=ax, node_color=[lighten(color)],
                                   node_size=sizes, edgecolors=color, linewidths=1.8, alpha=0.95)
    coll.set_picker(5)
    node_collections[cat] = (coll, nodes_in_cat)

label_objects = nx.draw_networkx_labels(G, pos, ax=ax, font_size=11, font_family=FONT, font_weight="bold")

ax.axis("off")

fig.suptitle("PPI Network of HG-Downregulated DEPs by Category", fontsize=19, fontweight="bold",
             fontfamily=FONT, y=0.97)

legend_handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=lighten(color),
                              markersize=12, markeredgecolor=color, markeredgewidth=1.8, label=cat)
                   for cat, color in CATEGORY_COLORS.items() if any(c == cat for c in categories.values())]
# Legend placed in the figure's reserved left column, entirely outside the
# square network axes (which starts at x=0.24) — no overlap regardless of
# how nodes are dragged.
fig.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(0.01, 0.68),
           frameon=False, fontsize=12, prop={"family": FONT}, labelspacing=1.4, borderpad=1.0)

# -- Shared helpers (used by both the de-overlap pass below and interactive drag) --
node_lookup = {}
for cat, (coll, nodes_in_cat) in node_collections.items():
    for i, n in enumerate(nodes_in_cat):
        node_lookup[n] = (cat, i)

def redraw_node(node):
    cat, i = node_lookup[node]
    coll, nodes_in_cat = node_collections[cat]
    offsets = coll.get_offsets()
    offsets[i] = pos[node]
    coll.set_offsets(offsets)
    label_objects[node].set_position(pos[node])

def redraw_edges():
    segments = [(pos[u], pos[v]) for u, v in edge_list]
    edge_artist.set_segments(segments)

# -- Automatic de-overlap pass ------------------------------------------------------
# Spring layout alone still leaves the densely-interconnected hub cluster
# overlapping. This pass works entirely in screen PIXELS (matching what
# ax.transData.transform gives, and what get_window_extent reports natively
# — no unit mixing) and iteratively pushes any two nodes whose circle+label
# footprints collide apart, so the starting layout needs only minor manual
# touch-up rather than a full rearrangement.
fig.canvas.draw()  # realize the layout once so transforms/label extents are valid
renderer = fig.canvas.get_renderer()
px_per_point = fig.dpi / 72.0

node_list = list(G.nodes())
radius_px = {}
for n in node_list:
    node_size = 200 + degrees[n] * 60
    circle_r_px = 0.5 * (node_size ** 0.5) * px_per_point  # node_size is in points^2
    bbox = label_objects[n].get_window_extent(renderer=renderer)  # already in pixels
    label_r_px = 0.5 * (bbox.width ** 2 + bbox.height ** 2) ** 0.5
    radius_px[n] = max(circle_r_px, label_r_px) + 4  # small pixel buffer

initial_pts = np.array([ax.transData.transform(pos[n]) for n in node_list], dtype=float)
radii_arr = np.array([radius_px[n] for n in node_list])
N = len(node_list)
rng = np.random.default_rng(42)

positions_arr = initial_pts.copy()
for iteration in range(500):
    diff = positions_arr[:, None, :] - positions_arr[None, :, :]  # (N, N, 2)
    dist = np.linalg.norm(diff, axis=2)
    np.fill_diagonal(dist, np.inf)
    zero_mask = dist < 1e-6
    if zero_mask.any():
        jitter = rng.uniform(-1, 1, size=(N, N, 2))
        diff = np.where(zero_mask[:, :, None], jitter, diff)
        dist = np.where(zero_mask, np.linalg.norm(jitter, axis=2), dist)
    min_dist = (radii_arr[:, None] + radii_arr[None, :]) * 1.08
    overlap = np.clip(min_dist - dist, 0, None)
    if overlap.max() < 0.5:  # converged to sub-pixel precision
        break
    direction = diff / dist[:, :, None]
    disp = (overlap[:, :, None] / 2 * direction).sum(axis=1)
    positions_arr += disp * 0.5  # damping to avoid oscillation

print(f"De-overlap pass: stopped after {iteration + 1} iterations "
      f"(max remaining overlap: {overlap.max():.1f}px)")

inv = ax.transData.inverted()
for i, n in enumerate(node_list):
    pos[n] = tuple(inv.transform(positions_arr[i]))
    redraw_node(n)
redraw_edges()
fig.canvas.draw_idle()

# -- Interactive node dragging ----------------------------------------------------
xs = np.array([p[0] for p in pos.values()])
ys = np.array([p[1] for p in pos.values()])
click_radius = 0.03 * max(xs.max() - xs.min(), ys.max() - ys.min())

dragging = {"node": None}

def find_nearest_node(x, y):
    best, best_dist = None, click_radius
    for n, (nx_, ny_) in pos.items():
        dist = ((nx_ - x) ** 2 + (ny_ - y) ** 2) ** 0.5
        if dist < best_dist:
            best, best_dist = n, dist
    return best

def on_press(event):
    if event.inaxes != ax or event.xdata is None:
        return
    dragging["node"] = find_nearest_node(event.xdata, event.ydata)

def on_motion(event):
    if dragging["node"] is None or event.inaxes != ax or event.xdata is None:
        return
    pos[dragging["node"]] = (event.xdata, event.ydata)
    redraw_node(dragging["node"])
    redraw_edges()
    fig.canvas.draw_idle()

def on_release(event):
    dragging["node"] = None

fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("motion_notify_event", on_motion)
fig.canvas.mpl_connect("button_release_event", on_release)

print("Drag any node to reposition it. Use the toolbar's save icon to export once done.")
plt.show()
