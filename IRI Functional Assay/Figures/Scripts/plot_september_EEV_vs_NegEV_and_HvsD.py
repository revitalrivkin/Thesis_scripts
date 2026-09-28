"""
4.1.3 IRI Functional Assay Results - September: EEV vs. NegEV, Healthy vs. Diabetic EEV
================================================================================
Two-panel figure, September 2026 IRI experiment, mirroring July's structure
exactly (plot_july_EEV_vs_NegEV_and_HvsD.py):
  (A) IRI control / EEV / NegEV fold rescue, all donors pooled by type.
  (B) IRI control / Healthy EEV / Diabetic EEV fold rescue.

Same 3-group/3-bracket-per-panel design, same uncorrected pairwise
Mann-Whitney, same p<0.05/"p < 0.001"/"ns" labeling, same bar-color roles
(Control=pale grey, NegEV=darker grey, EEV=green, Healthy=blue,
Diabetic=red) as July.

Structural exclusions (see plot_september_raw_overview.py and note 12's
3.8 candidate material for full detail): Diabetic donors 3T and P4 excluded
entirely (pipetting error); usable Diabetic donors are P5, P6, P7, 1T only.
Diabetic IRI control row A excluded (pipetting/edge error). Healthy plate's
own IRI control not used (IRI had no effect there) - the diabetic plate's
control (rows B-H, n=7, mean 7.1957M) is the shared reference for both
plates. Healthy H10 (col3) and H1 (col11) each lose their row-H well
(dropped, not reallocated - Revital's call, 2026-09-21). All other
exclusions: generic Tukey's-fence per donor column.

Donor palette: EXTENDS July's per-donor-shade-within-group-family palette
rather than replacing it. H10 and P5 (the two donors appearing in both
experiments alongside 3T, which isn't usable in September) keep their
EXACT July colors, for visual continuity across the two experiments'
figures. New donors (H11, H2, H3, H4, H1 - Healthy; P6, P7, 1T - Diabetic)
get new shades slotted into the same blue/warm-family progressions.

Data source: ../../Data/September/06.09.26 Viabilty Healthy plate.xlsx,
             ../../Data/September/06.09.26 Viabilty Diabetic plate.xlsx

Output: September_EEV_vs_NegEV_and_HvsD.png/.svg (one level up, in Figures/)
"""

import os
import numpy as np
import openpyxl
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT = "Arial"
BG = "#fcfcfb"
GRID = "#e1e0d9"
SPINE = "#c3c2b7"
TITLE_COLOR = "#0b0b0b"
LABEL_COLOR = "#52514e"

C_CONTROL = "#c7cbd1"
C_CONTROL_DOT = "#9aa0a8"
C_NEG = "#7d838c"
C_EEV = "#3ba272"
C_HEALTHY, C_DIABETIC = "#2a78d6", "#e34948"

# Donor colors: H10 and D5 MUST match July's plot_july_EEV_vs_NegEV_and_HvsD.py
# exactly (same donors, cross-experiment visual continuity). New donors get
# new shades within the same blue (Healthy) / warm (Diabetic) families,
# ordered light-to-dark alongside the existing anchors. Diabetic codes use
# the unified D-naming convention (1T->D1, 2T->D2, 3T->D3, P5->D5, P6->D6,
# P7->D7) - Healthy stays H.
DONOR_COLORS = {
    "H8": "#123c69", "H11": "#1c5cab", "H9": "#2a78d6", "H2": "#4a9de8",
    "H4": "#6fb8f0", "H10": "#8ec4f0", "H3": "#a8d5ff", "H1": "#c3e4ff",
    "D2": "#7a1f2b", "D1": "#a83246", "D3": "#e34948", "D6": "#eb6b3d",
    "D7": "#f2914f", "D5": "#e8a33d",
}
DONOR_ORDER = ["H11", "H10", "H2", "H3", "H4", "H1", "D5", "D6", "D7", "D1"]
HEALTHY_DONORS = ["H11", "H10", "H2", "H3", "H4", "H1"]
DIABETIC_DONORS = ["D5", "D6", "D7", "D1"]

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
DATA_DIR = os.path.join(FIG_DIR, "..", "Data", "September")
XLSX_HEALTHY = os.path.join(DATA_DIR, "06.09.26 Viabilty Healthy plate.xlsx")
XLSX_DIABETIC = os.path.join(DATA_DIR, "06.09.26 Viabilty Diabetic plate.xlsx")

RNG_SEED = 42
ROWS = list("ABCDEFGH")


def read_plate(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Luminescence 1_01"]
    raw = {}
    for ri, r in enumerate(ROWS):
        for col in range(1, 13):
            raw[(r, col)] = ws.cell(row=11 + ri, column=col + 1).value
    return raw


def get_wells_M(raw, col, drop_rows=()):
    vals = []
    for r in ROWS:
        if r in drop_rows:
            continue
        v = raw.get((r, col))
        if v is not None:
            vals.append(v / 1e6)
    return np.array(vals)


def tukey_filter(vals):
    a = np.asarray(vals, dtype=float)
    if len(a) < 4:
        return a
    q1, q3 = np.percentile(a, 25), np.percentile(a, 75)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return a[(a >= lo) & (a <= hi)]


healthy_raw = read_plate(XLSX_HEALTHY)
diabetic_raw = read_plate(XLSX_DIABETIC)

iri_ctrl_raw = get_wells_M(diabetic_raw, 12, drop_rows={"A"})
iri_ctrl_vals = tukey_filter(iri_ctrl_raw)
IRI_CTRL_MEAN = iri_ctrl_vals.mean()
iri_ctrl_fold = iri_ctrl_vals / IRI_CTRL_MEAN

DONOR_COLS = {
    "H11": {"plate": healthy_raw, "EEV": 1, "NegEV": 2},
    "H10": {"plate": healthy_raw, "EEV": (3, {"H"}), "NegEV": 4},
    "H2": {"plate": healthy_raw, "EEV": 5, "NegEV": 6},
    "H3": {"plate": healthy_raw, "EEV": 7, "NegEV": 8},
    "H4": {"plate": healthy_raw, "EEV": 9, "NegEV": 10},
    "H1": {"plate": healthy_raw, "EEV": (11, {"H"})},
    "D5": {"plate": diabetic_raw, "EEV": 5, "NegEV": 6},
    "D6": {"plate": diabetic_raw, "EEV": 7, "NegEV": 8},
    "D7": {"plate": diabetic_raw, "EEV": 9, "NegEV": 10},
    "D1": {"plate": diabetic_raw, "EEV": 11},
}

donor_data = {}
for donor, spec in DONOR_COLS.items():
    donor_data[donor] = {}
    plate = spec["plate"]
    for typ in ("EEV", "NegEV"):
        if typ not in spec:
            continue
        col_spec = spec[typ]
        if isinstance(col_spec, tuple):
            col, drop_rows = col_spec
        else:
            col, drop_rows = col_spec, set()
        raw_vals = get_wells_M(plate, col, drop_rows=drop_rows)
        donor_data[donor][typ] = tukey_filter(raw_vals) / IRI_CTRL_MEAN


def pooled(typ, donors):
    return np.concatenate([donor_data[d][typ] for d in donors if typ in donor_data[d]])


eev_all = pooled("EEV", DONOR_ORDER)
negev_all = pooled("NegEV", DONOR_ORDER)
eev_healthy = pooled("EEV", HEALTHY_DONORS)
eev_diabetic = pooled("EEV", DIABETIC_DONORS)


def sig_label(p):
    if p >= 0.05:
        return "ns"
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


def add_panel_letter(ax, letter):
    ax.text(-0.10, 1.10, letter, transform=ax.transAxes, fontsize=16,
            fontweight="bold", color=TITLE_COLOR, fontfamily=FONT,
            va="bottom", ha="left")


def sig_bracket(ax, x1, x2, y, label):
    tick = 0.10
    ax.plot([x1, x1], [y - tick, y], color=LABEL_COLOR, linewidth=1.2, zorder=6)
    ax.plot([x2, x2], [y - tick, y], color=LABEL_COLOR, linewidth=1.2, zorder=6)
    ax.plot([x1, x2], [y, y], color=LABEL_COLOR, linewidth=1.2, zorder=6)
    ax.text((x1 + x2) / 2, y + 0.05, label, ha="center", va="bottom",
             fontsize=10.5, fontfamily=FONT, color=LABEL_COLOR)


def style_axes(ax):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.tick_params(axis="y", labelsize=12, labelcolor=LABEL_COLOR)
    ax.tick_params(axis="x", bottom=False)
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(SPINE)
    ax.spines["left"].set_linewidth(0.8)


def plot_panel(ax, groups, bar_colors, letter, title):
    ax.set_facecolor(BG)
    rng = np.random.default_rng(RNG_SEED)
    x = np.arange(3)
    means = [g["pool"].mean() for g in groups]
    sds = [g["pool"].std(ddof=1) for g in groups]
    ns = [len(g["pool"]) for g in groups]

    ax.bar(x, means, yerr=sds, color=bar_colors, edgecolor=bar_colors,
           linewidth=0.8, error_kw=dict(elinewidth=1.4, capsize=5, capthick=1.4,
           ecolor=LABEL_COLOR), width=0.55, zorder=3, alpha=0.65)

    for xi, g in zip(x, groups):
        if g["by_donor"] is None:
            vals = g["pool"]
            jit = rng.uniform(-0.16, 0.16, len(vals))
            ax.scatter(np.full(len(vals), xi) + jit, vals, color=C_CONTROL_DOT,
                       edgecolors="white", linewidths=0.6, s=48, zorder=5)
        else:
            for donor in DONOR_ORDER:
                if donor not in g["by_donor"]:
                    continue
                vals = g["by_donor"][donor]
                if len(vals) == 0:
                    continue
                jit = rng.uniform(-0.16, 0.16, len(vals))
                ax.scatter(np.full(len(vals), xi) + jit, vals, color=DONOR_COLORS[donor],
                           edgecolors="white", linewidths=0.6, s=48, zorder=5)

    ax.axhline(1.0, color=LABEL_COLOR, linestyle=(0, (4, 3)), linewidth=1.1, zorder=2)

    data_top = max(max(g["pool"]) for g in groups)
    bar_top = max(m + s for m, s in zip(means, sds))
    base = max(data_top, bar_top)

    p01 = stats.mannwhitneyu(groups[0]["pool"], groups[1]["pool"], alternative="two-sided").pvalue
    p12 = stats.mannwhitneyu(groups[1]["pool"], groups[2]["pool"], alternative="two-sided").pvalue
    p02 = stats.mannwhitneyu(groups[0]["pool"], groups[2]["pool"], alternative="two-sided").pvalue

    h01, h12, h02 = base * 1.15, base * 1.30, base * 1.45
    sig_bracket(ax, x[0], x[1], h01, sig_label(p01))
    sig_bracket(ax, x[1], x[2], h12, sig_label(p12))
    sig_bracket(ax, x[0], x[2], h02, sig_label(p02))

    print(f"  {groups[0]['label']} vs {groups[1]['label']}: p={p01:.5f}")
    print(f"  {groups[1]['label']} vs {groups[2]['label']}: p={p12:.5f}")
    print(f"  {groups[0]['label']} vs {groups[2]['label']}: p={p02:.5f}")

    ax.set_ylim(0, h02 * 1.12)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g['label']}\n(n={n})" for g, n in zip(groups, ns)],
                        fontsize=12, fontfamily=FONT, color=TITLE_COLOR)
    ax.set_xlim(-0.6, 2.6)
    style_axes(ax)
    ax.set_title(title, fontsize=13.5, fontweight="bold", fontfamily=FONT,
                 color=TITLE_COLOR, pad=14, loc="center")
    add_panel_letter(ax, letter)


fig, axes = plt.subplots(1, 2, figsize=(14, 7.5))
fig.patch.set_facecolor(BG)

print("Panel A (IRI ctrl / EEV / NegEV):")
panel_a_groups = [
    dict(label="IRI Control", pool=iri_ctrl_fold, by_donor=None),
    dict(label="EEV", pool=eev_all, by_donor={d: donor_data[d]["EEV"] for d in DONOR_ORDER}),
    dict(label="NegEV", pool=negev_all, by_donor={d: donor_data[d]["NegEV"] for d in DONOR_ORDER if "NegEV" in donor_data[d]}),
]
plot_panel(axes[0], panel_a_groups, [C_CONTROL, C_EEV, C_NEG], "A", "EEV vs. NegEV Fold Rescue - September")

print("Panel B (IRI ctrl / Healthy EEV / Diabetic EEV):")
panel_b_groups = [
    dict(label="IRI Control", pool=iri_ctrl_fold, by_donor=None),
    dict(label="Healthy EEV", pool=eev_healthy, by_donor={d: donor_data[d]["EEV"] for d in HEALTHY_DONORS}),
    dict(label="Diabetic EEV", pool=eev_diabetic, by_donor={d: donor_data[d]["EEV"] for d in DIABETIC_DONORS}),
]
plot_panel(axes[1], panel_b_groups, [C_CONTROL, C_HEALTHY, C_DIABETIC], "B", "EEV Fold Rescue - Healthy vs. Diabetic - September")

for ax in axes:
    ax.set_ylabel("Fold Rescue (relative to IRI control)", fontsize=12.5,
                  fontfamily=FONT, color=LABEL_COLOR, labelpad=10)

legend_handles = [plt.scatter([], [], color=DONOR_COLORS[d], edgecolors="white",
                               linewidths=0.6, s=60, label=d) for d in DONOR_ORDER]
_leg = fig.legend(handles=legend_handles, loc="lower center", ncol=10, fontsize=11,
                  frameon=False, bbox_to_anchor=(0.5, -0.06), handletextpad=0.5,
                  columnspacing=1.2, title="Donor", title_fontsize=12,
                  prop={"family": FONT})
_leg.get_title().set_fontfamily(FONT)

fig.suptitle("September 2026 IRI Experiment - Functional Rescue", fontsize=17,
             fontweight="bold", color=TITLE_COLOR, y=1.02, fontfamily=FONT)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"September_EEV_vs_NegEV_and_HvsD.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print("\nSaved: September_EEV_vs_NegEV_and_HvsD.png/.svg")
