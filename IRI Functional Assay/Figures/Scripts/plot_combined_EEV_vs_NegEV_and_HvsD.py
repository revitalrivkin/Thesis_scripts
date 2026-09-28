"""
4.1.3 IRI Functional Assay Results - Combined (July + September): EEV vs. NegEV,
Healthy vs. Diabetic EEV
================================================================================
THREE-panel figure pooling well-level fold-rescue data from BOTH the July
and September IRI experiments, per [[11.3 - IRI Combined Analysis (July +
September 2026)]]:
  (A) EEV vs. NegEV, all donors pooled.
  (B) Healthy EEV vs. Diabetic EEV.
  (C) Healthy EEV / Healthy NegEV / Diabetic EEV / Diabetic NegEV, grouped
      by disease group so each within-group EEV-vs-NegEV bracket sits over
      the pair it tests (same dark/light shading convention as the NTA
      TEV-vs-EEV figure) - added 2026-09-21, per Revital's request, to show
      whether the EEV-over-NegEV effect holds within each donor group
      separately, not just pooled.

Same uncorrected pairwise Mann-Whitney, same labeling conventions as the
July and September figures. Uses only 2 groups/1 bracket in panels A and B
(the "IRI Control" bar was REMOVED 2026-09-21, per Revital's explicit call)
- see the normalization note below for why a shared/pooled control bar was
potentially confusing here specifically.

***IMPORTANT - fold-rescue normalization across experiments***: every value
plotted anywhere in this figure (EEV, NegEV, Healthy, Diabetic - in all
three panels) is EACH ALREADY EXPRESSED AS A FOLD RELATIVE TO ITS OWN
EXPERIMENT'S MEAN (July wells / July's own IRI control mean = 0.2371M;
September wells / September's own IRI control mean = 7.1957M) BEFORE being
pooled together. July and September have very different absolute RLU
scales (a protocol/instrument difference between runs, not biological), so
combining raw values directly would be meaningless - normalizing each
experiment to its own control first puts both onto the same unitless scale.
This must be stated explicitly in the caption/text (Revital's call,
2026-09-21) - it's the single most important methodological fact about this
figure, not an aside.

Pooling is literal concatenation of each experiment's own fold-rescue
values - verified 2026-09-21 that this reproduces note 11.3's finalized
n's and stats exactly (EEV n=115=44+71, NegEV n=89=27+62, Healthy EEV
n=63=20+43, Diabetic EEV n=52=24+28).

Overlapping donors (per note 11.3's "Overlapping Donors" decisions):
  - H10 (Healthy) and P5 (Diabetic): BOTH experiments' wells included,
    tracked separately by marker shape (matches note 11.3's own convention:
    circle = September, triangle = July), same donor color in both.
  - 3T (Diabetic): JULY ONLY - September's 3T data is the same
    pipetting-error-affected run already excluded from the September
    figure entirely; not re-included here via a different route.
  - P4: excluded from both experiments (never had usable July data; September
    data is pipetting-error-affected).

Donor palette: identical DONOR_COLORS dict to the September script (already
built to cover every donor across both experiments) - H10 and P5 keep one
color regardless of which experiment a given well came from; marker shape
alone distinguishes the experiment.

Data source: raw plate-reader files already copied into Data/July/ and
Data/September/ for the individual-experiment figures.

Output: Combined_EEV_vs_NegEV_and_HvsD.png/.svg (one level up, in Figures/)
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

# Diabetic codes use the unified D-naming convention (1T->D1, 2T->D2,
# 3T->D3, P5->D5, P6->D6, P7->D7) - Healthy stays H.
DONOR_COLORS = {
    "H8": "#123c69", "H11": "#1c5cab", "H9": "#2a78d6", "H2": "#4a9de8",
    "H4": "#6fb8f0", "H10": "#8ec4f0", "H3": "#a8d5ff", "H1": "#c3e4ff",
    "D2": "#7a1f2b", "D1": "#a83246", "D3": "#e34948", "D6": "#eb6b3d",
    "D7": "#f2914f", "D5": "#e8a33d",
}
DONOR_ORDER = ["H8", "H11", "H10", "H9", "H2", "H3", "H4", "H1",
               "D2", "D1", "D3", "D6", "D7", "D5"]
HEALTHY_DONORS = ["H8", "H11", "H10", "H9", "H2", "H3", "H4", "H1"]
DIABETIC_DONORS = ["D2", "D1", "D3", "D6", "D7", "D5"]

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
JULY_DIR = os.path.join(FIG_DIR, "..", "Data", "July")
SEPT_DIR = os.path.join(FIG_DIR, "..", "Data", "September")

RNG_SEED = 42
ROWS = list("ABCDEFGH")


def tukey_filter(vals):
    a = np.asarray(vals, dtype=float)
    if len(a) < 4:
        return a
    q1, q3 = np.percentile(a, 25), np.percentile(a, 75)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    return a[(a >= lo) & (a <= hi)]


# ═══════════════════════════ JULY ═══════════════════════════════════════════
XLSX_JULY_IRI = os.path.join(JULY_DIR, "Viability_post_IRI 12.07.26 read 1.xlsx")
JULY_PRE_EXCLUDED = {("A", 1), ("B", 1), ("C", 1), ("A", 2), ("B", 2), ("E", 5)}

wb_j = openpyxl.load_workbook(XLSX_JULY_IRI)
ws_j = wb_j.active
july_raw = {}
for ri, r in enumerate(ROWS):
    for col in range(1, 13):
        july_raw[(r, col)] = ws_j.cell(row=11 + ri, column=col + 1).value


def july_wells_M(col):
    vals = []
    for r in ROWS:
        if (r, col) in JULY_PRE_EXCLUDED:
            continue
        v = july_raw.get((r, col))
        if v is not None:
            vals.append(v / 1e6)
    return np.array(vals)


july_iri_ctrl_vals = july_wells_M(11)
JULY_IRI_CTRL_MEAN = july_iri_ctrl_vals.mean()
july_iri_ctrl_fold = tukey_filter(july_iri_ctrl_vals) / JULY_IRI_CTRL_MEAN

JULY_DONOR_COLS = {
    "H8": {"EEV": 5}, "H9": {"EEV": 2, "NegEV": 1}, "H10": {"EEV": 4, "NegEV": 3},
    "D2": {"EEV": 10}, "D3": {"EEV": 7, "NegEV": 6}, "D5": {"EEV": 9, "NegEV": 8},
}
july_donor_data = {}
for donor, cols in JULY_DONOR_COLS.items():
    july_donor_data[donor] = {}
    for typ, col in cols.items():
        july_donor_data[donor][typ] = tukey_filter(july_wells_M(col)) / JULY_IRI_CTRL_MEAN

# ═══════════════════════════ SEPTEMBER ══════════════════════════════════════
XLSX_SEPT_HEALTHY = os.path.join(SEPT_DIR, "06.09.26 Viabilty Healthy plate.xlsx")
XLSX_SEPT_DIABETIC = os.path.join(SEPT_DIR, "06.09.26 Viabilty Diabetic plate.xlsx")


def read_plate(path):
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb["Luminescence 1_01"]
    raw = {}
    for ri, r in enumerate(ROWS):
        for col in range(1, 13):
            raw[(r, col)] = ws.cell(row=11 + ri, column=col + 1).value
    return raw


def sept_wells_M(raw, col, drop_rows=()):
    vals = []
    for r in ROWS:
        if r in drop_rows:
            continue
        v = raw.get((r, col))
        if v is not None:
            vals.append(v / 1e6)
    return np.array(vals)


sept_healthy_raw = read_plate(XLSX_SEPT_HEALTHY)
sept_diabetic_raw = read_plate(XLSX_SEPT_DIABETIC)

sept_iri_ctrl_raw = sept_wells_M(sept_diabetic_raw, 12, drop_rows={"A"})
SEPT_IRI_CTRL_MEAN = tukey_filter(sept_iri_ctrl_raw).mean()
sept_iri_ctrl_fold = tukey_filter(sept_iri_ctrl_raw) / SEPT_IRI_CTRL_MEAN

SEPT_DONOR_COLS = {
    "H11": {"plate": sept_healthy_raw, "EEV": 1, "NegEV": 2},
    "H10": {"plate": sept_healthy_raw, "EEV": (3, {"H"}), "NegEV": 4},
    "H2": {"plate": sept_healthy_raw, "EEV": 5, "NegEV": 6},
    "H3": {"plate": sept_healthy_raw, "EEV": 7, "NegEV": 8},
    "H4": {"plate": sept_healthy_raw, "EEV": 9, "NegEV": 10},
    "H1": {"plate": sept_healthy_raw, "EEV": (11, {"H"})},
    "D5": {"plate": sept_diabetic_raw, "EEV": 5, "NegEV": 6},
    "D6": {"plate": sept_diabetic_raw, "EEV": 7, "NegEV": 8},
    "D7": {"plate": sept_diabetic_raw, "EEV": 9, "NegEV": 10},
    "D1": {"plate": sept_diabetic_raw, "EEV": 11},
}
sept_donor_data = {}
for donor, spec in SEPT_DONOR_COLS.items():
    sept_donor_data[donor] = {}
    plate = spec["plate"]
    for typ in ("EEV", "NegEV"):
        if typ not in spec:
            continue
        col_spec = spec[typ]
        col, drop_rows = col_spec if isinstance(col_spec, tuple) else (col_spec, set())
        sept_donor_data[donor][typ] = tukey_filter(sept_wells_M(plate, col, drop_rows=drop_rows)) / SEPT_IRI_CTRL_MEAN

# ═══════════════════════════ COMBINE ════════════════════════════════════════
# donor -> typ -> list of (values_array, experiment_tag) - D3 (3T) is July-only,
# per note 11.3 (September's D3/3T run is the pipetting-error-affected one,
# already excluded from the September figure/donor set).
combined = {}
for donor, data in july_donor_data.items():
    combined.setdefault(donor, {})
    for typ, vals in data.items():
        combined[donor].setdefault(typ, []).append((vals, "July"))
for donor, data in sept_donor_data.items():
    if donor == "D3":
        continue  # not applicable here (D3/3T has no September entry anyway)
    combined.setdefault(donor, {})
    for typ, vals in data.items():
        combined[donor].setdefault(typ, []).append((vals, "September"))


def pooled_combined(typ, donors):
    out = []
    for d in donors:
        for vals, _exp in combined.get(d, {}).get(typ, []):
            out.append(vals)
    return np.concatenate(out) if out else np.array([])


iri_ctrl_fold = np.concatenate([july_iri_ctrl_fold, sept_iri_ctrl_fold])
eev_all = pooled_combined("EEV", DONOR_ORDER)
negev_all = pooled_combined("NegEV", DONOR_ORDER)
eev_healthy = pooled_combined("EEV", HEALTHY_DONORS)
eev_diabetic = pooled_combined("EEV", DIABETIC_DONORS)

print(f"n check: EEV={len(eev_all)} (expect 115), NegEV={len(negev_all)} (expect 89), "
      f"Healthy EEV={len(eev_healthy)} (expect 63), Diabetic EEV={len(eev_diabetic)} (expect 52)")


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


def plot_panel(ax, groups, x_positions, bar_colors, bracket_pairs, letter, title):
    """groups: list of dicts {label, pool, by_donor}. x_positions: x for each
    group (allows a visual gap, e.g. [0,1,2.6,3.6]). bracket_pairs: list of
    (i, j) index pairs to test/annotate; non-overlapping x-ranges get drawn
    at the same (locally-scaled) height, since they can't visually collide."""
    ax.set_facecolor(BG)
    rng = np.random.default_rng(RNG_SEED)
    x = np.array(x_positions)
    means = [g["pool"].mean() for g in groups]
    sds = [g["pool"].std(ddof=1) for g in groups]
    ns = [len(g["pool"]) for g in groups]

    ax.bar(x, means, yerr=sds, color=bar_colors, edgecolor=bar_colors,
           linewidth=0.8, error_kw=dict(elinewidth=1.4, capsize=5, capthick=1.4,
           ecolor=LABEL_COLOR), width=0.55, zorder=3, alpha=0.65)

    for xi, g in zip(x, groups):
        for donor in DONOR_ORDER:
            if donor not in g["by_donor"]:
                continue
            for vals, exp in g["by_donor"][donor]:
                if len(vals) == 0:
                    continue
                marker = "^" if exp == "July" else "o"
                jit = rng.uniform(-0.16, 0.16, len(vals))
                ax.scatter(np.full(len(vals), xi) + jit, vals, color=DONOR_COLORS[donor],
                           edgecolors="white", linewidths=0.6, s=48, zorder=5, marker=marker)

    ax.axhline(1.0, color=LABEL_COLOR, linestyle=(0, (4, 3)), linewidth=1.1, zorder=2)

    max_bracket_top = 0
    for i, j in bracket_pairs:
        local_top = max(max(groups[i]["pool"]), max(groups[j]["pool"]),
                         means[i] + sds[i], means[j] + sds[j])
        y = local_top * 1.15
        p = stats.mannwhitneyu(groups[i]["pool"], groups[j]["pool"], alternative="two-sided").pvalue
        sig_bracket(ax, x[i], x[j], y, sig_label(p))
        print(f"  {groups[i]['label']} vs {groups[j]['label']}: p={p:.5f}")
        max_bracket_top = max(max_bracket_top, y)

    ax.set_ylim(0, max_bracket_top * 1.15)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g['label']}\n(n={n})" for g, n in zip(groups, ns)],
                        fontsize=14, fontfamily=FONT, color=TITLE_COLOR)
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)
    style_axes(ax)
    ax.set_title(title, fontsize=16, fontweight="bold", fontfamily=FONT,
                 color=TITLE_COLOR, pad=14, loc="center")
    add_panel_letter(ax, letter)


fig, axes = plt.subplots(1, 3, figsize=(19, 7.8))
fig.patch.set_facecolor(BG)

by_donor = lambda typ, donors: {d: combined[d][typ] for d in donors if typ in combined.get(d, {})}

print("Panel A (EEV / NegEV):")
panel_a_groups = [
    dict(label="EEV", pool=eev_all, by_donor=by_donor("EEV", DONOR_ORDER)),
    dict(label="NegEV", pool=negev_all, by_donor=by_donor("NegEV", DONOR_ORDER)),
]
plot_panel(axes[0], panel_a_groups, [0, 1], [C_EEV, C_NEG], [(0, 1)],
           "A", "EEV vs. NegEV Fold Rescue - Combined")

print("Panel B (Healthy EEV / Diabetic EEV):")
panel_b_groups = [
    dict(label="Healthy EEV", pool=eev_healthy, by_donor=by_donor("EEV", HEALTHY_DONORS)),
    dict(label="Diabetic EEV", pool=eev_diabetic, by_donor=by_donor("EEV", DIABETIC_DONORS)),
]
plot_panel(axes[1], panel_b_groups, [0, 1], [C_HEALTHY, C_DIABETIC], [(0, 1)],
           "B", "EEV Fold Rescue - Healthy vs. Diabetic - Combined")

print("Panel C (Healthy EEV / Healthy NegEV / Diabetic EEV / Diabetic NegEV):")
negev_healthy = pooled_combined("NegEV", HEALTHY_DONORS)
negev_diabetic = pooled_combined("NegEV", DIABETIC_DONORS)
panel_c_groups = [
    dict(label="Healthy\nEEV", pool=eev_healthy, by_donor=by_donor("EEV", HEALTHY_DONORS)),
    dict(label="Healthy\nNegEV", pool=negev_healthy, by_donor=by_donor("NegEV", HEALTHY_DONORS)),
    dict(label="Diabetic\nEEV", pool=eev_diabetic, by_donor=by_donor("EEV", DIABETIC_DONORS)),
    dict(label="Diabetic\nNegEV", pool=negev_diabetic, by_donor=by_donor("NegEV", DIABETIC_DONORS)),
]
# Grouped by disease group (Healthy pair, gap, Diabetic pair) so each
# within-group EEV-vs-NegEV bracket sits directly over the pair it tests -
# same TEV/EEV-style dark/light shading used for the NTA diameter figure.
plot_panel(axes[2], panel_c_groups, [0, 1, 2.6, 3.6],
           [C_HEALTHY, "#6aadff", C_DIABETIC, "#eb6834"], [(0, 1), (2, 3)],
           "C", "EEV vs. NegEV - Healthy vs. Diabetic - Combined")

for ax in axes:
    ax.set_ylabel("Fold Rescue (relative to IRI control)", fontsize=15.5,
                  fontfamily=FONT, color=LABEL_COLOR, labelpad=10)

donor_legend = [plt.scatter([], [], color=DONOR_COLORS[d], edgecolors="white",
                             linewidths=0.6, s=60, label=d) for d in DONOR_ORDER]
exp_legend = [
    plt.scatter([], [], color="#888888", edgecolors="white", linewidths=0.6, s=60, marker="^", label="July"),
    plt.scatter([], [], color="#888888", edgecolors="white", linewidths=0.6, s=60, marker="o", label="September"),
]

leg1 = fig.legend(handles=donor_legend, loc="lower center", ncol=14, fontsize=10,
                   frameon=False, bbox_to_anchor=(0.5, -0.10), handletextpad=0.4,
                   columnspacing=1.0, title="Donor", title_fontsize=11.5,
                   prop={"family": FONT})
leg1.get_title().set_fontfamily(FONT)
leg2 = fig.legend(handles=exp_legend, loc="lower center", ncol=2, fontsize=11,
                   frameon=False, bbox_to_anchor=(0.5, -0.20), handletextpad=0.5,
                   columnspacing=2.0, title="Experiment", title_fontsize=11.5,
                   prop={"family": FONT})
leg2.get_title().set_fontfamily(FONT)
fig.add_artist(leg1)

fig.suptitle("Combined (July + September) IRI Experiments - Functional Rescue",
             fontsize=22, fontweight="bold", color=TITLE_COLOR, y=1.02, fontfamily=FONT)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"Combined_EEV_vs_NegEV_and_HvsD.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print("\nSaved: Combined_EEV_vs_NegEV_and_HvsD.png/.svg")
