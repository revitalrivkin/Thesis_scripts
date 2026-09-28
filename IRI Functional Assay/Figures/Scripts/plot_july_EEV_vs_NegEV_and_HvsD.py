"""
4.1.3 IRI Functional Assay Results - July: EEV vs. NegEV, Healthy vs. Diabetic EEV
================================================================================
Two-panel figure, July 2026 IRI experiment only (September and Combined get
their own analogous figures later). Each panel now has THREE groups (IRI
control included) and three pairwise brackets, matching the structure of
the Faculty Retreat slides `Presentations\Faculty retreat\Viability results\
Viability_EEV_vs_NegEV_slide.png` and `Viability_H_vs_D_slide.png`
(Revital's reference, 2026-09-20):
  (A) IRI control / EEV / NegEV fold rescue, all donors pooled by type.
  (B) IRI control / Healthy EEV / Diabetic EEV fold rescue.

Y-axis is normalized fold rescue over the IRI control (RLU_well /
IRI_ctrl_mean), per Sections 3.7/3.8 - not raw RLU. IRI control = 1x by
definition, marked as a reference line AND shown as its own bar/group
(Control wells' own well-level spread around 1x), so it can be directly
statistically compared to EEV and NegEV.

Colors: Control = pale grey, NegEV = a darker grey (same family, since both
are "non-endothelial or no-treatment" comparators) - Revital's explicit
call, 2026-09-20 - EEV stays green, Healthy/Diabetic keep the established
blue/red pair.

Dots are colored per-donor, nested within group color family (blue shades
for Healthy, warm red/orange for Diabetic) - Control dots use a plain
neutral grey, since control wells aren't derived from any of the 6 donors.

Statistics: Mann-Whitney U (well-level), all three pairwise comparisons per
panel, UNCORRECTED for multiple comparisons - deliberately NOT using the
Bonferroni correction the retreat's H-vs-D slide applied, per Revital's
explicit call (2026-09-20): "do not complicate with corrections, leave as
it was" - matches every other statistical comparison in this thesis so far
(none of which apply a multiple-comparison correction). Exact p-value shown
when significant (p<0.05), "ns" otherwise - not asterisks.

Same Tukey's-fence (1.5x IQR per donor group) + hard exclusions as the raw
overview figure - IRI control itself is not Tukey-filtered, matching note
11.1.

Data source: ../../Data/July/Viability_post_IRI 12.07.26 read 1.xlsx

Output: July_EEV_vs_NegEV_and_HvsD.png/.svg (one level up, in Figures/)
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

C_CONTROL = "#c7cbd1"   # pale grey
C_CONTROL_DOT = "#9aa0a8"
C_NEG = "#7d838c"       # darker grey (same family as control)
C_EEV = "#3ba272"
C_HEALTHY, C_DIABETIC = "#2a78d6", "#e34948"

# Donor colors: per-donor shade, nested within group color family.
# Diabetic codes use the unified D-naming convention (D1-D7), matching the
# clinical/NTA data (1T->D1, 2T->D2, 3T->D3, P4->D4, P5->D5, P6->D6, P7->D7)
# - Healthy stays H.
DONOR_COLORS = {
    "H8": "#123c69", "H9": "#2a78d6", "H10": "#8ec4f0",
    "D2": "#7a1f2b", "D3": "#e34948", "D5": "#e8a33d",
}
DONOR_ORDER = ["H8", "H9", "H10", "D2", "D3", "D5"]
HEALTHY_DONORS = ["H8", "H9", "H10"]
DIABETIC_DONORS = ["D2", "D3", "D5"]

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
DATA_DIR = os.path.join(FIG_DIR, "..", "Data", "July")
XLSX_IRI = os.path.join(DATA_DIR, "Viability_post_IRI 12.07.26 read 1.xlsx")

RNG_SEED = 42
ROWS = list("ABCDEFGH")
PRE_EXCLUDED = {("A", 1), ("B", 1), ("C", 1), ("A", 2), ("B", 2), ("E", 5)}

wb = openpyxl.load_workbook(XLSX_IRI)
ws = wb.active
raw = {}
for ri, r in enumerate(ROWS):
    for col in range(1, 13):
        raw[(r, col)] = ws.cell(row=11 + ri, column=col + 1).value


def get_wells_M(col):
    vals = []
    for r in ROWS:
        if (r, col) in PRE_EXCLUDED:
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


iri_ctrl_vals = get_wells_M(11)
IRI_CTRL_MEAN = iri_ctrl_vals.mean()
iri_ctrl_fold = iri_ctrl_vals / IRI_CTRL_MEAN

DONOR_COLS = {
    "H8": {"EEV": 5}, "H9": {"EEV": 2, "NegEV": 1}, "H10": {"EEV": 4, "NegEV": 3},
    "D2": {"EEV": 10}, "D3": {"EEV": 7, "NegEV": 6}, "D5": {"EEV": 9, "NegEV": 8},
}

donor_data = {}
for donor, cols in DONOR_COLS.items():
    donor_data[donor] = {}
    for typ, col in cols.items():
        kept = tukey_filter(get_wells_M(col))
        donor_data[donor][typ] = kept / IRI_CTRL_MEAN


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
    """groups: list of 3 dicts {label, pool, by_donor (or None for Control)}"""
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

    # Three pairwise comparisons: (0,1), (1,2), (0,2) - uncorrected Mann-Whitney
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
plot_panel(axes[0], panel_a_groups, [C_CONTROL, C_EEV, C_NEG], "A", "EEV vs. NegEV Fold Rescue - July")

print("Panel B (IRI ctrl / Healthy EEV / Diabetic EEV):")
panel_b_groups = [
    dict(label="IRI Control", pool=iri_ctrl_fold, by_donor=None),
    dict(label="Healthy EEV", pool=eev_healthy, by_donor={d: donor_data[d]["EEV"] for d in HEALTHY_DONORS}),
    dict(label="Diabetic EEV", pool=eev_diabetic, by_donor={d: donor_data[d]["EEV"] for d in DIABETIC_DONORS}),
]
plot_panel(axes[1], panel_b_groups, [C_CONTROL, C_HEALTHY, C_DIABETIC], "B", "EEV Fold Rescue - Healthy vs. Diabetic - July")

for ax in axes:
    ax.set_ylabel("Fold Rescue (relative to IRI control)", fontsize=12.5,
                  fontfamily=FONT, color=LABEL_COLOR, labelpad=10)

# Shared donor legend, ordered Healthy (blue shades) then Diabetic (warm shades)
legend_handles = [plt.scatter([], [], color=DONOR_COLORS[d], edgecolors="white",
                               linewidths=0.6, s=60, label=d) for d in DONOR_ORDER]
_leg = fig.legend(handles=legend_handles, loc="lower center", ncol=6, fontsize=11.5,
                  frameon=False, bbox_to_anchor=(0.5, -0.06), handletextpad=0.5,
                  columnspacing=1.4, title="Donor", title_fontsize=12,
                  prop={"family": FONT})
_leg.get_title().set_fontfamily(FONT)

fig.suptitle("July 2026 IRI Experiment - Functional Rescue", fontsize=17,
             fontweight="bold", color=TITLE_COLOR, y=1.02, fontfamily=FONT)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"July_EEV_vs_NegEV_and_HvsD.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print("\nSaved: July_EEV_vs_NegEV_and_HvsD.png/.svg")
