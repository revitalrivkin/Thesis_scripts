"""
4.1.3 IRI Functional Assay Results - July Experiment, Raw Viability Overview
================================================================================
Straightforward raw-data overview for the July 2026 IRI experiment: baseline
(no IRI), IRI control (no EVs), and every donor's EEV/Negative EV wells,
before any fold-rescue normalization or statistical comparison. Intended as
an orienting figure ahead of the Healthy-vs-Diabetic and EEV-vs-NegEV
comparison figures that follow.

Deliberately simplified relative to the exploratory draft version
(`Plasma EVs\IRI\July IRI\Viability\plot_July_raw_all_wells.py`), which used
a custom Tukey-fence box+whisker visualization with outlier callouts - not
this thesis's house style. This version uses the standard house chart type
(bar = group mean, error bar = SD, individual wells as jittered dots) to
stay visually consistent with Figures 1-3.

Both the hard, non-statistical exclusions (pipetting errors, empty/seeding-
error wells) and statistical outliers (Tukey's fence, 1.5x IQR per donor
group) are removed here, matching the same exclusion criteria used for the
downstream fold-rescue comparisons (note 11.1) - added 2026-09-20, per
Revital's call, so group means shown here match what the later comparison
figures actually use, rather than showing a "rawer" version that would
disagree with its own downstream numbers. Baseline and IRI control are
NOT Tukey-filtered, matching note 11.1 (both reference means use all
wells, n=16/n=8).

Data source: ../../Data/July/Viability_post_IRI 12.07.26 read 1.xlsx,
             ../../Data/July/Baseline viability 11.07.26 read 1.xlsx
(copied from the draft workspace, Plasma EVs\IRI\July IRI\Viability\data\)

Output: July_Raw_Overview.png/.svg (one level up, in Figures/)
"""

import os
import numpy as np
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT = "Arial"
BG = "#fcfcfb"
GRID = "#e1e0d9"
SPINE = "#c3c2b7"
TITLE_COLOR = "#0b0b0b"
LABEL_COLOR = "#52514e"

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
DATA_DIR = os.path.join(FIG_DIR, "..", "Data", "July")
XLSX_IRI = os.path.join(DATA_DIR, "Viability_post_IRI 12.07.26 read 1.xlsx")
XLSX_BASE = os.path.join(DATA_DIR, "Baseline viability 11.07.26 read 1.xlsx")

RNG_SEED = 42
ROWS = list("ABCDEFGH")

# ── Hard exclusions only (pipetting errors / empty wells) -- not statistical
# outliers -- per note 11.1 "Exclusions (finalized)" ────────────────────────
PRE_EXCLUDED = {
    ("A", 1), ("B", 1), ("C", 1),   # H9 NegEV -- pipetting error
    ("A", 2), ("B", 2),             # H9 EEV -- empty well, seeding error
    ("E", 5),                       # H8 EEV -- empty well, seeding error
}

# ── Read raw plate data ──────────────────────────────────────────────────────
wb = openpyxl.load_workbook(XLSX_IRI)
ws = wb.active
raw = {}
for ri, r in enumerate(ROWS):
    for col in range(1, 13):
        raw[(r, col)] = ws.cell(row=11 + ri, column=col + 1).value

wb2 = openpyxl.load_workbook(XLSX_BASE)
ws2 = wb2.active
# Columns 8 and 9 are TWO PLATE-READER REPS OF THE SAME PHYSICAL WELL per row
# (rows A-H) - NOT 16 independent wells. Corrected 2026-09-20: previously
# pooled both columns as 16 separate values; per the established reference
# script (Plasma EVs\IRI\July IRI\Viability\plot_Viability_bars.py, comment:
# "BASELINE plate: cols 8 & 9 = 2 plate-reader reps per well (rows A-H)"),
# the two columns must be averaged per row first, giving n=8 (one value per
# physical well), matching every other group's well count in this figure.
baseline_vals = []
for ri, r in enumerate(ROWS):
    reps = [ws2.cell(row=11 + ri, column=col + 1).value for col in [8, 9]]
    reps = [v for v in reps if v is not None]
    if reps:
        baseline_vals.append(np.mean(reps) / 1e6)
baseline_vals = np.array(baseline_vals)


def tukey_filter(vals):
    """Drop statistical outliers via Tukey's fence (1.5x IQR), matching the
    exclusion method used for the downstream fold-rescue comparisons (note
    11.1). Applied only to donor EEV/NegEV groups, not baseline/IRI ctrl,
    consistent with note 11.1 (both reference means use all wells, n=16/n=8)."""
    a = np.asarray(vals, dtype=float)
    if len(a) < 4:
        return a, 0
    q1, q3 = np.percentile(a, 25), np.percentile(a, 75)
    iqr = q3 - q1
    lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    kept = a[(a >= lo) & (a <= hi)]
    return kept, len(a) - len(kept)


def get_wells(col):
    vals = []
    for r in ROWS:
        if (r, col) in PRE_EXCLUDED:
            continue
        v = raw.get((r, col))
        if v is not None:
            vals.append(v / 1e6)
    return np.array(vals)


BASELINE_MEAN = baseline_vals.mean()
iri_ctrl_vals = get_wells(11)
IRI_CTRL_MEAN = iri_ctrl_vals.mean()

# ── Group definitions, in display order ──────────────────────────────────────
# (label, column, color, dot_color)
C_BASE, C_BASE_DOT = "#3ba272", "#c9ead9"
C_CTRL, C_CTRL_DOT = "#7a8694", "#dadde1"
C_H_EEV, C_H_EEV_DOT = "#2a78d6", "#cde2fb"
C_H_NEG, C_H_NEG_DOT = "#6aadff", "#d9ebff"
C_D_EEV, C_D_EEV_DOT = "#e34948", "#f5b8b8"
C_D_NEG, C_D_NEG_DOT = "#eb6834", "#fbd2bb"

# Donor EEV/NegEV groups get Tukey-filtered; baseline/IRI ctrl do not
# (matches note 11.1: both reference means use all wells, n=16/n=8)

# Diabetic donor codes use the unified D-naming convention (D1-D7), matching
# the clinical/NTA data (1T->D1, 2T->D2, 3T->D3, P4->D4, P5->D5, P6->D6,
# P7->D7) - Healthy stays H. Old lab codes (2T/3T/P5) are used only as raw
# column lookups below; display labels use the new names.
donor_raw = {
    "H8\nEEV": get_wells(5), "H9\nEEV": get_wells(2), "H9\nNegEV": get_wells(1),
    "H10\nEEV": get_wells(4), "H10\nNegEV": get_wells(3),
    "D2\nEEV": get_wells(10), "D3\nEEV": get_wells(7), "D3\nNegEV": get_wells(6),
    "D5\nEEV": get_wells(9), "D5\nNegEV": get_wells(8),
}
donor_filtered = {}
excluded_log = []
for label, vals in donor_raw.items():
    kept, n_dropped = tukey_filter(vals)
    donor_filtered[label] = kept
    if n_dropped:
        excluded_log.append(f"  {label.replace(chr(10), ' ')}: {n_dropped} well(s) dropped (of {len(vals)})")

TOTAL_EXCLUDED = sum(len(v) - len(donor_filtered[k]) for k, v in donor_raw.items())

groups = [
    ("Baseline\n(no IRI)", baseline_vals, C_BASE, C_BASE_DOT),
    ("IRI ctrl\n(no EVs)", iri_ctrl_vals, C_CTRL, C_CTRL_DOT),
    ("H8\nEEV", donor_filtered["H8\nEEV"], C_H_EEV, C_H_EEV_DOT),
    ("H9\nEEV", donor_filtered["H9\nEEV"], C_H_EEV, C_H_EEV_DOT),
    ("H9\nNegEV", donor_filtered["H9\nNegEV"], C_H_NEG, C_H_NEG_DOT),
    ("H10\nEEV", donor_filtered["H10\nEEV"], C_H_EEV, C_H_EEV_DOT),
    ("H10\nNegEV", donor_filtered["H10\nNegEV"], C_H_NEG, C_H_NEG_DOT),
    ("D2\nEEV", donor_filtered["D2\nEEV"], C_D_EEV, C_D_EEV_DOT),
    ("D3\nEEV", donor_filtered["D3\nEEV"], C_D_EEV, C_D_EEV_DOT),
    ("D3\nNegEV", donor_filtered["D3\nNegEV"], C_D_NEG, C_D_NEG_DOT),
    ("D5\nEEV", donor_filtered["D5\nEEV"], C_D_EEV, C_D_EEV_DOT),
    ("D5\nNegEV", donor_filtered["D5\nNegEV"], C_D_NEG, C_D_NEG_DOT),
]

rng = np.random.default_rng(RNG_SEED)
fig, ax = plt.subplots(figsize=(15, 6.5))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

x = np.arange(len(groups))
means = [g[1].mean() for g in groups]
sds = [g[1].std(ddof=1) for g in groups]
bar_colors = [g[2] for g in groups]
dot_colors = [g[3] for g in groups]
ns = [len(g[1]) for g in groups]

ax.bar(x, means, yerr=sds,
       color=bar_colors, edgecolor=bar_colors, linewidth=0.8,
       error_kw=dict(elinewidth=1.4, capsize=4, capthick=1.4, ecolor=LABEL_COLOR),
       width=0.6, zorder=3, alpha=0.80)

for xi, (label, vals, _, dcolor) in zip(x, groups):
    jit = rng.uniform(-0.14, 0.14, len(vals))
    ax.scatter(np.full(len(vals), xi) + jit, vals, color=dcolor,
               edgecolors="white", linewidths=0.5, s=36, zorder=5)

# Reference lines
ax.axhline(BASELINE_MEAN, color=C_BASE, linestyle=(0, (4, 3)), linewidth=1.1, zorder=2)
ax.axhline(IRI_CTRL_MEAN, color=C_CTRL, linestyle=(0, (4, 3)), linewidth=1.1, zorder=2)

ax.set_axisbelow(True)
ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax.set_xticks(x)
ax.set_xticklabels([f"{g[0]}\n(n={n})" for g, n in zip(groups, ns)],
                    fontsize=14, fontfamily=FONT, color=TITLE_COLOR)
ax.set_xlim(-0.6, len(groups) - 0.4)
ax.set_ylabel("Raw Luminescence (millions)", fontsize=15, fontfamily=FONT,
              color=LABEL_COLOR, labelpad=10)
ax.tick_params(axis="y", labelsize=14, labelcolor=LABEL_COLOR)
ax.tick_params(axis="x", bottom=False)
for lbl in ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
for spine in ["top", "right", "bottom"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color(SPINE)
ax.spines["left"].set_linewidth(0.8)

# Group section labels (Healthy / Diabetic) below the donor x-tick labels
trans = ax.get_xaxis_transform()
ax.text(4, -0.27, "Healthy donors", ha="center", va="top", fontsize=15,
        fontweight="bold", color=C_H_EEV, fontfamily=FONT, transform=trans, clip_on=False)
ax.text(9, -0.27, "Diabetic donors", ha="center", va="top", fontsize=15,
        fontweight="bold", color=C_D_EEV, fontfamily=FONT, transform=trans, clip_on=False)

ax.set_title("July 2026 IRI Experiment - Raw Viability Signal, All Groups",
             fontsize=18, fontweight="bold", fontfamily=FONT,
             color=TITLE_COLOR, pad=18, loc="center")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"July_Raw_Overview.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)

print(f"Baseline mean: {BASELINE_MEAN:.3f} M (n={len(baseline_vals)})")
print(f"IRI ctrl mean: {IRI_CTRL_MEAN:.3f} M (n={len(iri_ctrl_vals)})")
print(f"\nTukey fence exclusions (donor groups only): {TOTAL_EXCLUDED} well(s) total")
for line in excluded_log:
    print(line)
print("\nSaved: July_Raw_Overview.png/.svg")
