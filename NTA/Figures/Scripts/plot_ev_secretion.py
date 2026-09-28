"""
4.1.2 Plasma-EV Characterization - EV Secretion and %EEV Summary Figure
================================================================================
Two-panel figure: (A) particle concentration ("EV secretion", particles per
mL plasma) for TEV and EEV, Healthy vs. Diabetic; (B) %EEVs from total EVs,
Healthy vs. Diabetic. Styled to match Figure 1 (Donor_Cohort_Summary): Arial,
off-white background, bar + mean+-SD error bar + jittered per-donor scatter
dots, bold titles, plain hyphens (no em dashes), bold panel letters (A/B).

Panel A is grouped by fraction (TEV pair, then EEV pair), not by disease
group, so each Healthy-vs-Diabetic "ns" significance bracket sits directly
over the pair it actually tests, rather than spanning across the other
fraction's bar (added 2026-09-16, Revital's request to show comparisons
directly on the figure rather than only in the caption/text).

Panel A uses a logarithmic y-axis: Diabetic TEV concentration includes one
donor (D2, old code 2T) with a concentration ~30-300x higher than the rest of
the cohort, independently confirmed against the original ZetaView reports as
real (see note 16, "2T TEV Outlier"), not a data error. On a linear axis this
outlier collapses every other bar to near-zero and forces a mean-SD error bar
below zero (not physically meaningful for a particle count) - see the
superseded Draft_EVSecretion_v1_linear.png for the direct comparison.
Mean +/- SD is kept for all panels (including TEV) for consistency with the
rest of the thesis's NTA figures, per Revital's explicit choice 2026-09-16,
despite TEV's skew (Diabetic TEV Shapiro-Wilk p=0.0001) - only the upper half
of the SD whisker is drawn on the log axis, since a two-sided linear SD
whisker is not meaningful on a log scale.

Data source: ../../Plasma EV sample records 12.8 (D-labeled).csv

Output: EV_Secretion_and_Percent_EEV.png/.svg (one level up, in Figures/)
"""

import os
import numpy as np
import pandas as pd
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

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
DATA = os.path.join(FIG_DIR, "..", "Plasma EV sample records 12.8 (D-labeled).csv")

RNG_SEED = 42

df = pd.read_csv(DATA)
df.columns = [c.strip() for c in df.columns]
df["Group_label"] = df["Group"].map({"Control": "Healthy", "Diabetic": "Diabetic"})

groups = ["Healthy", "Diabetic"]


def get_vals(col, group_label):
    return df.loc[df["Group_label"] == group_label, col].dropna().values


tev = [get_vals("Normalized TEV [total Evs/mL plasma]", g) for g in groups]
eev = [get_vals("EEVs per mL plasma", g) for g in groups]
pct = [get_vals("%EEVs from total EVs", g) for g in groups]

# ── Stats (for reference / thesis text; not rendered on the figure) ────────
log_lines = []


def report(name, a, b):
    log_lines.append(f"\n{name}")
    for label, vals in zip(groups, [a, b]):
        sh = stats.shapiro(vals)
        log_lines.append(f"  {label}: n={len(vals)} mean={vals.mean():.3e} sd={vals.std(ddof=1):.3e} "
                          f"median={np.median(vals):.3e} Shapiro p={sh.pvalue:.4f}")
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    log_lines.append(f"  Mann-Whitney U={u.statistic:.1f} p={u.pvalue:.4f}")


report("TEV concentration (particles/mL plasma), Healthy vs. Diabetic", tev[0], tev[1])
report("EEV concentration (particles/mL plasma), Healthy vs. Diabetic", eev[0], eev[1])
report("%EEVs from total EVs, Healthy vs. Diabetic", pct[0], pct[1])

with open(os.path.join(HERE, "ev_secretion_stats.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))


def add_panel_letter(ax, letter):
    ax.text(-0.10, 1.10, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", color=TITLE_COLOR, fontfamily=FONT,
            va="bottom", ha="left")


def sig_bracket_log(ax, x1, x2, y_top, label, pad_mult=1.35, tick_mult=1.10):
    """Draw a significance bracket for a log-scale axis (multiplicative padding)."""
    y = y_top * pad_mult
    tick_low = y / tick_mult
    ax.plot([x1, x1], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x2, x2], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x1, x2], [y, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.text((x1 + x2) / 2, y * tick_mult, label, ha="center", va="bottom",
             fontsize=9.5, fontfamily=FONT, color=LABEL_COLOR)
    return y * tick_mult * 1.15


def conc_panel(ax, rng):
    ax.set_facecolor(BG)
    # Grouped by fraction (TEV pair, then EEV pair) so each Healthy-vs-
    # Diabetic significance bracket sits directly over the pair it tests,
    # rather than spanning across the other fraction's bar.
    x = np.array([0, 1, 2.6, 3.6])
    data4 = [tev[0], tev[1], eev[0], eev[1]]
    bar_colors4 = ["#2a78d6", "#e34948", "#2a78d6", "#e34948"]
    dot_colors4 = ["#cde2fb", "#f5b8b8", "#cde2fb", "#f5b8b8"]
    labels4 = ["TEV\nHealthy", "TEV\nDiabetic", "EEV\nHealthy", "EEV\nDiabetic"]

    means4 = [d.mean() for d in data4]
    sds4 = [d.std(ddof=1) for d in data4]
    ns4 = [len(d) for d in data4]

    ax.bar(x, means4, color=bar_colors4, edgecolor=bar_colors4, linewidth=0.8,
           width=0.55, zorder=3, alpha=0.80)

    # Upper-half SD whisker only (log axis - a symmetric linear SD whisker
    # is not meaningful here; see module docstring)
    for xi, m, s in zip(x, means4, sds4):
        upper = m + s
        ax.plot([xi, xi], [m, upper], color=LABEL_COLOR, linewidth=1.4, zorder=4)
        ax.plot([xi - 0.08, xi + 0.08], [upper, upper], color=LABEL_COLOR, linewidth=1.4, zorder=4)

    for xi, vals, dcolor in zip(x, data4, dot_colors4):
        jit = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(xi + jit, vals, color=dcolor, edgecolors="white",
                   linewidths=0.5, s=42, zorder=5)

    # ns brackets: TEV Healthy-vs-Diabetic (p=1.000), EEV Healthy-vs-Diabetic (p=0.596)
    tev_top = max(np.max(data4[0]), np.max(data4[1]), means4[0] + sds4[0], means4[1] + sds4[1])
    eev_top = max(np.max(data4[2]), np.max(data4[3]), means4[2] + sds4[2], means4[3] + sds4[3])
    label_top = sig_bracket_log(ax, x[0], x[1], tev_top, "ns")
    sig_bracket_log(ax, x[2], x[3], eev_top, "ns")

    ax.set_axisbelow(True)
    ax.set_yscale("log")
    ax.set_ylim(top=max(label_top, ax.get_ylim()[1]))
    ax.yaxis.grid(True, which="major", color=GRID, linewidth=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{l}\n(n={n})" for l, n in zip(labels4, ns4)],
                        fontsize=9.5, fontfamily=FONT, color=TITLE_COLOR)
    ax.set_ylabel("Particles per mL Plasma (log scale)", fontsize=10, fontfamily=FONT,
                  color=LABEL_COLOR, labelpad=8)
    ax.tick_params(axis="y", labelsize=9, labelcolor=LABEL_COLOR)
    ax.tick_params(axis="x", bottom=False)
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(SPINE)
    ax.spines["left"].set_linewidth(0.8)
    ax.set_title("EV Secretion - TEV vs. EEV, Healthy vs. Diabetic",
                 fontsize=11, fontweight="bold", fontfamily=FONT,
                 color=TITLE_COLOR, pad=12, loc="center")
    add_panel_letter(ax, "A")


def sig_bracket_linear(ax, x1, x2, y_top, label, pad=6, tick=3):
    """Draw a significance bracket for a linear-scale axis (additive padding)."""
    y = y_top + pad
    tick_low = y - tick
    ax.plot([x1, x1], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x2, x2], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x1, x2], [y, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.text((x1 + x2) / 2, y + 1, label, ha="center", va="bottom",
             fontsize=9.5, fontfamily=FONT, color=LABEL_COLOR)
    return y + 8


def pct_panel(ax, rng):
    ax.set_facecolor(BG)
    x = np.array([0, 1])
    means = [d.mean() for d in pct]
    sds = [d.std(ddof=1) for d in pct]
    ns = [len(d) for d in pct]
    colors = ["#2a78d6", "#e34948"]
    dcolors = ["#cde2fb", "#f5b8b8"]

    ax.bar(x, means, yerr=sds,
           color=colors, edgecolor=colors, linewidth=0.8,
           error_kw=dict(elinewidth=1.4, capsize=5, capthick=1.4, ecolor=LABEL_COLOR),
           width=0.55, zorder=3, alpha=0.80)
    for i, (vals, dcolor) in enumerate(zip(pct, dcolors)):
        jit = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jit, vals, color=dcolor,
                   edgecolors="white", linewidths=0.5, s=42, zorder=5)

    pct_top = max(np.max(pct[0]), np.max(pct[1]), means[0] + sds[0], means[1] + sds[1])
    label_top = sig_bracket_linear(ax, x[0], x[1], pct_top, "ns")

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, ns)],
                        fontsize=10, fontfamily=FONT, color=TITLE_COLOR)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, max(100, label_top))
    ax.set_ylabel("%EEVs from Total EVs", fontsize=10, fontfamily=FONT,
                  color=LABEL_COLOR, labelpad=8)
    ax.tick_params(axis="y", labelsize=9, labelcolor=LABEL_COLOR)
    ax.tick_params(axis="x", bottom=False)
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(SPINE)
    ax.spines["left"].set_linewidth(0.8)
    ax.set_title("%EEVs from Total EVs - Healthy vs. Diabetic",
                 fontsize=11, fontweight="bold", fontfamily=FONT,
                 color=TITLE_COLOR, pad=12, loc="center")
    add_panel_letter(ax, "B")


rng = np.random.default_rng(RNG_SEED)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.8))
fig.patch.set_facecolor(BG)
conc_panel(axes[0], rng)
rng2 = np.random.default_rng(RNG_SEED)
pct_panel(axes[1], rng2)

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"EV_Secretion_and_Percent_EEV.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)

print("Saved: EV_Secretion_and_Percent_EEV.png/.svg")
