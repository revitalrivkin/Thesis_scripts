"""
4.1.1 Donor Cohort Characteristics - Summary Figure (Age, BMI, HbA1c, Sex)
================================================================================
4-panel (2x2) summary of the donor cohort's clinical characteristics, Healthy
vs. Diabetic, built to accompany Table 1 (see ../Table1_Donor_Cohort_
Characteristics_MedianRange.xlsx).

Restyled 2026-09-16 to match the established house style from
`NTA\Particle concentration\plot_EEV_concentration_H_vs_D.py`: bar + mean+-SD
error bars + jittered per-donor scatter dots, same colors/fonts/spines, for
visual consistency across the thesis's donor-comparison figures.

Panels:
  - Age, BMI, HbA1c: bar (mean) + SD error bar + every donor as a jittered
    dot (dot color = pale group tint, not per-donor).
  - HbA1c panel additionally marks the 6.5% diagnostic threshold (dashed
    line) -- visually confirms the two groups are cleanly separated by the
    actual criterion defining them (Section 3.5.1), not just numerically
    claimed.
  - Sex: bar chart, % male per group (categorical, no dots).

Input:  ../../Clinical parameters of donor cohort.csv
Output: Donor_Cohort_Summary.png/.svg (one level up, in Figures/)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT = "Arial"
BG = "#fcfcfb"
GRID = "#e1e0d9"
SPINE = "#c3c2b7"
TITLE_COLOR = "#0b0b0b"
LABEL_COLOR = "#52514e"
BAR_COLORS = {"Healthy": "#2a78d6", "Diabetic": "#e34948"}
DOT_COLORS = {"Healthy": "#cde2fb", "Diabetic": "#f5b8b8"}

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.dirname(HERE)
DATA = os.path.join(FIG_DIR, "..", "Clinical parameters of donor cohort.csv")

RNG_SEED = 42

df = pd.read_csv(DATA)
df.columns = [c.strip() for c in df.columns]
df["Group_label"] = df["Group"].map({"Control": "Healthy", "Diabetic": "Diabetic"})

for c in ["Age", "BMI", "HbA₁c (%)"]:
    df[c] = pd.to_numeric(df[c].astype(str).str.strip().replace("X", np.nan), errors="coerce")

groups = ["Healthy", "Diabetic"]

fig, axes = plt.subplots(2, 2, figsize=(10, 9))
fig.patch.set_facecolor(BG)
rng = np.random.default_rng(RNG_SEED)

def add_panel_letter(ax, letter):
    ax.text(-0.14, 1.12, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", color=TITLE_COLOR, fontfamily=FONT,
            va="bottom", ha="left")

def bar_dot_panel(ax, col, title, ylabel, threshold=None, threshold_label=None, letter=None):
    ax.set_facecolor(BG)
    x = np.array([0, 1])
    data = [df.loc[df["Group_label"] == g, col].dropna().values for g in groups]
    means = [d.mean() for d in data]
    sds = [d.std(ddof=1) for d in data]
    ns = [len(d) for d in data]

    ax.bar(x, means, yerr=sds,
           color=[BAR_COLORS[g] for g in groups], edgecolor=[BAR_COLORS[g] for g in groups],
           linewidth=0.8, error_kw=dict(elinewidth=1.4, capsize=5, capthick=1.4, ecolor="#52514e"),
           width=0.55, zorder=3, alpha=0.80)

    for i, (g, vals) in enumerate(zip(groups, data)):
        jit = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jit, vals, color=DOT_COLORS[g],
                   edgecolors="white", linewidths=0.5, s=42, zorder=5)

    if threshold is not None:
        ax.axhline(threshold, color=LABEL_COLOR, linestyle=(0, (4, 3)), linewidth=1.2, zorder=1)
        ax.text(1.45, threshold, threshold_label, fontsize=8.5, color=LABEL_COLOR,
                fontfamily=FONT, va="center", ha="left")

    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, ns)],
                        fontsize=10, fontfamily=FONT, color=TITLE_COLOR)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel(ylabel, fontsize=10, fontfamily=FONT, color=LABEL_COLOR, labelpad=8)
    ax.tick_params(axis="y", labelsize=9, labelcolor=LABEL_COLOR)
    ax.tick_params(axis="x", bottom=False)
    for lbl in ax.get_yticklabels():
        lbl.set_fontfamily(FONT)
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(SPINE)
    ax.spines["left"].set_linewidth(0.8)
    ax.set_title(title, fontsize=11, fontweight="bold", fontfamily=FONT, color=TITLE_COLOR, pad=12, loc="center")
    if letter:
        add_panel_letter(ax, letter)

bar_dot_panel(axes[0, 0], "Age", "Age - Healthy vs. Diabetic Donors", "Years", letter="A")
bar_dot_panel(axes[0, 1], "BMI", "BMI - Healthy vs. Diabetic Donors", "kg/m²", letter="B")
bar_dot_panel(axes[1, 0], "HbA₁c (%)", "HbA1c - Healthy vs. Diabetic Donors", "%",
              threshold=6.5, threshold_label="Diagnostic\nthreshold (6.5%)", letter="C")

# Sex panel: % male per group (categorical, no dots/error bars)
ax_sex = axes[1, 1]
ax_sex.set_facecolor(BG)
pct_male, n_labels, ns_sex = [], [], []
for g in groups:
    sub = df.loc[df["Group_label"] == g, "Sex"].dropna()
    n_total = len(sub)
    n_male = (sub.astype(str).str.strip() == "M").sum()
    pct_male.append(100 * n_male / n_total if n_total else 0)
    n_labels.append(f"{n_male}/{n_total} male")
    ns_sex.append(n_total)

ax_sex.bar([0, 1], pct_male, color=[BAR_COLORS[g] for g in groups],
           edgecolor=[BAR_COLORS[g] for g in groups], linewidth=0.8,
           width=0.55, zorder=3, alpha=0.80)
for i, (pct, lbl) in enumerate(zip(pct_male, n_labels)):
    ax_sex.text(i, pct + 3, f"{pct:.0f}%\n({lbl})", ha="center", va="bottom",
                fontsize=9, color=TITLE_COLOR, fontfamily=FONT)
ax_sex.set_axisbelow(True)
ax_sex.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax_sex.set_xticks([0, 1])
ax_sex.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(groups, ns_sex)],
                        fontsize=10, fontfamily=FONT, color=TITLE_COLOR)
ax_sex.set_xlim(-0.5, 1.5)
ax_sex.set_ylim(0, 100)
ax_sex.set_ylabel("% male", fontsize=10, fontfamily=FONT, color=LABEL_COLOR, labelpad=8)
ax_sex.tick_params(axis="y", labelsize=9, labelcolor=LABEL_COLOR)
ax_sex.tick_params(axis="x", bottom=False)
for lbl in ax_sex.get_yticklabels():
    lbl.set_fontfamily(FONT)
for spine in ["top", "right", "bottom"]:
    ax_sex.spines[spine].set_visible(False)
ax_sex.spines["left"].set_color(SPINE)
ax_sex.spines["left"].set_linewidth(0.8)
ax_sex.set_title("Sex - Healthy vs. Diabetic Donors", fontsize=11, fontweight="bold", fontfamily=FONT,
                  color=TITLE_COLOR, pad=12, loc="center")
add_panel_letter(ax_sex, "D")

fig.suptitle("Donor Cohort Characteristics", fontsize=15, fontweight="bold",
             color=TITLE_COLOR, y=0.99, fontfamily=FONT)

fig.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"Donor_Cohort_Summary.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print("Saved: Donor_Cohort_Summary.png/.svg")
