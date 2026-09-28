"""
4.1.2 Plasma-EV Characterization - Particle Diameter Summary Figure
================================================================================
Single-panel figure: median particle diameter (nm) for TEV and EEV, Healthy
vs. Diabetic. Styled to match Figure 1 and Figure 2: Arial, off-white
background, bar + mean+-SD error bar + jittered per-donor scatter dots, bold
title, plain hyphens (no em dashes). Grouped by fraction (TEV pair, then EEV
pair - not by disease group), matching Figure 2 Panel A's layout, so each
Healthy-vs-Diabetic "ns" significance bracket sits directly over the pair it
tests. Single panel - no panel letter needed (per feedback_panel_lettering,
letters are for multi-panel figures only).

Unlike Figure 2 (concentration), diameter values fall within one order of
magnitude across all four groups (63-203 nm) - no extreme outlier forcing a
log axis or an asymmetric whisker, so this figure uses a standard linear
axis with a full two-sided mean+-SD error bar.

Statistics (recomputed here from the D-labeled CSV; matches note 16's
original values, computed under the old donor codes):
  - TEV, Healthy vs. Diabetic: Mann-Whitney (ns)
  - EEV, Healthy vs. Diabetic: Mann-Whitney (ns)
  - TEV vs. EEV, paired within donor, all 18 donors pooled (Healthy +
    Diabetic together): Wilcoxon signed-rank + paired t-test (ns). This is a
    whole-cohort paired comparison, not a per-group one, so it does not map
    onto a single bracket between two of the four bars - reported in the
    caption/text instead of as a figure annotation, to avoid implying a
    comparison between two specific bars that isn't the one actually tested.

Data source: ../../Plasma EV sample records 12.8 (D-labeled).csv

Output: Particle_Diameter_TEV_EEV_H_vs_D.png/.svg (one level up, in Figures/)
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


tev = [get_vals("Median TEV particle diameter [nm]", g) for g in groups]
eev = [get_vals("Median EEV particle diameter [nm]", g) for g in groups]

# ── Stats ────────────────────────────────────────────────────────────────
log_lines = []


def report(name, a, b):
    log_lines.append(f"\n{name}")
    for label, vals in zip(groups, [a, b]):
        sh = stats.shapiro(vals)
        log_lines.append(f"  {label}: n={len(vals)} mean={vals.mean():.1f} sd={vals.std(ddof=1):.1f} "
                          f"median={np.median(vals):.1f} Shapiro p={sh.pvalue:.4f}")
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    log_lines.append(f"  Mann-Whitney U={u.statistic:.1f} p={u.pvalue:.4f}")


report("TEV diameter (nm), Healthy vs. Diabetic", tev[0], tev[1])
report("EEV diameter (nm), Healthy vs. Diabetic", eev[0], eev[1])

# Paired TEV vs. EEV, all 18 donors pooled
all_tev = df["Median TEV particle diameter [nm]"].astype(float).values
all_eev = df["Median EEV particle diameter [nm]"].astype(float).values
diffs = all_tev - all_eev
w = stats.wilcoxon(all_tev, all_eev)
t = stats.ttest_rel(all_tev, all_eev)
sh_diff = stats.shapiro(diffs)
log_lines.append(f"\nTEV vs. EEV, paired within donor, all n={len(all_tev)} donors pooled")
log_lines.append(f"  median diff (TEV-EEV) = {np.median(diffs):.2f} nm; "
                  f"{(diffs > 0).sum()}/{len(diffs)} donors TEV>EEV")
log_lines.append(f"  Shapiro (differences) p={sh_diff.pvalue:.4f}")
log_lines.append(f"  Wilcoxon signed-rank W={w.statistic:.1f} p={w.pvalue:.4f}")
log_lines.append(f"  Paired t-test p={t.pvalue:.4f}")

with open(os.path.join(HERE, "particle_diameter_stats.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines))


def sig_bracket_linear(ax, x1, x2, y_top, label, pad=8, tick=3):
    y = y_top + pad
    tick_low = y - tick
    ax.plot([x1, x1], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x2, x2], [tick_low, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.plot([x1, x2], [y, y], color=LABEL_COLOR, linewidth=1.1, zorder=6)
    ax.text((x1 + x2) / 2, y + 1.5, label, ha="center", va="bottom",
             fontsize=9.5, fontfamily=FONT, color=LABEL_COLOR)
    return y + 12


rng = np.random.default_rng(RNG_SEED)
fig, ax = plt.subplots(figsize=(7.5, 5.8))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

# Grouped by fraction (TEV pair, then EEV pair), matching Figure 2 Panel A
x = np.array([0, 1, 2.6, 3.6])
data4 = [tev[0], tev[1], eev[0], eev[1]]
bar_colors4 = ["#2a78d6", "#e34948", "#2a78d6", "#e34948"]
dot_colors4 = ["#cde2fb", "#f5b8b8", "#cde2fb", "#f5b8b8"]
labels4 = ["TEV\nHealthy", "TEV\nDiabetic", "EEV\nHealthy", "EEV\nDiabetic"]

means4 = [d.mean() for d in data4]
sds4 = [d.std(ddof=1) for d in data4]
ns4 = [len(d) for d in data4]

ax.bar(x, means4, yerr=sds4,
       color=bar_colors4, edgecolor=bar_colors4, linewidth=0.8,
       error_kw=dict(elinewidth=1.4, capsize=5, capthick=1.4, ecolor=LABEL_COLOR),
       width=0.55, zorder=3, alpha=0.80)

for xi, vals, dcolor in zip(x, data4, dot_colors4):
    jit = rng.uniform(-0.12, 0.12, len(vals))
    ax.scatter(xi + jit, vals, color=dcolor, edgecolors="white",
               linewidths=0.5, s=42, zorder=5)

# ns brackets: TEV Healthy-vs-Diabetic, EEV Healthy-vs-Diabetic
tev_top = max(np.max(data4[0]), np.max(data4[1]), means4[0] + sds4[0], means4[1] + sds4[1])
eev_top = max(np.max(data4[2]), np.max(data4[3]), means4[2] + sds4[2], means4[3] + sds4[3])
label_top1 = sig_bracket_linear(ax, x[0], x[1], tev_top, "ns")
label_top2 = sig_bracket_linear(ax, x[2], x[3], eev_top, "ns")

ax.set_axisbelow(True)
ax.set_ylim(0, max(label_top1, label_top2) + 5)
ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
ax.set_xticks(x)
ax.set_xticklabels([f"{l}\n(n={n})" for l, n in zip(labels4, ns4)],
                    fontsize=9.5, fontfamily=FONT, color=TITLE_COLOR)
ax.set_ylabel("Median Particle Diameter (nm)", fontsize=10, fontfamily=FONT,
              color=LABEL_COLOR, labelpad=8)
ax.tick_params(axis="y", labelsize=9, labelcolor=LABEL_COLOR)
ax.tick_params(axis="x", bottom=False)
for lbl in ax.get_yticklabels():
    lbl.set_fontfamily(FONT)
for spine in ["top", "right", "bottom"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color(SPINE)
ax.spines["left"].set_linewidth(0.8)
ax.set_title("Particle Diameter - TEV vs. EEV, Healthy vs. Diabetic",
             fontsize=11.5, fontweight="bold", fontfamily=FONT,
             color=TITLE_COLOR, pad=12, loc="center")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(FIG_DIR, f"Particle_Diameter_TEV_EEV_H_vs_D.{ext}"), dpi=200,
                bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)

print("Saved: Particle_Diameter_TEV_EEV_H_vs_D.png/.svg")
