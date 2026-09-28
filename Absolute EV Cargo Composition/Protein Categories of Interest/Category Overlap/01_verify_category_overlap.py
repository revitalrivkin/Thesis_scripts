"""
Category overlap - re-verification of counts (2026-09-26)
==============================================================================
Recomputes, from EVs_df_filtered.csv + category_genesets.csv, how many detected
proteins belong to 2+ of the 5 categories, and the pairwise overlaps, for (a) all
detected proteins and (b) DEPs only (HG and/or HL; p < 0.05 and |log2FC| > 0.85).
Also counts, per volcano DEP set, how many DEPs are in no category.

Output: category_overlap_check_report.txt
"""

import itertools
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")
SIG, FC = 0.05, 0.85

df = pd.read_csv(DATA)
detected = set(df["Genes"].dropna())
gs = pd.read_csv(GENESETS)
cats = ["Adhesion / Docking / Uptake", "Cardioprotective", "ECM Organization", "Immune-evasion", "EV Biogenesis"]
sets = {c: set(gs.loc[gs["Category"] == c, "Gene"]) & detected for c in cats}

d = df[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
        "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d.columns = ["Genes", "pHG", "HG", "pHL", "HL"]
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")
deps = {
    "HG up": set(d.index[(d.pHG < SIG) & (d.HG > FC)]), "HG down": set(d.index[(d.pHG < SIG) & (d.HG < -FC)]),
    "HL up": set(d.index[(d.pHL < SIG) & (d.HL > FC)]), "HL down": set(d.index[(d.pHL < SIG) & (d.HL < -FC)]),
}
deps["any DEP"] = set().union(*[v for k, v in deps.items()])

rep = []


def block(label, universe):
    rows = [f"=== {label} (n = {len(universe)}) ==="]
    counts = {g: sum(g in sets[c] for c in cats) for g in universe}
    n1 = sum(v >= 1 for v in counts.values())
    n2 = sum(v >= 2 for v in counts.values())
    rows.append(f"in >=1 category: {n1} ({n1 / len(universe):.1%}); in >=2: {n2} "
                f"({(n2 / n1 if n1 else 0):.1%} of categorized); in none: {len(universe) - n1}")
    for c in cats:
        rows.append(f"  {c}: {len(sets[c] & universe)}")
    for a, b in itertools.combinations(cats, 2):
        rows.append(f"  overlap {a} x {b}: {len(sets[a] & sets[b] & universe)}")
    return rows


rep += block("All detected proteins", detected)
for k, v in deps.items():
    rep += block(f"DEPs: {k}", v & detected)
open(os.path.join(HERE, "category_overlap_check_report.txt"), "w", encoding="utf-8").write("\n".join(rep))
print("\n".join(rep))
