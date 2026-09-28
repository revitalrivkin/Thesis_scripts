"""
Targeted Enrichment of Adhesion, Docking & Uptake Proteins Among DEPs
======================================================================
Step 3 (corrected 2026-08-19): the STRING network input is the FULL
downregulated DEP list per condition (not just the curated-set overlap).
The curated adhesion/docking/uptake genes are a separate "highlight" list —
paste the full list into STRING, then use the highlight list to mark/color
those specific nodes once the network is built.

Input:  ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
        ../Protein Categories of Interest/Adhesion Docking Uptake/targeted_enrichment_results.csv
Output (2026-08-21, updated: written into the direction-specific subfolders
— HG downregulated/, HL downregulated/ — created alongside HG upregulated/
and HL upregulated/, which are empty until that direction is built):
        HG downregulated/HG_down_ALL_DEPs_STRING_input.txt   (full network input)
        HG downregulated/HG_down_adhesion_docking_uptake_HIGHLIGHT.txt (subset to mark)
        HL downregulated/HL_down_ALL_DEPs_STRING_input.txt
        HL downregulated/HL_down_adhesion_docking_uptake_HIGHLIGHT.txt
"""

import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")

SIG_THRESHOLD, LOG2FC_THRESHOLD = 0.05, 0.85

df = pd.read_csv(DATA)

def get_down(pval_col, fc_col):
    d = df[["Genes", pval_col, fc_col]].dropna()
    return sorted(d.loc[(d[pval_col] < SIG_THRESHOLD) & (d[fc_col] < -LOG2FC_THRESHOLD), "Genes"])

hg_down_all = get_down("Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal")
hl_down_all = get_down("Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal")

results = pd.read_csv(os.path.join(HERE, "..", "Protein Categories of Interest", "Adhesion Docking Uptake",
                                    "targeted_enrichment_results.csv"))
hg_highlight = sorted(g.strip() for g in results.loc[results["DEP_category"] == "HG down", "genes"].iloc[0].split(","))
hl_highlight = sorted(g.strip() for g in results.loc[results["DEP_category"] == "HL down", "genes"].iloc[0].split(","))

TARGETS = {
    "HG": {"folder": "HG downregulated", "all_deps": hg_down_all, "highlight": hg_highlight},
    "HL": {"folder": "HL downregulated", "all_deps": hl_down_all, "highlight": hl_highlight},
}

for label, info in TARGETS.items():
    folder = os.path.join(HERE, info["folder"])
    os.makedirs(folder, exist_ok=True)

    all_path = os.path.join(folder, f"{label}_down_ALL_DEPs_STRING_input.txt")
    with open(all_path, "w", encoding="utf-8") as f:
        f.write("\n".join(info["all_deps"]))
    print(f"{label} down — full STRING input: {len(info['all_deps'])} genes -> {all_path}")

    hl_path = os.path.join(folder, f"{label}_down_adhesion_docking_uptake_HIGHLIGHT.txt")
    with open(hl_path, "w", encoding="utf-8") as f:
        f.write("\n".join(info["highlight"]))
    print(f"{label} down — highlight subset: {len(info['highlight'])} genes -> {hl_path}")
