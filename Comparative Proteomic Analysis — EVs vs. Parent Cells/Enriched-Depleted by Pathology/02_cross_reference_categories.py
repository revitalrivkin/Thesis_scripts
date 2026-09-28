"""
4.2.3.5 Enriched/Depleted by Pathology — Cross-Reference Against Categories of Interest
================================================================================
Checks which of the pathology-depleted cargo proteins (HG, HL — from
01_build_pathology_depleted_cargo_list.py) belong to the 5 finalized
categories of interest (category_genesets.csv), using the same priority-
resolution order as the 4.2.2 figures (most specific category wins for
genes matching more than one).

Input:  01_pathology_depleted_cargo_HG.xlsx, 01_pathology_depleted_cargo_HL.xlsx
        ../../Absolute EV Cargo Composition/Protein Categories of Interest/category_genesets.csv
Output: printed summary + 02_pathology_depleted_categories.xlsx
"""

import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GENESETS = (r"C:\Users\רויטל\Desktop\ISF\Thesis\Results\Absolute EV Cargo Composition"
            r"\Protein Categories of Interest\category_genesets.csv")

CATEGORY_PRIORITY = [
    "EV Biogenesis",
    "ECM Organization",
    "Immune-evasion",
    "Cardioprotective",
    "Adhesion / Docking / Uptake",
]

genesets = pd.read_csv(GENESETS)
gene_to_cats = genesets.groupby("Gene")["Category"].apply(set).to_dict()

def categorize(gene):
    cats = gene_to_cats.get(gene, set())
    for cat in CATEGORY_PRIORITY:
        if cat in cats:
            return cat
    return "Other"

all_rows = []
for cond in ["HG", "HL"]:
    df = pd.read_excel(os.path.join(HERE, f"01_pathology_depleted_cargo_{cond}.xlsx"))
    df["Category"] = df["Gene"].apply(categorize)
    df["Condition"] = cond
    print(f"\n=== {cond} pathology-depleted cargo (n={len(df)}) ===")
    counts = df["Category"].value_counts()
    for cat in CATEGORY_PRIORITY + ["Other"]:
        n = counts.get(cat, 0)
        genes = sorted(df.loc[df["Category"] == cat, "Gene"]) if n else []
        print(f"  {cat}: {n}" + (f" — {', '.join(genes)}" if n else ""))
    all_rows.append(df)

out_df = pd.concat(all_rows, ignore_index=True)
out_path = os.path.join(HERE, "02_pathology_depleted_categories.xlsx")
out_df.to_excel(out_path, index=False)
print("\nSaved:", out_path)
