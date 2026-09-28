"""
4.2.3.5 Enriched/Depleted by Pathology — Cross-Reference Against Categories of Interest, HG
================================================================================
Checks which of the HG pathology-depleted cargo proteins (from
01_build_pathology_depleted_cargo_HG.py) belong to the 5 finalized
categories of interest (category_genesets.csv), using the same priority-
resolution order as the 4.2.2 figures.

Input:  01_pathology_depleted_cargo_HG.xlsx
        ../../../../Absolute EV Cargo Composition/Protein Categories of Interest/category_genesets.csv
Output: printed summary + 02_pathology_depleted_categories_HG.xlsx
"""

import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
COND = "HG"
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

df = pd.read_excel(os.path.join(HERE, f"01_pathology_depleted_cargo_{COND}.xlsx"))
df["Category"] = df["Gene"].apply(categorize)

print(f"=== {COND} pathology-depleted cargo (n={len(df)}) ===")
counts = df["Category"].value_counts()
for cat in CATEGORY_PRIORITY + ["Other"]:
    n = counts.get(cat, 0)
    genes = sorted(df.loc[df["Category"] == cat, "Gene"]) if n else []
    print(f"  {cat}: {n}" + (f" — {', '.join(genes)}" if n else ""))

out_path = os.path.join(HERE, f"02_pathology_depleted_categories_{COND}.xlsx")
df.to_excel(out_path, index=False)
print("\nSaved:", out_path)
