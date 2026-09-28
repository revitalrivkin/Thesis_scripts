"""
Targeted Enrichment of Adhesion, Docking & Uptake Proteins Among DEPs
======================================================================
Step 2: Test whether the curated gene set (from 01_build_curated_geneset.py)
is enriched among each DEP category (Fisher's exact, one-sided), using the
same DEP definitions as the volcano plot (Proteome Overview/Volcano).

Input:  curated_geneset.csv (from step 1)
        ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv
Output: targeted_enrichment_results.csv
"""

import os
import pandas as pd
from scipy.stats import fisher_exact

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
OUT = os.path.dirname(os.path.abspath(__file__))
CURATED = os.path.join(OUT, "curated_geneset.csv")

SIG_THRESHOLD, LOG2FC_THRESHOLD = 0.05, 0.85

# Background: all human genes carrying a BP annotation (same universe used to
# build the curated set) — approximated here via the curated-set build step's
# gene2go_direct keys is not saved, so we re-derive background size directly
# from the GAF for consistency.
import gzip
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"
background = set()
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 10 or cols[8] != "P" or "NOT" in cols[3]:
            continue
        background.add(cols[2])

curated_genes = set(pd.read_csv(CURATED)["Gene"])
n_bg = len(background)
n_curated = len(curated_genes & background)
print(f"Background: {n_bg} genes. Curated set (in background): {n_curated} genes.")

# -- DEP lists (same logic/thresholds as the volcano plot) ---------------------
df = pd.read_csv(DATA)

def get_deps(pval_col, fc_col):
    d = df[["Genes", pval_col, fc_col]].dropna()
    up = set(d.loc[(d[pval_col] < SIG_THRESHOLD) & (d[fc_col] > LOG2FC_THRESHOLD), "Genes"])
    down = set(d.loc[(d[pval_col] < SIG_THRESHOLD) & (d[fc_col] < -LOG2FC_THRESHOLD), "Genes"])
    return up, down

hg_up, hg_down = get_deps("Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal")
hl_up, hl_down = get_deps("Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal")
dep_categories = {"HG up": hg_up, "HG down": hg_down, "HL up": hl_up, "HL down": hl_down}

# -- Fisher's exact test ---------------------------------------------------------
rows = []
for label, genes in dep_categories.items():
    query = {g for g in genes if g in background}
    n_query = len(query)
    overlap = sorted(query & curated_genes)
    k = len(overlap)
    table = [[k, n_query - k], [n_curated - k, n_bg - n_curated - (n_query - k)]]
    _, pval = fisher_exact(table, alternative="greater")
    rows.append({"DEP_category": label, "n_DEPs": n_query, "n_in_curated_set": k,
                 "fraction": round(k / n_query, 4) if n_query else 0,
                 "fisher_p": pval, "genes": ", ".join(overlap)})
    print(f"{label}: {k}/{n_query} DEPs in curated set (Fisher p={pval:.3g})")

pd.DataFrame(rows).to_csv(os.path.join(OUT, "targeted_enrichment_results.csv"), index=False)
print("\nSaved: targeted_enrichment_results.csv")
