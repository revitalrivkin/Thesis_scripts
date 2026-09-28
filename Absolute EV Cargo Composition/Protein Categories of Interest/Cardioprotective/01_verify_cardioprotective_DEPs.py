"""
Cardioprotective category - verification of DEP counts and enrichment (2026-09-26)
==============================================================================
1. Recomputes, from EVs_df_filtered.csv + category_genesets.csv, how many
   Cardioprotective proteins are detected / are DEPs, overall and per
   subcategory (subcategories overlap, so per-subcategory counts are
   non-exclusive), and their direction in HG and HL.
2. Tests whether DEPs (up / down, HG and HL) are over- or under-represented in
   the category and in each subcategory relative to the DETECTED proteome
   (one-sided Fisher's exact, both tails reported; same fair background chosen
   for the Adhesion category), with Benjamini-Hochberg across all tests shown.

DEP definition = volcano: Student's t-test p < 0.05 and |log2FC| > 0.85.

Output: cardioprotective_DEP_check.csv, cardioprotective_DEP_check_report.txt
"""

import os

import pandas as pd
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")
SIG, FC = 0.05, 0.85
SUBS = ["Ca2+ homeostasis", "AMPK signaling", "Antioxidants/Redox", "Glycolysis", "UPR", "HSPs"]

df = pd.read_csv(DATA)
detected = set(df["Genes"].dropna())
N = len(detected)
gs = pd.read_csv(GENESETS)
cardio = gs[gs["Category"] == "Cardioprotective"]
sets = {"Cardioprotective (all)": set(cardio["Gene"]) & detected}
for s in SUBS:
    sets[s] = set(cardio.loc[cardio["Subcategory"] == s, "Gene"]) & detected

d = df[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
        "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d.columns = ["Genes", "pHG", "HG", "pHL", "HL"]
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")

dep = {}
for comp, p, f in (("HG", "pHG", "HG"), ("HL", "pHL", "HL")):
    dep[f"{comp} up"] = set(d.index[(d[p] < SIG) & (d[f] > FC)])
    dep[f"{comp} down"] = set(d.index[(d[p] < SIG) & (d[f] < -FC)])
any_dep = set().union(*dep.values())

rep = [f"Detected proteome: {N} proteins; DEPs in HG and/or HL: {len(any_dep)}"]
rep.append("\n=== Detected / DEP counts (subcategories overlap) ===")
for name, g in sets.items():
    both = {x for x in g if x in dep["HG up"] | dep["HG down"]} & {x for x in g if x in dep["HL up"] | dep["HL down"]}
    rep.append(f"{name}: detected {len(g)}, DEP in HG or HL {len(g & any_dep)}, DEP in both {len(both)}"
               f" | HG up {len(g & dep['HG up'])} down {len(g & dep['HG down'])}"
               f" | HL up {len(g & dep['HL up'])} down {len(g & dep['HL down'])}")

rows = []
for name, g in sets.items():
    for dname, dset in dep.items():
        q = dset  # all DEPs in that set are in the detected universe
        k = len(g & q)
        n_q, m = len(q), len(g)
        tbl = [[k, n_q - k], [m - k, N - m - (n_q - k)]]
        p_over = fisher_exact(tbl, alternative="greater")[1]
        p_under = fisher_exact(tbl, alternative="less")[1]
        rows.append({"Set": name, "DEP_set": dname, "n_DEPs": n_q, "k_in_set": k, "n_set_detected": m,
                     "observed_frac": round(k / n_q, 4) if n_q else None, "expected_frac": round(m / N, 4),
                     "p_enriched": p_over, "p_depleted": p_under})
res = pd.DataFrame(rows)


def bh(p):
    order = sorted(range(len(p)), key=lambda i: p[i])
    n, adj, prev = len(p), [None] * len(p), 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, p[i] * n / rank)
        adj[i] = prev
    return adj


res["q_enriched_BH"] = bh(list(res["p_enriched"]))
res["q_depleted_BH"] = bh(list(res["p_depleted"]))
res.to_csv(os.path.join(HERE, "cardioprotective_DEP_check.csv"), index=False)

rep.append("\n=== Enrichment/depletion vs detected proteome (Fisher, one-sided each way) ===")
for _, r in res.iterrows():
    rep.append(f"{r['Set']} | {r['DEP_set']}: {r['k_in_set']}/{r['n_DEPs']} ({r['observed_frac']:.1%} obs vs "
               f"{r['expected_frac']:.1%} exp) p_enr={r['p_enriched']:.3g} (q={r['q_enriched_BH']:.3g}) "
               f"p_dep={r['p_depleted']:.3g} (q={r['q_depleted_BH']:.3g})")

with open(os.path.join(HERE, "cardioprotective_DEP_check_report.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(rep))
print("\n".join(rep))
