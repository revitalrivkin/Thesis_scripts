"""
EV Biogenesis category - verification of DEP counts and enrichment (2026-09-26)
==============================================================================
Same design as Cardioprotective/01_verify_cardioprotective_DEPs.py: counts
detected proteins and DEPs (HG/HL vs Normal; p < 0.05 and |log2FC| > 0.85) in
the EV Biogenesis category, lists the DEPs with their fold changes, and
tests over/under-representation of each DEP set against the DETECTED
proteome (2,866 proteins) with one-sided Fisher's exact tests (both tails).

Output: ev_biogenesis_DEP_check_report.txt, ev_biogenesis_DEPs_table.csv
"""

import os

import pandas as pd
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")
CATEGORY = "EV Biogenesis"
SIG, FC = 0.05, 0.85

df = pd.read_csv(DATA)
detected = set(df["Genes"].dropna())
N = len(detected)
gs = pd.read_csv(GENESETS)
cat = set(gs.loc[gs["Category"] == CATEGORY, "Gene"]) & detected

d = df[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
        "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d.columns = ["Genes", "pHG", "HG", "pHL", "HL"]
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")
d["DEP_HG"] = (d.pHG < SIG) & (d.HG.abs() > FC)
d["DEP_HL"] = (d.pHL < SIG) & (d.HL.abs() > FC)

sets = {
    "HG up": set(d.index[(d.pHG < SIG) & (d.HG > FC)]), "HG down": set(d.index[(d.pHG < SIG) & (d.HG < -FC)]),
    "HL up": set(d.index[(d.pHL < SIG) & (d.HL > FC)]), "HL down": set(d.index[(d.pHL < SIG) & (d.HL < -FC)]),
}
any_dep = set(d.index[d.DEP_HG | d.DEP_HL])
rep = [f"{CATEGORY}: detected {len(cat)}; DEP in HG or HL {len(cat & any_dep)}; DEP in both "
       f"{len(cat & set(d.index[d.DEP_HG & d.DEP_HL]))}; all-proteome DEPs {len(any_dep)}/{N} ({len(any_dep)/N:.1%})"]
for k, s in sets.items():
    kk, q, m = len(cat & s), len(s), len(cat)
    tbl = [[kk, q - kk], [m - kk, N - m - (q - kk)]]
    rep.append(f"{k}: {kk}/{q} ({kk/q:.1%} obs vs {m/N:.1%} exp) p_enr={fisher_exact(tbl, alternative='greater')[1]:.3g} "
               f"p_dep={fisher_exact(tbl, alternative='less')[1]:.3g}")

t = d.loc[sorted(cat & any_dep), ["HG", "pHG", "HL", "pHL", "DEP_HG", "DEP_HL"]].round(3)
t["mx"] = t[["HG", "HL"]].abs().max(axis=1)
t = t.sort_values("mx", ascending=False).drop(columns="mx")
t.to_csv(os.path.join(HERE, "ev_biogenesis_DEPs_table.csv"))
rep.append("\nDEPs (sorted by largest |log2FC|):\n" + t.to_string())
up_hg, dn_hg = (t.HG[t.DEP_HG] > 0).sum(), (t.HG[t.DEP_HG] < 0).sum()
up_hl, dn_hl = (t.HL[t.DEP_HL] > 0).sum(), (t.HL[t.DEP_HL] < 0).sum()
rep.append(f"\nDirection among EV Biogenesis DEPs: HG up {up_hg} down {dn_hg}; HL up {up_hl} down {dn_hl}")
open(os.path.join(HERE, "ev_biogenesis_DEP_check_report.txt"), "w", encoding="utf-8").write("\n".join(rep))
print("\n".join(rep))
