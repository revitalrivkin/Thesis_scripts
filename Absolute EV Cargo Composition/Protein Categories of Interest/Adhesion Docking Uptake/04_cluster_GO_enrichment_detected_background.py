"""
Adhesion / Docking / Uptake - GO BP enrichment of the 3 direction clusters
==============================================================================
Option A: the background (universe) is the DETECTED EV proteome
(EVs_df_filtered.csv, 2,866 proteins), not the whole human genome.

Clusters (same definition as plot_adhesion_docking_uptake_clusters.py):
DEPs (p < 0.05 and |log2FC| > 0.85 in HG and/or HL vs Normal) that belong to the
Adhesion / Docking / Uptake category, split by direction of change:
  both > 0 = Consistently Up; both < 0 = Consistently Down; otherwise Discordant.

Method: one-sided Fisher's exact test per GO BP term, annotations propagated up
the GO graph (is_a + part_of, from go-basic.obo), terms with 5-500 background
proteins tested, Benjamini-Hochberg FDR within each cluster.
Only proteins with at least one BP annotation are used (query and background).

Output: cluster_GO_enrichment_detected_background.csv (all terms with >= 2 hits)
"""

import gzip
import os
from collections import defaultdict

import pandas as pd
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "..", "category_genesets.csv")
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"

SIG, FC = 0.05, 0.85
CATEGORY = "Adhesion / Docking / Uptake"
MIN_BG, MAX_BG = 5, 500

# -- GO graph (is_a + part_of) and term names ------------------------------------
parents, names = defaultdict(set), {}
cur = None
with open(os.path.join(GO_DATA, "go-basic.obo"), encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if line in ("[Term]", "[Typedef]"):
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif cur and line.startswith("name: "):
            names[cur] = line[6:]
        elif cur and line.startswith("is_a: "):
            parents[cur].add(line.split()[1])
        elif cur and line.startswith("relationship: part_of "):
            parents[cur].add(line.split()[2])

_anc = {}


def ancestors(t):
    if t in _anc:
        return _anc[t]
    out = {t}
    for p in parents.get(t, ()):
        out |= ancestors(p)
    _anc[t] = out
    return out


# -- gene -> propagated BP terms ------------------------------------------------
direct = defaultdict(set)
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10 or c[8] != "P" or "NOT" in c[3]:
            continue
        direct[c[2]].add(c[4])


def terms_of(gene_field):
    t = set()
    for g in str(gene_field).split(";"):
        for go in direct.get(g, ()):
            t |= ancestors(go)
    return t


# -- data -------------------------------------------------------------------------
df = pd.read_csv(DATA)
gene_terms = {g: terms_of(g) for g in df["Genes"].dropna()}
background = {g for g, t in gene_terms.items() if t}

d = df[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
        "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d.columns = ["Genes", "p_HG", "HG", "p_HL", "HL"]
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")
dep = ((d["p_HG"] < SIG) & (d["HG"].abs() > FC)) | ((d["p_HL"] < SIG) & (d["HL"].abs() > FC))
cat_genes = set(pd.read_csv(GENESETS).query("Category == @CATEGORY")["Gene"])
d = d[dep & d.index.isin(cat_genes)]

clusters = {
    "Consistently Up": d[(d["HG"] > 0) & (d["HL"] > 0)],
    "Consistently Down": d[(d["HG"] < 0) & (d["HL"] < 0)],
    "Discordant": d[~(((d["HG"] > 0) & (d["HL"] > 0)) | ((d["HG"] < 0) & (d["HL"] < 0)))],
}

# term -> background genes
term_bg = defaultdict(set)
for g in background:
    for t in gene_terms[g]:
        term_bg[t].add(g)
tested = {t: s for t, s in term_bg.items() if MIN_BG <= len(s) <= MAX_BG}
N = len(background)
print(f"Background: {N} detected proteins with a BP annotation; {len(tested)} testable GO terms")

rows = []
for cname, sub in clusters.items():
    q = set(sub.index) & background
    print(f"{cname}: {len(sub)} DEPs, {len(q)} annotated")
    res = []
    for t, bg in tested.items():
        k = len(q & bg)
        if k < 2:
            continue
        m = len(bg)
        p = fisher_exact([[k, len(q) - k], [m - k, N - m - (len(q) - k)]], alternative="greater")[1]
        res.append([cname, t, names.get(t, t), k, len(q), m, N, p])
    res.sort(key=lambda r: r[-1])
    n = len(res)
    prev = 1.0
    adj = [None] * n
    for rank in range(n, 0, -1):
        prev = min(prev, res[rank - 1][-1] * n / rank)
        adj[rank - 1] = prev
    for r, a in zip(res, adj):
        rows.append(r + [a])

out = pd.DataFrame(rows, columns=["Cluster", "GO_ID", "Term", "k_in_cluster", "n_cluster_annotated",
                                  "n_term_in_background", "N_background", "p", "FDR_BH"])
out.to_csv(os.path.join(HERE, "cluster_GO_enrichment_detected_background.csv"), index=False)

for cname in clusters:
    print(f"\n--- {cname} (top 8 by p) ---")
    print(out[out.Cluster == cname].head(8)[["Term", "k_in_cluster", "n_term_in_background", "p", "FDR_BH"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))
