"""
Specificity-Weighting Check — Categories of Interest
==============================================================================
Exploratory methodology check (not a figure). For each detected protein in
each category, computes:
  - N_matched_terms: how many of the category's curated top-level GO terms
    the gene matches (via descendant propagation) — multi-term = more
    confidently "core" to the category
  - min_distance: shortest path (in GO-graph hops) from the gene's own
    DIRECT annotation to the nearest matched curated top-level term —
    0 = gene is directly annotated to the curated term itself ("Direct"),
    larger = only connected via a deep, more tangential descendant term

Genes are then bucketed: Direct (distance=0) / Multi-term (matches 2+
curated terms) / Single-distant (matches exactly 1 term, only via a
descendant — the most "loose/tangential" bucket).

Only genes actually detected in the filtered EV dataset are checked (the
"all detected" universe used by the categories heatmap), since that's the
practically relevant set.

Output: specificity_weighting_report.csv (per-gene detail)
        printed summary table (per category)
"""

import gzip
import os
from collections import defaultdict, deque

import pandas as pd

GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(BASE, "Absolute EV Cargo Composition", "Volcano", "Volcano with categories", "category_genesets.csv")
OUT = os.path.dirname(os.path.abspath(__file__))

# Category -> curated top-level GO terms (same as the reference table)
CATEGORY_TERMS = {
    "Adhesion / Docking / Uptake": ["GO:0007155", "GO:0006897", "GO:0006898", "GO:0016192",
                                      "GO:0006906", "GO:0006904", "GO:0140029", "GO:0032940"],
    "EV Biogenesis": ["GO:0140112", "GO:0036258"],
    "Immune-evasion": ["GO:0050777", "GO:0045916", "GO:0050765", "GO:0002698", "GO:0001915", "GO:0045953"],
    "ECM Organization": ["GO:0030198"],
    "Cardioprotective": ["GO:0006874", "GO:0070588", "GO:0031669", "GO:0006112", "GO:0045454",
                          "GO:0098869", "GO:0006979", "GO:0006096", "GO:0034976", "GO:0006986",
                          "GO:0034975", "GO:0006487", "GO:0006457", "GO:0034605"],
}


def parse_obo(file_path):
    parents = defaultdict(list)
    terms, current = {}, None
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line in ["[Term]", "[Typedef]"]:
                if current:
                    terms[current["id"]] = current
                current = {}
            elif line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("is_a: "):
                parents[current["id"]].append(line.split()[1])
            elif line.startswith("relationship: part_of "):
                parents[current["id"]].append(line.split()[2])
    if current:
        terms[current["id"]] = current
    return terms, parents


print("Parsing go-basic.obo ...")
terms, parents = parse_obo(os.path.join(GO_DATA, "go-basic.obo"))

# Build children map (reverse of parents) for BFS distance from a top-level term
children = defaultdict(list)
for child, plist in parents.items():
    for p in plist:
        children[p].append(child)


def bfs_distances(root):
    """Distance (hops) from root to every descendant, including itself (distance 0)."""
    dist = {root: 0}
    q = deque([root])
    while q:
        cur = q.popleft()
        for ch in children.get(cur, []):
            if ch not in dist:
                dist[ch] = dist[cur] + 1
                q.append(ch)
    return dist


print("Parsing goa_human.gaf.gz (direct annotations) ...")
gene2go_direct = defaultdict(set)
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 10 or cols[8] != "P" or "NOT" in cols[3]:
            continue
        gene2go_direct[cols[2]].add(cols[4])

# -- Detected genes (practically relevant universe) -----------------------------
df = pd.read_csv(DATA)
detected = set(df["Genes"])

genesets = pd.read_csv(GENESETS)

rows = []
for category, top_terms in CATEGORY_TERMS.items():
    cat_genes = sorted(set(genesets.loc[genesets["Category"] == category, "Gene"]) & detected)
    print(f"\n{category}: {len(cat_genes)} detected genes to check")

    # distance maps for each curated top-level term
    term_dist_maps = {t: bfs_distances(t) for t in top_terms}

    for gene in cat_genes:
        direct_terms = gene2go_direct.get(gene, set())
        matched = []  # (curated_term, min_distance_for_this_term)
        for t, dmap in term_dist_maps.items():
            dists = [dmap[dt] for dt in direct_terms if dt in dmap]
            if dists:
                matched.append((t, min(dists)))
        if not matched:
            continue  # shouldn't happen, but guard
        n_matched_terms = len(matched)
        min_distance = min(d for _, d in matched)
        if min_distance == 0:
            bucket = "Direct"
        elif n_matched_terms >= 2:
            bucket = "Multi-term"
        else:
            bucket = "Single-distant"
        rows.append({"Gene": gene, "Category": category, "N_matched_curated_terms": n_matched_terms,
                      "Min_distance_hops": min_distance, "Confidence_bucket": bucket})

report = pd.DataFrame(rows)
report.to_csv(os.path.join(OUT, "specificity_weighting_report.csv"), index=False)

print("\n" + "=" * 70)
print("SUMMARY — Confidence bucket breakdown per category")
print("=" * 70)
summary = report.groupby(["Category", "Confidence_bucket"]).size().unstack(fill_value=0)
summary = summary.reindex(columns=["Direct", "Multi-term", "Single-distant"], fill_value=0)
summary["Total"] = summary.sum(axis=1)
for col in ["Direct", "Multi-term", "Single-distant"]:
    summary[f"{col}_%"] = (summary[col] / summary["Total"] * 100).round(1)
print(summary[["Total", "Direct", "Direct_%", "Multi-term", "Multi-term_%", "Single-distant", "Single-distant_%"]])

print("\nMedian distance-to-nearest-curated-term per category (higher = more tangential on average):")
print(report.groupby("Category")["Min_distance_hops"].median())

print("\nSaved: specificity_weighting_report.csv")
