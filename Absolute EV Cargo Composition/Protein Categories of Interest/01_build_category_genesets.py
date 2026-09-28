"""
Categories of Interest — Build All Remaining GO-Derived Gene Sets
==============================================================================
Builds gene sets for EV Biogenesis, Immune-evasion, ECM Organization, and
Cardioprotective (union of all 6 subcategories), using Ophir's descendant-
gathering method (is_a + part_of) — same method as the Adhesion/Docking/
Uptake set (Absolute EV Cargo Composition/DEPs/Targeted Enrichment -
Adhesion Docking Uptake/01_build_curated_geneset.py), reused here rather
than rebuilt.

**EV Biogenesis marker supplement (2026-08-21)**: tested adding Cellular
Component GO terms (GO:0070062 extracellular exosome, GO:1903561, etc.) to
fix the category's known gap (0 DEP hits, missing CD63/CD82/CD9 — all
CC-annotated, not BP-annotated). That fix worked but WAY overcorrected —
GO:0070062 alone is a promiscuous CC term ("detected in some exosome MS
study"), ballooning the category from 55 to 2,346 genes, mostly generic
abundant proteins (ribosomes, proteasome, histones) with no real biogenesis
role. Reverted; instead the 4 canonical tetraspanin markers (CD9, CD63,
CD81, CD82) are hard-coded onto the original 55-gene BP-only set below —
precise fix, no reopening to thousands of unrelated proteins. See
[[12 - Thesis Writing]] for the full test writeup.

See [[12 - Thesis Writing]] for the terms table and coverage validation
against the PI's cardioprotective gene list.

Output: category_genesets.csv — Gene, Category, Subcategory (long format,
one row per gene-category membership; a gene can appear more than once if
it matches multiple categories)
"""

import gzip
import os
from collections import defaultdict

import pandas as pd

GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"
ADHESION_SET = (r"C:\Users\רויטל\Desktop\ISF\Thesis\Results\Absolute EV Cargo Composition"
                 r"\Protein Categories of Interest\Adhesion Docking Uptake\curated_geneset.csv")
OUT = os.path.dirname(os.path.abspath(__file__))

# Canonical exosome/EV surface markers, hard-coded onto EV Biogenesis's
# BP-only gene set — precise fix for known CC-vs-BP annotation gap, not a
# blanket CC-term expansion (see docstring above for why).
EV_BIOGENESIS_MARKER_SUPPLEMENT = {"CD9", "CD63", "CD81", "CD82"}

# category -> subcategory -> [GO IDs]
CATEGORY_TERMS = {
    "EV Biogenesis": {"—": ["GO:0140112", "GO:0036258"]},
    "Immune-evasion": {"—": ["GO:0050777", "GO:0045916", "GO:0050765", "GO:0002698", "GO:0001915", "GO:0045953"]},
    "ECM Organization": {"—": ["GO:0030198"]},
    "Cardioprotective": {
        "Ca2+ homeostasis": ["GO:0006874", "GO:0070588"],
        "AMPK signaling": ["GO:0031669", "GO:0006112"],
        "Antioxidants/Redox": ["GO:0045454", "GO:0098869", "GO:0006979"],
        "Glycolysis": ["GO:0006096"],
        "UPR": ["GO:0034976", "GO:0006986", "GO:0034975", "GO:0006487"],
        "HSPs": ["GO:0006457", "GO:0034605"],
    },
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


def get_descendants(term_id, parents):
    stack, descendants = [term_id], set()
    while stack:
        cur = stack.pop()
        if cur not in descendants:
            descendants.add(cur)
            stack.extend([c for c, pl in parents.items() if cur in pl])
    return descendants


print("Parsing go-basic.obo (is_a + part_of) ...")
terms, parents = parse_obo(os.path.join(GO_DATA, "go-basic.obo"))

print("Parsing goa_human.gaf.gz (direct annotations) ...")
gene2go = defaultdict(set)
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 10 or cols[8] != "P" or "NOT" in cols[3]:
            continue
        gene2go[cols[2]].add(cols[4])

rows = []
for category, subcats in CATEGORY_TERMS.items():
    for subcat, go_ids in subcats.items():
        desc = set()
        for gid in go_ids:
            desc |= get_descendants(gid, parents)
        genes = {g for g, terms_g in gene2go.items() if terms_g & desc}
        if category == "EV Biogenesis":
            added = EV_BIOGENESIS_MARKER_SUPPLEMENT - genes
            genes |= EV_BIOGENESIS_MARKER_SUPPLEMENT
            print(f"{category} / {subcat}: {len(genes)} genes (BP-derived + {len(added)} marker supplement: {sorted(added)})")
        else:
            print(f"{category} / {subcat}: {len(genes)} genes")
        for g in sorted(genes):
            rows.append({"Gene": g, "Category": category, "Subcategory": subcat})

# Adhesion/Docking/Uptake: reuse the already-built curated set
adhesion_genes = pd.read_csv(ADHESION_SET)["Gene"]
for g in adhesion_genes:
    rows.append({"Gene": g, "Category": "Adhesion / Docking / Uptake", "Subcategory": "—"})
print(f"Adhesion / Docking / Uptake: {len(adhesion_genes)} genes (reused from existing build)")

out_df = pd.DataFrame(rows)
out_df.to_csv(os.path.join(OUT, "category_genesets.csv"), index=False)
print(f"\nSaved: category_genesets.csv ({len(out_df)} gene-category rows, "
      f"{out_df['Gene'].nunique()} unique genes)")
