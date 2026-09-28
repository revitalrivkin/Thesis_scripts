"""
Targeted Enrichment of Adhesion, Docking & Uptake Proteins Among DEPs
======================================================================
Step 1: Build the curated candidate gene set.

Adapted from Ophir's original descendant-gathering method
(Data analysis projects/Ophir's original script -Gather_Descendants.py):
parses go-basic.obo directly (is_a AND part_of relationships) to get every
descendant of each curated top-level GO BP term, then collects every gene
DIRECTLY annotated (GOA, no propagation needed since we already expanded to
the full descendant subtree) to any of those descendant terms.

This supersedes the earlier exploratory version, which propagated each
gene's own annotations UP to parents using goatools' default GODag (is_a
only, no part_of) — less complete than Ophir's top-down approach.

Curated top-level terms (chosen to match the PI's own network claim:
"downregulation of adhesion and docking proteins, as well as uptake
clusters"):
  GO:0007155  cell adhesion
  GO:0006897  endocytosis
  GO:0006898  receptor-mediated endocytosis
  GO:0016192  vesicle-mediated transport
  GO:0006906  vesicle fusion
  GO:0006904  vesicle docking involved in exocytosis
  GO:0140029  exocytic process
  GO:0032940  secretion by cell

Output: curated_geneset.csv (gene, matched top-level term(s))
        curated_geneset_term_sizes.csv (per-term gene counts, for reporting)
"""

import gzip
import os
from collections import defaultdict

import pandas as pd

GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"
OUT = os.path.dirname(os.path.abspath(__file__))

CURATED_TERMS = {
    "GO:0007155": "cell adhesion",
    "GO:0006897": "endocytosis",
    "GO:0006898": "receptor-mediated endocytosis",
    "GO:0016192": "vesicle-mediated transport",
    "GO:0006906": "vesicle fusion",
    "GO:0006904": "vesicle docking involved in exocytosis",
    "GO:0140029": "exocytic process",
    "GO:0032940": "secretion by cell",
}


# -- Ophir's OBO parser (adapted: is_a + part_of relationships) ----------------
def parse_obo(file_path):
    """Parse OBO file and build parent-child relationships (is_a + part_of)."""
    parents = defaultdict(list)
    terms = {}
    current_term = None

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if line in ["[Term]", "[Typedef]"]:
                if current_term:
                    terms[current_term["id"]] = current_term
                current_term = {}
            elif line.startswith("id: "):
                current_term["id"] = line[4:]
            elif line.startswith("is_a: "):
                parent_id = line.split()[1]
                parents[current_term["id"]].append(parent_id)
            elif line.startswith("relationship: part_of "):
                parent_id = line.split()[2]
                parents[current_term["id"]].append(parent_id)
    if current_term:
        terms[current_term["id"]] = current_term

    return terms, parents


def get_descendants(term_id, parents):
    """Recursively find all descendants of a given term (includes the term itself)."""
    stack = [term_id]
    descendants = set()
    while stack:
        current = stack.pop()
        if current not in descendants:
            descendants.add(current)
            stack.extend([child for child, parent_list in parents.items() if current in parent_list])
    return descendants


print("Parsing go-basic.obo (is_a + part_of) ...")
terms, parents = parse_obo(os.path.join(GO_DATA, "go-basic.obo"))

# -- Descendant set per curated top-level term ----------------------------------
term_descendants = {}
for gid, name in CURATED_TERMS.items():
    desc = get_descendants(gid, parents)
    term_descendants[gid] = desc
    print(f"  {gid} {name}: {len(desc)} descendant terms (incl. self)")

all_curated_term_ids = set().union(*term_descendants.values())
print(f"Union of all curated descendant terms: {len(all_curated_term_ids)}")

# -- Direct gene -> GO BP annotations (no propagation; we already expanded terms) --
print("Parsing goa_human.gaf.gz (direct annotations only) ...")
gene2go_direct = defaultdict(set)
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 10 or cols[8] != "P" or "NOT" in cols[3]:
            continue
        gene2go_direct[cols[2]].add(cols[4])

# -- Which genes fall under each curated top-level term (and the union) --------
term_gene_sets = {}
for gid, name in CURATED_TERMS.items():
    desc = term_descendants[gid]
    genes = {g for g, go_ids in gene2go_direct.items() if go_ids & desc}
    term_gene_sets[gid] = genes
    print(f"  {gid} {name}: {len(genes)} genes")

curated_genes = set().union(*term_gene_sets.values())
print(f"\nCurated gene set total (union across all 8 terms): {len(curated_genes)} genes")

# -- Save outputs ----------------------------------------------------------------
gene_rows = []
for gene in sorted(curated_genes):
    matched_terms = [CURATED_TERMS[gid] for gid, genes in term_gene_sets.items() if gene in genes]
    gene_rows.append({"Gene": gene, "Matched_terms": "; ".join(matched_terms)})
pd.DataFrame(gene_rows).to_csv(os.path.join(OUT, "curated_geneset.csv"), index=False)

size_rows = [{"GO_ID": gid, "Term": name, "N_genes": len(term_gene_sets[gid])}
             for gid, name in CURATED_TERMS.items()]
size_rows.append({"GO_ID": "UNION", "Term": "All curated terms (union)", "N_genes": len(curated_genes)})
pd.DataFrame(size_rows).to_csv(os.path.join(OUT, "curated_geneset_term_sizes.csv"), index=False)

print("\nSaved: curated_geneset.csv")
print("Saved: curated_geneset_term_sizes.csv")
