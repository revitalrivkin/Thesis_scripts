"""
HG-upregulated PPI network: which proteins are splicing factors / RNA-binding / ribosomal?
==============================================================================
Classes come from GO annotations (goa_human.gaf.gz + go-basic.obo, is_a + part_of, all
descendant terms, no NOT annotations), not from gene-name patterns:
  Splicing         = BP  GO:0008380  RNA splicing (and descendants)
  RNA-binding      = MF  GO:0003723  RNA binding (and descendants)
  Ribosomal        = MF  GO:0003735  structural constituent of ribosome (and descendants)
Nodes = proteins with at least one STRING interaction in string_interactions_HG_upreg.tsv.

Output: HG_up_network_protein_classes.csv (Gene, degree, categories, flags)
"""

import gzip
import os
from collections import defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"

import sys
sys.path.insert(0, os.path.join(HERE, ".."))
from plot_PPI_networks_pie import node_categories  # noqa: E402

CLASSES = {"Splicing": ("P", "GO:0008380"), "RNA-binding": ("F", "GO:0003723"),
           "Ribosomal": ("F", "GO:0003735")}

children = defaultdict(set)
cur = None
with open(os.path.join(GO_DATA, "go-basic.obo"), encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if line in ("[Term]", "[Typedef]"):
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif cur and line.startswith("is_a: "):
            children[line.split()[1]].add(cur)
        elif cur and line.startswith("relationship: part_of "):
            children[line.split()[2]].add(cur)


def descendants(t):
    out, stack = set(), [t]
    while stack:
        x = stack.pop()
        if x not in out:
            out.add(x)
            stack.extend(children.get(x, ()))
    return out


desc = {k: descendants(v[1]) for k, v in CLASSES.items()}
direct = {"P": defaultdict(set), "F": defaultdict(set)}
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10 or "NOT" in c[3] or c[8] not in direct:
            continue
        direct[c[8]][c[2]].add(c[4])

t = pd.read_csv(os.path.join(HERE, "string_interactions_HG_upreg.tsv"), sep="\t")
deg = pd.concat([t.iloc[:, 0], t["node2"]]).value_counts()
rows = []
for g, d in deg.items():
    row = {"Gene": g, "degree": int(d), "categories": "; ".join(node_categories(g))}
    for k, (aspect, _) in CLASSES.items():
        row[k] = bool(direct[aspect].get(g, set()) & desc[k])
    rows.append(row)
out = pd.DataFrame(rows).sort_values("degree", ascending=False)
out.to_csv(os.path.join(HERE, "HG_up_network_protein_classes.csv"), index=False)

n = len(out)
print(f"network proteins: {n}")
for k in CLASSES:
    print(f"{k}: {int(out[k].sum())}")
print("\nSplicing:", ", ".join(out.loc[out["Splicing"], "Gene"]))
rb = out[out["RNA-binding"] & ~out["Splicing"] & ~out["Ribosomal"]]
print(f"\nRNA-binding (not splicing, not ribosomal): {len(rb)}\n", ", ".join(rb["Gene"]))
print("\nRibosomal:", ", ".join(out.loc[out["Ribosomal"], "Gene"]))
print("\nRNA-binding in total (incl. splicing/ribosomal):", int(out["RNA-binding"].sum()))
print("None of the three classes:", int((~out[["Splicing", "RNA-binding", "Ribosomal"]].any(axis=1)).sum()))
