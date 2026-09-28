"""
Fetch STRING interaction networks for the UP-regulated DEP sets (HG up, HL up)
==============================================================================
The down-regulated networks were exported manually from the STRING website.
For the up-regulated sets the same query is made through the STRING web
service (API, https://string-db.org/api/tsv/network): human (9606), the DEP
gene list only (no extra nodes added), combined-score cutoff 0.4 (required_score
= 400), same as the manual exports. Revital approved using the API (2026-09-26).

As a control, the HG-down list is also fetched and compared with the existing
manual export (node/edge counts and edge overlap) to confirm the API returns the
same network.

Outputs (written next to the existing conditions folders):
  HG upregulated/HG_up_ALL_DEPs_STRING_input.txt, string_interactions_HG_upreg.tsv
  HL upregulated/HL_up_ALL_DEPs_STRING_input.txt, string_interactions_HL_upreg.tsv
  HG downregulated/string_interactions_HG_downreg_API_check.tsv (control only)
The TSVs use the web export's column names (node1, node2, combined_score, ...),
so plot_PPI_networks_pie.py reads them unchanged.
"""

import io
import os
import time

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
SIG, FC = 0.05, 0.85
URL = "https://string-db.org/api/tsv/network"

df = pd.read_csv(DATA)


def deps(pcol, fcol, direction):
    d = df[["Genes", pcol, fcol]].dropna()
    m = (d[pcol] < SIG) & ((d[fcol] > FC) if direction == "up" else (d[fcol] < -FC))
    return sorted(d.loc[m, "Genes"])


LISTS = {
    "HG up": deps("Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal", "up"),
    "HL up": deps("Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal", "up"),
    "HG down": deps("Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal", "down"),
}


def fetch(genes):
    r = requests.post(URL, data={"identifiers": "\r".join(genes), "species": 9606, "required_score": 400,
                                 "caller_identity": "thesis_PPI_networks"}, timeout=120)
    r.raise_for_status()
    t = pd.read_csv(io.StringIO(r.text), sep="\t")
    # STRING returns each pair once per direction-free edge; convert to the web-export column names
    out = pd.DataFrame({
        "#node1": t["preferredName_A"], "node2": t["preferredName_B"],
        "node1_string_id": t["stringId_A"], "node2_string_id": t["stringId_B"],
        "neighborhood_on_chromosome": t["nscore"], "gene_fusion": t["fscore"],
        "phylogenetic_cooccurrence": t["pscore"], "homology": 0, "coexpression": t["ascore"],
        "experimentally_determined_interaction": t["escore"], "database_annotated": t["dscore"],
        "automated_textmining": t["tscore"], "combined_score": t["score"],
    })
    return out.drop_duplicates(subset=["#node1", "node2"])


def edge_set(tsv):
    return {frozenset((a, b)) for a, b in zip(tsv.iloc[:, 0], tsv["node2"])}


for label, genes in LISTS.items():
    print(f"{label}: {len(genes)} genes sent")
    net = fetch(genes)
    nodes = set(net["#node1"]) | set(net["node2"])
    print(f"  returned: {len(net)} edges, {len(nodes)} nodes with at least one edge")
    if label == "HG down":
        old = pd.read_csv(os.path.join(HERE, "HG downregulated", "string_interactions_HG_downreg.tsv"), sep="\t")
        old_e, new_e = edge_set(old), edge_set(net)
        print(f"  control vs manual export: manual {len(old_e)} edges, API {len(new_e)} edges, "
              f"shared {len(old_e & new_e)}")
        net.to_csv(os.path.join(HERE, "HG downregulated", "string_interactions_HG_downreg_API_check.tsv"),
                   sep="\t", index=False)
    else:
        cond, _ = label.split()
        folder = os.path.join(HERE, f"{cond} upregulated")
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, f"{cond}_up_ALL_DEPs_STRING_input.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(genes))
        net.to_csv(os.path.join(folder, f"string_interactions_{cond}_upreg.tsv"), sep="\t", index=False)
        print(f"  saved to {folder}")
    time.sleep(2)
