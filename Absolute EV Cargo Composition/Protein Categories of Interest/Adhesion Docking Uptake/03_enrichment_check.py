"""
Adhesion / Docking / Uptake - Enrichment Re-check (2026-09-25)
==============================================================================
Re-verifies the targeted-enrichment result from 02_targeted_enrichment_test.py
before it is written up, answering two questions:

 1. Reconcile the DEP denominators: the volcano reports 141 (HG down) / 237
    (HL down) DEPs, but 02_targeted_enrichment_test.py reported 128 / 219.
    Where do the missing DEPs go? (Step-by-step accounting below.)
 2. Background choice: 02_targeted_enrichment_test.py tests against ALL human
    genes carrying a BP annotation in GOA (genome-wide), NOT the detected EV
    proteome. Because the detected proteome is itself enriched for this
    category (detected fraction >> genome fraction), a genome-wide background
    can inflate significance. This script re-runs the same one-sided Fisher's
    exact test with the DETECTED proteome (2,866 filtered proteins) as the
    universe, side by side with the original, plus Benjamini-Hochberg
    correction across the 4 tests for each variant.

Same DEP definitions as the volcano: Student's t-test p < 0.05 and
|log2FC| > 0.85, on EVs_df_filtered.csv.

Input:  curated_geneset.csv, ../../../MS proteomics data/Filtered Datasets/EVs_df_filtered.csv,
        GOA human GAF (for the original genome-wide background only)
Output: enrichment_check.csv, enrichment_check_report.txt (this folder)
"""

import gzip
import os

import pandas as pd
from scipy.stats import fisher_exact

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
DATA = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
HERE = os.path.dirname(os.path.abspath(__file__))
CURATED = os.path.join(HERE, "curated_geneset.csv")
GO_DATA = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\05-2026_Comparative_Analysis\08-GO_enrichment\GO_data"

SIG, FC = 0.05, 0.85

df = pd.read_csv(DATA)
curated = set(pd.read_csv(CURATED)["Gene"])
detected = set(df["Genes"].dropna())

# Original background: all human genes with a BP annotation in GOA
bg_genome = set()
with gzip.open(os.path.join(GO_DATA, "goa_human.gaf.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        cols = line.strip().split("\t")
        if len(cols) < 10 or cols[8] != "P" or "NOT" in cols[3]:
            continue
        bg_genome.add(cols[2])

report = []
report.append(f"Filtered proteome: {len(df)} rows, {len(detected)} unique gene names "
              f"({df['Genes'].isna().sum()} rows with missing gene name)")
report.append(f"Genome-wide BP-annotated background (original): {len(bg_genome)} genes")
report.append(f"Curated Adhesion/Docking/Uptake set: {len(curated)} genes total, "
              f"{len(curated & detected)} detected in filtered proteome, "
              f"{len(curated & bg_genome)} in genome background")

COMPARISONS = {
    "HG": ("Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal"),
    "HL": ("Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"),
}

rows = []
report.append("\n=== DEP denominator reconciliation ===")
for comp, (pcol, fccol) in COMPARISONS.items():
    d = df[["Genes", pcol, fccol]].dropna(subset=[pcol, fccol])
    for direction in ("up", "down"):
        m = (d[pcol] < SIG) & ((d[fccol] > FC) if direction == "up" else (d[fccol] < -FC))
        sel = d[m]
        n_rows = len(sel)
        genes_all = set(sel["Genes"].dropna())
        n_missing_name = int(sel["Genes"].isna().sum())
        in_bg = genes_all & bg_genome
        no_bp = sorted(genes_all - bg_genome)
        label = f"{comp} {direction}"

        # Variant A: original (genome-wide BP-annotated background, query restricted to it)
        n_bg = len(bg_genome)
        n_cur_bg = len(curated & bg_genome)
        k_a = len(in_bg & curated)
        q_a = len(in_bg)
        tbl_a = [[k_a, q_a - k_a], [n_cur_bg - k_a, n_bg - n_cur_bg - (q_a - k_a)]]
        p_a = fisher_exact(tbl_a, alternative="greater")[1]

        # Variant B: detected filtered proteome as the universe (all DEPs kept as query)
        n_u = len(detected)
        n_cur_u = len(curated & detected)
        k_b = len(genes_all & curated)
        q_b = len(genes_all)
        tbl_b = [[k_b, q_b - k_b], [n_cur_u - k_b, n_u - n_cur_u - (q_b - k_b)]]
        p_b = fisher_exact(tbl_b, alternative="greater")[1]

        report.append(f"{label}: {n_rows} DEP rows | {n_missing_name} missing gene name | "
                      f"{len(genes_all)} unique genes | {len(in_bg)} carry a GOA BP annotation "
                      f"(original query) | {len(no_bp)} do not: {', '.join(no_bp)}")
        rows.append({
            "DEP_category": label, "n_DEP_rows_volcano": n_rows,
            "n_in_original_query": q_a, "k_in_category_original": k_a,
            "expected_fraction_genome_bg": round(n_cur_bg / n_bg, 4),
            "p_original_genome_background": p_a,
            "n_query_detected_universe": q_b, "k_in_category_detected": k_b,
            "expected_fraction_detected_bg": round(n_cur_u / n_u, 4),
            "observed_fraction": round(k_b / q_b, 4) if q_b else None,
            "p_detected_proteome_background": p_b,
        })

res = pd.DataFrame(rows)


def bh(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    n = len(pvals)
    adj = [None] * n
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, pvals[i] * n / rank)
        adj[i] = prev
    return adj


res["q_original_BH"] = bh(list(res["p_original_genome_background"]))
res["q_detected_BH"] = bh(list(res["p_detected_proteome_background"]))
res.to_csv(os.path.join(HERE, "enrichment_check.csv"), index=False)

report.append("\n=== Fisher's exact (one-sided, greater), original vs. detected-proteome background ===")
report.append(f"Category share of genome background: {len(curated & bg_genome)}/{len(bg_genome)} = "
              f"{len(curated & bg_genome) / len(bg_genome):.1%}; "
              f"share of detected proteome: {len(curated & detected)}/{len(detected)} = "
              f"{len(curated & detected) / len(detected):.1%}")
for _, r in res.iterrows():
    report.append(
        f"{r['DEP_category']}: original {r['k_in_category_original']}/{r['n_in_original_query']} "
        f"p={r['p_original_genome_background']:.3g} (q={r['q_original_BH']:.3g}) | "
        f"detected-bg {r['k_in_category_detected']}/{r['n_query_detected_universe']} "
        f"({r['observed_fraction']:.1%} observed vs {r['expected_fraction_detected_bg']:.1%} expected) "
        f"p={r['p_detected_proteome_background']:.3g} (q={r['q_detected_BH']:.3g})")

with open(os.path.join(HERE, "enrichment_check_report.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(report))
print("done")
