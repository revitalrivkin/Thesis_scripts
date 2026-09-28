"""
Categories of Interest - GO Terms Reference Table (Excel)
==============================================================================
Reference table for 4.2.2.2: summarizes the GO BP terms used to build each
protein category/subcategory gene set, one row per GO term (long format),
plus a per-category/subcategory summary sheet with GO-set size AND the
number of those genes actually detected in the filtered EV proteome
(MS proteomics data/Filtered Datasets/EVs_df_filtered.csv, n=2,866 proteins)
- the detected count is what actually matters for the thesis table, not the
genome-wide GO-set size alone.

Updated 2026-09-25 (Revital's request, verified against category_genesets.csv
and EVs_df_filtered.csv directly, not just carried over from prior notes):
  - Added N_detected_filtered column (computed fresh from the two source
    files each run, not hard-coded) to the Category Summary sheet.
  - Replaced all em dashes with plain hyphens (house style).
  - Added the 3 NEW Adhesion/Docking/Uptake subcategories adopted from the
    Faculty Retreat deck (Presentations\Faculty retreat\Proteomics\
    plot_adhesion_docking_uptake_clusters_slide.py) - Consistently Up,
    Consistently Down, Discordant. These are fundamentally different from
    every other row in this table: they are NOT GO-derived (no defining GO
    term - membership is empirical, based on direction of log2FC change
    under HG/HL), and their scope is DEPs only (127 total Adhesion/Docking/
    Uptake DEPs), not all detected proteins in the category - "direction of
    change" isn't a meaningful label for a non-significant protein. The GO
    term shown for these 3 rows is the POST-HOC enrichment result used to
    name the cluster after it was found, not a term that defines membership.

See [[12 - Thesis Writing]] for the underlying discussion, coverage
validation against the PI's cardioprotective list, and the subcategories'
full derivation (fcluster tested first, didn't give clean breaks; the
simple up/down/discordant split is what was kept).

IMPORTANT: the "GO Terms" sheet lists TOP-LEVEL SEED terms only (for the 4
GO-derived categories + Cardioprotective's 6 subcategories). Each is expanded
to its full descendant subtree (is_a + part_of, Ophir's method - Data
analysis projects\Ophir's original script -Gather_Descendants.py) before
matching genes - individual child terms are not enumerated here (some
subtrees run 250+ terms deep).

Output: categories_GO_reference_table.xlsx (this folder)
"""

import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results"
FILTERED = os.path.join(BASE, "MS proteomics data", "Filtered Datasets", "EVs_df_filtered.csv")
GENESETS = os.path.join(HERE, "category_genesets.csv")

# -- One row per GO term (long format) -------------------------------------------
TERM_ROWS = [
    ("Adhesion / Docking / Uptake", "-", "GO:0007155", "cell adhesion"),
    ("Adhesion / Docking / Uptake", "-", "GO:0006897", "endocytosis"),
    ("Adhesion / Docking / Uptake", "-", "GO:0006898", "receptor-mediated endocytosis"),
    ("Adhesion / Docking / Uptake", "-", "GO:0016192", "vesicle-mediated transport"),
    ("Adhesion / Docking / Uptake", "-", "GO:0006906", "vesicle fusion"),
    ("Adhesion / Docking / Uptake", "-", "GO:0006904", "vesicle docking involved in exocytosis"),
    ("Adhesion / Docking / Uptake", "-", "GO:0140029", "exocytic process"),
    ("Adhesion / Docking / Uptake", "-", "GO:0032940", "secretion by cell"),

    ("Adhesion / Docking / Uptake", "Consistently Up (DEPs only, empirical)",
     "(not GO-derived)", "Post-hoc GO enrichment theme: Golgi / secretory vesicle transport (FDR=4e-13)"),
    ("Adhesion / Docking / Uptake", "Consistently Down (DEPs only, empirical)",
     "(not GO-derived)", "Post-hoc GO enrichment theme: cell-matrix adhesion + receptor-mediated endocytosis (FDR=3e-12 / 5e-8)"),
    ("Adhesion / Docking / Uptake", "Discordant (DEPs only, empirical)",
     "(not GO-derived)", "Post-hoc GO enrichment theme: ECM / basement membrane organization"),

    ("EV Biogenesis", "-", "GO:0140112", "extracellular vesicle biogenesis"),
    ("EV Biogenesis", "-", "GO:0036258", "multivesicular body assembly"),
    ("EV Biogenesis", "-", "(hard-coded, not a GO term)", "Marker supplement: CD9, CD63, CD81, CD82"),

    ("Immune-evasion", "-", "GO:0050777", "negative regulation of immune response"),
    ("Immune-evasion", "-", "GO:0045916", "negative regulation of complement activation"),
    ("Immune-evasion", "-", "GO:0050765", "negative regulation of phagocytosis"),
    ("Immune-evasion", "-", "GO:0002698", "negative regulation of immune effector process"),
    ("Immune-evasion", "-", "GO:0001915", "negative regulation of T cell-mediated cytotoxicity"),
    ("Immune-evasion", "-", "GO:0045953", "negative regulation of NK cell-mediated cytotoxicity"),

    ("ECM Organization", "-", "GO:0030198", "extracellular matrix organization"),

    ("Cardioprotective", "Ca2+ homeostasis", "GO:0006874", "intracellular calcium ion homeostasis"),
    ("Cardioprotective", "Ca2+ homeostasis", "GO:0070588", "calcium ion transmembrane transport"),

    ("Cardioprotective", "AMPK signaling", "GO:0031669", "cellular response to nutrient levels"),
    ("Cardioprotective", "AMPK signaling", "GO:0006112", "energy reserve metabolic process"),

    ("Cardioprotective", "Antioxidants / Redox", "GO:0045454", "cell redox homeostasis"),
    ("Cardioprotective", "Antioxidants / Redox", "GO:0098869", "cellular oxidant detoxification"),
    ("Cardioprotective", "Antioxidants / Redox", "GO:0006979", "response to oxidative stress"),

    ("Cardioprotective", "Glycolysis", "GO:0006096", "glycolytic process"),

    ("Cardioprotective", "UPR", "GO:0034976", "response to endoplasmic reticulum stress"),
    ("Cardioprotective", "UPR", "GO:0006986", "response to unfolded protein"),
    ("Cardioprotective", "UPR", "GO:0034975", "protein folding in endoplasmic reticulum"),
    ("Cardioprotective", "UPR", "GO:0006487", "protein N-linked glycosylation"),

    ("Cardioprotective", "HSPs", "GO:0006457", "protein folding"),
    ("Cardioprotective", "HSPs", "GO:0034605", "cellular response to heat"),

    ("Cardioprotective", "Preconditioning-associated", "-", "No clean GO term - dropped as a GO-derived subcategory"),
]
terms_df = pd.DataFrame(TERM_ROWS, columns=["Category", "Subcategory", "GO_ID", "GO_Term_Name"])

# -- Detected-in-filtered-data counts, computed fresh (not hard-coded) -----------
evs = pd.read_csv(FILTERED)
detected_genes = set(evs["Genes"].dropna())

gs = pd.read_csv(GENESETS)
detected_by_subcat = {}
for (cat, subcat), sub in gs.groupby(["Category", "Subcategory"]):
    genes = set(sub["Gene"])
    detected_by_subcat[(cat, subcat)] = len(genes & detected_genes)

cardio_genes = set(gs[gs["Category"] == "Cardioprotective"]["Gene"])
cardio_detected_unique = len(cardio_genes & detected_genes)

ev_bio_genes = set(gs[gs["Category"] == "EV Biogenesis"]["Gene"])
markers = {"CD9", "CD63", "CD81", "CD82"}
ev_bio_detected = len(ev_bio_genes & detected_genes)
ev_bio_marker_detected = len(markers & detected_genes)

# -- Adhesion/Docking/Uptake direction-of-change DEP subclusters (empirical) -----
SIG_THRESHOLD, LOG2FC_THRESHOLD = 0.05, 0.85
d = evs[["Genes", "Student's T-test p-value HG_Normal", "Student's T-test Difference HG_Normal",
         "Student's T-test p-value HL_Normal", "Student's T-test Difference HL_Normal"]].copy()
d = d.rename(columns={"Student's T-test Difference HG_Normal": "HG",
                       "Student's T-test Difference HL_Normal": "HL",
                       "Student's T-test p-value HG_Normal": "p_HG",
                       "Student's T-test p-value HL_Normal": "p_HL"})
d = d.dropna(subset=["HG", "HL"]).set_index("Genes")
is_dep_hg = (d["p_HG"] < SIG_THRESHOLD) & (d["HG"].abs() > LOG2FC_THRESHOLD)
is_dep_hl = (d["p_HL"] < SIG_THRESHOLD) & (d["HL"].abs() > LOG2FC_THRESHOLD)
d = d[(is_dep_hg.fillna(False)) | (is_dep_hl.fillna(False))]

adhesion_genes = set(gs[gs["Category"] == "Adhesion / Docking / Uptake"]["Gene"])
d = d[d.index.isin(adhesion_genes)]
n_adhesion_deps = len(d)
n_up = int(((d["HG"] > 0) & (d["HL"] > 0)).sum())
n_down = int(((d["HG"] < 0) & (d["HL"] < 0)).sum())
n_discordant = n_adhesion_deps - n_up - n_down

print(f"Detected genes in filtered EV dataset: {len(detected_genes)}")
print(f"Adhesion/Docking/Uptake DEPs: {n_adhesion_deps} (Up={n_up}, Down={n_down}, Discordant={n_discordant})")

# -- Per-category/subcategory summary --------------------------------------------
SUMMARY_ROWS = [
    ("Adhesion / Docking / Uptake", "-", 2169, detected_by_subcat.get(("Adhesion / Docking / Uptake", "-"), None),
     "Enrichment-validated: significantly enriched among downregulated DEPs "
     "(HG down 34/128 p=4.6e-12, HL down 67/219 p=3.3e-26)"),
    ("Adhesion / Docking / Uptake", "Consistently Up (DEPs only, empirical)", None, n_up,
     f"NOT GO-derived - empirical direction-of-change cluster among Adhesion/Docking/Uptake DEPs (n={n_adhesion_deps} total DEPs, not all-detected). "
     "Largely HL-specific (MICALL1, EHBP1L1, PNN, TMED5, TOR1A, CLASP1 near-zero in HG, strongly induced in HL)."),
    ("Adhesion / Docking / Uptake", "Consistently Down (DEPs only, empirical)", None, n_down,
     "NOT GO-derived - empirical direction-of-change cluster; near word-for-word match to the PI's own prior mechanistic claim."),
    ("Adhesion / Docking / Uptake", "Discordant (DEPs only, empirical)", None, n_discordant,
     "NOT GO-derived - empirical direction-of-change cluster; HG and HL move in opposite directions."),
    ("EV Biogenesis", "-", 59, ev_bio_detected,
     f"55 GO-derived (BP terms) + 4 hard-coded markers (CD9/CD63/CD81/CD82, all {ev_bio_marker_detected}/4 detected). "
     "Blanket CC-term expansion was tested to fix the marker gap but overcorrected badly (55->2,346 genes, "
     "mostly generic abundant proteins); reverted to a precise hard-coded supplement instead. "
     "Has 2 DEP genes: CD82 (HG+HL down), CD63 (HL down, newly recovered) - up from 0."),
    ("Immune-evasion", "-", 242, detected_by_subcat.get(("Immune-evasion", "-"), None), "Not yet validated against DEPs"),
    ("ECM Organization", "-", 236, detected_by_subcat.get(("ECM Organization", "-"), None),
     "Matches existing RA and unbiased-GO findings (ECM-dominated cargo)"),
    ("Cardioprotective", "(all subcategories, unique genes)", None, cardio_detected_unique,
     "Unique detected genes across all 6 subcategories combined (some genes belong to more than one)"),
    ("Cardioprotective", "Ca2+ homeostasis", None, detected_by_subcat.get(("Cardioprotective", "Ca2+ homeostasis"), None),
     "2/4 PI-list genes covered (misses S100A10, SDF4)"),
    ("Cardioprotective", "AMPK signaling", None, detected_by_subcat.get(("Cardioprotective", "AMPK signaling"), None),
     "1/3 PI-list genes covered (misses ACACA, PPP2R2A)"),
    ("Cardioprotective", "Antioxidants / Redox", None, detected_by_subcat.get(("Cardioprotective", "Antioxidants/Redox"), None),
     "12/12 PI-list genes covered"),
    ("Cardioprotective", "Glycolysis", None, detected_by_subcat.get(("Cardioprotective", "Glycolysis"), None),
     "10/11 PI-list genes covered (misses TALDO1 - pentose phosphate pathway, not glycolysis)"),
    ("Cardioprotective", "UPR", None, detected_by_subcat.get(("Cardioprotective", "UPR"), None), "8/8 PI-list genes covered"),
    ("Cardioprotective", "HSPs", None, detected_by_subcat.get(("Cardioprotective", "HSPs"), None), "7/7 PI-list genes covered"),
    ("Cardioprotective", "Preconditioning-associated", None, None, "No clean GO term - dropped as a GO-derived subcategory"),
]
summary_df = pd.DataFrame(SUMMARY_ROWS, columns=["Category", "Subcategory", "N_genes_in_GO_set", "N_detected_filtered", "Validation_notes"])

# -- Write Excel with basic formatting ------------------------------------------
out_path = os.path.join(HERE, "categories_GO_reference_table.xlsx")
with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
    terms_df.to_excel(writer, sheet_name="GO Terms (long format)", index=False)
    summary_df.to_excel(writer, sheet_name="Category Summary", index=False)

    from openpyxl.styles import Font, PatternFill, Alignment
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill(start_color="333333", end_color="333333", fill_type="solid")

    for sheet_name, df in [("GO Terms (long format)", terms_df), ("Category Summary", summary_df)]:
        ws = writer.sheets[sheet_name]
        for cell in ws[1]:
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(vertical="center")
        for col_cells in ws.columns:
            length = max(len(str(c.value)) if c.value is not None else 0 for c in col_cells)
            ws.column_dimensions[col_cells[0].column_letter].width = min(max(length + 2, 12), 60)
        ws.freeze_panes = "A2"

print("Saved:", out_path)
