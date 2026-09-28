"""
Filter single-peptide-detected proteins — Cells & EVs
======================================================
Reads the two Perseus-derived CSVs (Cells, EVs) and removes proteins
identified by only 1 peptide ("reliable" = False), keeping the dataset
consistent with the manually-filtered files already used for 4.2.1/4.2.2
(HUVEC-EVs_df.csv / PCA_new_run.csv).

No detection-threshold filtering applied (checked 2026-08-13: 0/3,274 EV
proteins fail it — Perseus already excludes never-detected proteins).

Input:
  ../Cells_df-report 90777-88new_runs-perseus.csv
  ../EVs_df-report 91540-8 without HG3-perseus.csv

Output (../Filtered Datasets/):
  Cells_df_filtered.csv
  EVs_df_filtered.csv
"""

import os
import pandas as pd

BASE = r"C:\Users\רויטל\Desktop\ISF\Thesis\Results\MS proteomics data"
OUT_DIR = os.path.join(BASE, "Filtered Datasets")
os.makedirs(OUT_DIR, exist_ok=True)

CELLS_IN = os.path.join(BASE, "Cells_df-report 90777-88new_runs-perseus.csv")
EVS_IN   = os.path.join(BASE, "EVs_df-report 91540-8 without HG3-perseus.csv")

# -- Cells --------------------------------------------------------------
cells = pd.read_csv(CELLS_IN)
cells["N_detected_peptides"] = pd.to_numeric(cells["number of peptides"], errors="coerce")
cells["reliable"] = cells["N_detected_peptides"] >= 2

n_total = len(cells)
n_removed = (~cells["reliable"]).sum()
cells_filtered = cells[cells["reliable"]].copy()
cells_filtered.to_csv(os.path.join(OUT_DIR, "Cells_df_filtered.csv"), index=False)

print(f"Cells: {n_total} total | {n_removed} removed (1-peptide) | {len(cells_filtered)} kept")

# -- EVs ------------------------------------------------------------------
evs = pd.read_csv(EVS_IN)
evs["N_detected_peptides"] = pd.to_numeric(evs["peptides"], errors="coerce")
evs["reliable"] = evs["N_detected_peptides"] >= 2

n_total_e = len(evs)
n_removed_e = (~evs["reliable"]).sum()
evs_filtered = evs[evs["reliable"]].copy()
evs_filtered.to_csv(os.path.join(OUT_DIR, "EVs_df_filtered.csv"), index=False)

print(f"EVs  : {n_total_e} total | {n_removed_e} removed (1-peptide) | {len(evs_filtered)} kept")
