import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# 1. Store your file path in its own variable
file_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\07-SSMD\Barrier_and_Angiogensis_SSMD_results.csv"

# 2. Use the file_path variable to load the DataFrame
df = pd.read_csv(file_path)

# 3. Use the file_path variable AGAIN to get the name
file_name = Path(file_path).stem

# --- Define groups in the order: HG, HL, Combo ---
# Using a list of tuples to preserve order
groups = [
    ("SSMD_HG_vs_Normal", (0.99, 0.57, 0)),        # Orange - HG first
    ("SSMD_HL_vs_Normal", (0.99, 0.14, 0)),          # Red - HL second
    ("SSMD_Combo_vs_Normal", (0.57, 0.125, 0.57)),   # Purple - Combo third
]

# Optional: sort by Combo values
df = df.sort_values("SSMD_Combo_vs_Normal")

n_assays = len(df)
bar_height = 0.25
y_positions = np.arange(n_assays)

# --- FIX 1: Removed duplicate fig, ax line ---
fig, ax = plt.subplots(figsize=(10, max(6, n_assays * 0.4)))

# --- FIX 2: Initialize lists for legend handles and labels ---
handles = []
labels = []

# Plot each SSMD group in the correct order
for i, (col, color) in enumerate(groups):
    label_text = col.replace("SSMD_", "").replace("_vs_Normal", "")
    
    # Store the output of barh (the "handle") in the 'bar' variable
    bar = ax.barh(
        y_positions + (1-i) * bar_height,
        df[col],
        height=bar_height,
        color=color,
        label=label_text
    )
    
    # --- FIX 3: Append the handle and label AFTER the bar is created ---
    handles.append(bar)
    labels.append(label_text)

ax.set_title(f"SSMD Scores for Barrier and Angiogenic function", fontsize=16)
ax.set_yticks(y_positions)
ax.set_yticklabels(df["Assay"])
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("SSMD", fontsize=12)
ax.set_xticks(np.arange(-12, 13, 2))

# --- FIX 4: Pass the ordered handles and labels to the legend ---
ax.legend(handles, labels, title="Groups", loc='best')

# Add reference lines for SSMD thresholds
for threshold in [2]:
    ax.axvline(threshold, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(-threshold, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)

# Add grid for easier reading
# ax.grid(axis='x', alpha=0.3, linestyle=':')

plt.tight_layout()
plt.show()