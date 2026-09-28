import sys
sys.path.append(r"C:\Users\רויטל\Desktop\ISF\Proteomics")
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage
from matplotlib.lines import Line2D

# Set default font to Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']


def create_volcano_plot(df, group1_name, group2_name, pval_col, log2fc_col, protein_col, significance_threshold=0.05, log2fc_threshold=0.85):
    """
    Generates a volcano plot comparing two groups from a DataFrame.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group1_name: str, name of the first group for labeling.
    - group2_name: str, name of the second group for labeling.
    - pval_col: str, name of the column containing p-values.
    - log2fc_col: str, name of the column containing log2 fold change values.
    - protein_col: str, name of the column containing protein names/identifiers.
    - significance_threshold: float, p-value threshold for significance.
    - log2fc_threshold: float, log2 fold change threshold for biological significance.
    """
    # Calculate -log10 of p-values
    df['-log10(p-value)'] = -np.log10(df[pval_col])

    # Determine significance based on both p-value and log2FC thresholds
    df['Significance'] = 'Not Significant'
    df.loc[(df[pval_col] < significance_threshold) & (df[log2fc_col] > log2fc_threshold), 'Significance'] = 'Upregulated'
    df.loc[(df[pval_col] < significance_threshold) & (df[log2fc_col] < -log2fc_threshold), 'Significance'] = 'Downregulated'

    
    # Create and print lists of upregulated and downregulated proteins
    upregulated_proteins = df[df['Significance'] == 'Upregulated'][protein_col].tolist()
    downregulated_proteins = df[df['Significance'] == 'Downregulated'][protein_col].tolist()
    
    print(f"\nUpregulated proteins ({len(upregulated_proteins)}):")
    for protein in upregulated_proteins:
        print(f"{protein}")
        
    print(f"\nDownregulated proteins ({len(downregulated_proteins)}):")
    for protein in downregulated_proteins:
        print(f"{protein}")
    
    # Return these lists for further use if needed
    result_dict = {
        'upregulated': upregulated_proteins,
        'downregulated': downregulated_proteins
    }

    # Plot
    # plt.figure(figsize=(8, 8))
    
    # Use scatter directly for annotation purposes
    fig, ax = plt.subplots(figsize=(7, 7))
    
    # Define RGB colors for each category
    rgb_colors = {
        'Upregulated': (0.57, 0.125, 0.57),        # Purple
        'Downregulated': (0.57, 0.125, 0.57),        # Purple
        # 'Downregulated': (0.78, 0.55, 0.77),      # LightPurple
        'Not Significant': (0.7, 0.7, 0.7)     # Gray
    }
    # Create scatter plot by significance category
    for category, color in rgb_colors.items():
        subset = df[df['Significance'] == category]
        scatter = ax.scatter(
            subset[log2fc_col], 
            subset['-log10(p-value)'], 
            color=[color],
            alpha=0.7,
            picker=True  # Enable picking for hover functionality
        )

    # Add lines for thresholds
    p_threshold_line = ax.axhline(y=-np.log10(significance_threshold), color='gray', linestyle='--', linewidth=1)
    fc_threshold_lines = [
        ax.axvline(x=log2fc_threshold, color='gray', linestyle='--', linewidth=1),
        ax.axvline(x=-log2fc_threshold, color='gray', linestyle='--', linewidth=1)
    ]

    # Customize legend
    # legend_elements = [
    #     Line2D([0], [0], marker='o', color='w', label='Upregulated', markersize=8, markerfacecolor=rgb_colors['Upregulated']),
    #     Line2D([0], [0], marker='o', color='w', label='Downregulated', markersize=8, markerfacecolor=rgb_colors['Downregulated']),
    #     Line2D([0], [0], color='gray', lw=1, linestyle='--', label=f'P-value threshold = {significance_threshold}'),
    #     Line2D([0], [0], color='gray', lw=1, linestyle='--', label=f'Log2FC threshold = {log2fc_threshold}')
    # ]
    # ax.legend(handles=legend_elements, loc='upper right')
    
    # Add annotation capabilities
    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"),
                        fontfamily='Arial')
    annot.set_visible(False)
    
    # Create a dictionary mapping coordinates to protein names for hover functionality
    point_to_protein = {}
    for _, row in df.iterrows():
        x, y = row[log2fc_col], -np.log10(row[pval_col])
        point_key = (round(x, 4), round(y, 4))  # Round to reduce floating point issues
        point_to_protein[point_key] = row[protein_col]
    
    def update_annot(ind, scatter_obj):
        """Update the annotation with the protein name."""
        pos = scatter_obj.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        # Find the closest point in our dictionary
        x, y = pos
        closest_point = min(point_to_protein.keys(), 
                           key=lambda p: abs(p[0]-x) + abs(p[1]-y))
        text = point_to_protein[closest_point]
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.8)
    
    def hover(event):
        """Handle hover events to show protein name."""
        if event.inaxes == ax:
            for scatter_obj in ax.collections:
                cont, ind = scatter_obj.contains(event)
                if cont:
                    update_annot(ind, scatter_obj)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
            # If we get here, we're not over a point
            if annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

    # Labels and title
    ax.set_xlabel(f'Log2 Fold Change ({group1_name} / {group2_name})', weight='bold',size='13')
    ax.set_ylabel('-Log10(p-value)', weight='bold',size='13')
    ax.set_title(f'{group1_name} vs {group2_name}', x=0.5, y=1.0, size='17', weight='bold')
    
    # Make tick labels bold and one size larger
    ax.tick_params(axis='both', which='major', labelsize=12, width=1.5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_weight('bold')
        label.set_fontsize(12)

    # Connect the hover event to the figure
    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    # Show plot
    plt.tight_layout()
    plt.show()

    # Return the dictionary with upregulated and downregulated protein lists
    return result_dict

if __name__ == "__main__":  
    # Load data using pandas
    df1 = pd.read_csv(r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\03-volcano\new_runs_df.csv")
    
    # Example usage with a column named 'log2FC' that contains pre-calculated log2 fold change values
    # Replace 'log2FC' with the actual name of your log2 fold change column
    create_volcano_plot(
        df=df1, 
        group1_name='Combo', 
        group2_name='Normal', 
        pval_col="Student T-test p-value Combo_Normal", 
        log2fc_col='Student T-test Difference log2 (FC) Combo_Normal',  # Replace with your actual log2FC column name
        protein_col='Genes',  # Replace with your actual protein name column
        significance_threshold=0.05,
        log2fc_threshold=0.85  # Adjust this threshold for your desired level of biological significance
    )
