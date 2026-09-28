import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def GO_histogram(df, FDR_col, count_col, Term_col='Term', colormap='Spectral_r'):
    """
    Create a horizontal bar chart for GO terms with color based on FDR values
    
    Parameters:
    -----------
    df : str or pandas.DataFrame
        Either a path to a CSV file or a pandas DataFrame containing the data
    FDR_col : str
        Column name containing the FDR/adjusted p-values
    count_col : str
        Column name containing the count values (used for bar lengths)
    Term_col : str, default='Term'
        Column name containing the GO term descriptions
    colormap : str, default='Spectral_r'
        Matplotlib colormap to use for the bars (default is reversed Spectral)
    """
    # Load data if df is a file path
    if isinstance(df, str):
        df = pd.read_csv(df)
    
    # Sort by Count for better visualization
    df_sorted = df.sort_values(count_col, ascending=True)
    
    # Create a custom colormap based on FDR values
    fdr_values = df_sorted[FDR_col].values
    min_fdr = fdr_values.min()
    max_fdr = fdr_values.max()
    
    # Create a colormap using the specified colormap
    cmap = plt.cm.get_cmap(colormap)
    
    # Using -log10(FDR) for better visualization (higher values = more significant)
    log_fdr = -np.log10(fdr_values)
    min_log_fdr = log_fdr.min()
    max_log_fdr = log_fdr.max()
    
    # Normalize -log10(FDR) values to [0, 1] for coloring
    norm_fdr = (log_fdr - min_log_fdr) / (max_log_fdr - min_log_fdr) if max_log_fdr > min_log_fdr else np.zeros_like(log_fdr)
    colors = [cmap(nf) for nf in norm_fdr]
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    bars = ax.barh(df_sorted[Term_col], df_sorted[count_col], color=colors)
    
    # Add a colorbar to show the FDR scale (-log10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min_log_fdr, max_log_fdr))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label='-log10(FDR)')
    
    # Set Arial font for colorbar label
    cbar.set_label('-log10(FDR)', fontname='Arial', fontsize=12)
    
    # Set Arial font for colorbar tick labels
    for label in cbar.ax.get_yticklabels():
        label.set_fontname('Arial')
    
    # Annotate with FDR values
    # for bar, fdr in zip(bars, df_sorted[FDR_col]):
    #     ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
    #             f"FDR={fdr:.1e}", va='center', fontsize=8)
    
    # Add labels and title
    ax.set_xlabel(count_col, fontname='Arial', fontsize=12)
    ax.set_ylabel(Term_col, fontname='Arial', fontsize=12)
    ax.set_title(f"Enriched GO Terms for Metabolism DEPs", fontname='Arial', fontsize=17, fontweight='bold')
    
    # Set Arial font for axis tick labels
    for label in ax.get_xticklabels():
        label.set_fontname('Arial')
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    
    plt.tight_layout()
    plt.show()
    
    return fig  # Return the figure object for further customization if needed


if __name__ == "__main__":
    
    df_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\06-GO histograms\Metabolism_GO.csv"
    
    GO_histogram(
        df=df_path, 
        FDR_col="FDR", 
        count_col='Count',
        Term_col='Term',
        colormap='RdYlBu_r'  # Try different colormaps like 'viridis', 'plasma', 'Spectral_r', 'RdYlBu_r'
    )