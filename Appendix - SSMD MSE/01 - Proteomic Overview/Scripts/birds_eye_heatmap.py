import sys
sys.path.append(r"C:\Users\רויטל\Desktop\ISF\Proteomics")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import zscore
import warnings

# Try to import fastcluster which is more efficient for large matrices
try:
    import fastcluster
    linkage_function = fastcluster.linkage
    print("Using fastcluster for enhanced performance")
except ImportError:
    from scipy.cluster.hierarchy import linkage as scipy_linkage
    linkage_function = scipy_linkage
    print("Consider installing fastcluster: pip install fastcluster")

# Increase recursion limit to handle large datasets
sys.setrecursionlimit(10000)


def create_heatmap_for_merged_df(df, gene_col, level_cols, sample_size=None, cluster_cols=False):
    """
    Creates a heatmap with hierarchical clustering and a dendrogram for the merged_df DataFrame.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    gene_col (str): The column name for genes.
    level_cols (list): The list of columns to include in clustering.
    sample_size (int, optional): Number of rows to sample if dataset is large. None means use all rows.
    cluster_cols (bool): Whether to cluster columns or keep original order.

    Returns:
    None
    """
    # Make a copy of the DataFrame to avoid modifying the original
    df_copy = df.copy()
    
    # Ensure the level_cols are numeric
    df_copy[level_cols] = df_copy[level_cols].apply(pd.to_numeric, errors='coerce')
    
    # Reorder columns to match the specified order
    df_copy = df_copy[['Genes'] + level_cols]
    
    # Set the gene column as the index
    df_data = df_copy.set_index(gene_col)[level_cols]
    
    # Drop rows with any NaN values to ensure clean clustering
    df_data = df_data.dropna()
    
    # If dataset is large, sample rows to avoid excessive computation
    if sample_size and len(df_data) > sample_size:
        print(f"Dataset is large ({len(df_data)} rows). Sampling {sample_size} rows for clustering.")
        df_data = df_data.sample(sample_size, random_state=42)
    
    try:
        # Suppress specific warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Clustering large matrix with scipy")
            
            # Create a clustermap with row dendrogram and optional column dendrogram
            g = sns.clustermap(
                df_data,
                figsize=(12, 10),
                cmap='coolwarm',
                z_score=0,  # Z-score normalize rows
                col_cluster=cluster_cols,  # Only cluster columns if specified
                row_cluster=True,
                method='average',  # Use average linkage for better clustering
                metric='euclidean',
                xticklabels=True,
                yticklabels=False,  # No gene labels as requested
                dendrogram_ratio=(0.1, 0.2),  # Adjust dendrogram size
                cbar_pos=(0.02, 0.8, 0.05, 0.18),  # Position the colorbar
                tree_kws={'linewidths': 0.5}  # Make dendrograms more visible
            )
            
            # Set title
            plt.suptitle('Gene Expression Heatmap', fontsize=16, y=0.98)
            
            # Format column labels
            plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right', fontsize=10)
            
            # Set colorbar title
            g.cax.set_title('Z-Score', fontsize=10)
            
            # Show the plot
            plt.tight_layout(rect=[0, 0, 1, 0.97])  # Leave space for title
            plt.show()
        
    except Exception as e:
        print(f"Error in heatmap generation: {str(e)}")
        
        # If recursion error or memory error, try with a smaller sample
        if isinstance(e, (RecursionError, MemoryError)) and (sample_size is None or sample_size > 100):
            new_sample = 100 if sample_size is None else min(sample_size // 2, 100)
            print(f"Trying again with a smaller sample size: {new_sample}")
            create_heatmap_for_merged_df(df, gene_col, level_cols, sample_size=new_sample, cluster_cols=cluster_cols)
        else:
            # Last resort: basic heatmap without clustering
            print("Falling back to basic heatmap without clustering...")
            plt.figure(figsize=(10, 8))
            
            # Compute z-scores for the visualization
            df_z = pd.DataFrame(data=np.zeros(df_data.shape), 
                              index=df_data.index, 
                              columns=df_data.columns)
                              
            for idx in df_data.index:
                row = df_data.loc[idx]
                row_std = row.std()
                if row_std != 0:
                    df_z.loc[idx] = (row - row.mean()) / row_std
            
            # Create a simple heatmap
            ax = sns.heatmap(df_z, cmap='coolwarm', xticklabels=True, yticklabels=False)
            
            plt.title('Gene Expression Heatmap (Basic Version)')
            plt.xticks(rotation=45, ha='right', fontsize=10)
            
            plt.tight_layout()
            plt.show()


if __name__ == "__main__":
    # === USER: Set these variables as needed ===
    # csv_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\04-PCA\new_runs_df.csv"  # Path to your CSV file
    csv_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\05-heatmaps\filtered_significant_results.csv"
    gene_col = "Genes"  # Change if your gene column is named differently
    level_cols = ['Normal-1','Normal-2','Normal-3', 'HG-1','HG-2','HG-3','HL-1', 'HL-2','HL-3','Combo-1','Combo-2','Combo-3']  # Change as needed

    try:
        df = pd.read_csv(csv_path)
        # Check for required columns
        missing = [col for col in [gene_col] + level_cols if col not in df.columns]
        if missing:
            raise ValueError(f"Missing columns in CSV: {missing}")
        # Set sample_size if large
        sample_size = 500 if len(df) > 500 else None
        create_heatmap_for_merged_df(
            df,
            gene_col,
            level_cols,
            sample_size=sample_size,
            cluster_cols=False
        )
    except Exception as e:
        print(f"Error: {e}")
