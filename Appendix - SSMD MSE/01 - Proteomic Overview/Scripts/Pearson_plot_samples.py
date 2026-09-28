import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors


def create_pairwise_correlation_map2(df):
    """
    Creates a pairwise correlation map showing only the diagonal and lower triangular part,
    with a horizontal colorbar at the bottom.

    Args:
        df: The pandas DataFrame containing the data.
    
    Returns:
        None (displays the correlation map directly).
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Calculate the correlation matrix
    corr_matrix = df.corr()
    
    # Create a mask for the upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    
    # Set up the matplotlib figure
    plt.figure(figsize=(10, 8))
    
    # Set Arial font for all text
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']

    # Create a monochromatic colormap (e.g., blue)
    # cmap = plt.cm.Purples  # You can choose Blues, Reds, Greens, Purples, Greys, etc.
    cmap = plt.cm.Greens
    
    # Draw the heatmap with the mask
    cbar_kws = {"orientation": "vertical", "shrink": 0.8, "aspect": 20, "pad": 0.12}
    sns.heatmap(corr_matrix, mask=mask, annot=True, cmap=cmap, 
                square=True, linewidths=.5, cbar_kws=cbar_kws)
    plt.xticks(rotation=45, ha='right',fontweight='bold')
    plt.yticks(fontweight='bold')
    
    plt.title('Pairwise Correlation Map', fontsize=17, fontweight='bold')
    # plt.tight_layout()
    plt.show()
    




if __name__ == "__main__":
    df = pd.read_csv(r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\02-Pearson\pearson_df_n.csv")
    
    # Specify the sample columns to use for correlation
    sample_cols = ['Normal-1','Normal-2','Normal-3', 'HG-1','HG-2','HG-3','HL-1', 'HL-2','HL-3','Combo-1','Combo-2','Combo-3']
    
    # Filter the dataframe to only include these columns
    df_filtered = df[sample_cols]
    
    # Create correlation map with filtered data
    create_pairwise_correlation_map2(df_filtered)