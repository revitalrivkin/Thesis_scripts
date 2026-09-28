import pandas as pd
import numpy as np
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage
import scipy.cluster.hierarchy as sch
import os


def merge_gaf_with_df(gaf_path, df_path):
    # Load the DataFrame
    df = pd.read_csv(df_path)

    # Initialize lists to store GAF data
    gaf_data = []

    # Read the GAF file
    with gzip.open(gaf_path, 'rt') as f:
        for line in f:
            if line.startswith('!'):
                continue  # Skip comment lines
            parts = line.strip().split('\t')
            if len(parts) > 2:
                uniprot_id = parts[1]
                gene_name = parts[2]
                gaf_data.append((uniprot_id, gene_name, parts))

    # Convert GAF data to DataFrame
    gaf_df = pd.DataFrame(gaf_data, columns=['UniProtID', 'Genes', 'GAF_Data'])

    # Merge with the original DataFrame on 'UniProtID' and 'Gene'
    merged_df = pd.merge(df, gaf_df, how='inner', on=['UniProtID', 'Genes'])

    # Print column names and the first row of the merged DataFrame
    # print("Column Names:", merged_df.columns.tolist())
    # if not merged_df.empty:
    #     print("First Row:\n", merged_df.iloc[1])
    return merged_df

def get_go_term_descendants(go_terms_dict):
    """
    Get all descendants for a dictionary of GO terms.
    
    Parameters:
    go_terms_dict (dict): Dictionary with labels as keys and GO IDs as values
    
    Returns:
    dict: Enhanced dictionary with descendants for each GO term
    """
    from collections import defaultdict
    
    def parse_obo(file_path):
        """Parse OBO file and build parent-child relationships."""
        parents = defaultdict(list)
        terms = {}
        current_term = None
            
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line in ["[Term]", "[Typedef]"]:
                    if current_term:
                        terms[current_term['id']] = current_term
                    current_term = {}
                elif line.startswith("id: "):
                    current_term = {'id': line[4:]}
                elif line.startswith("is_a: ") and current_term:
                    parent_id = line.split()[1]
                    parents[current_term['id']].append(parent_id)
                elif line.startswith("relationship: part_of ") and current_term:
                    parent_id = line.split()[2]
                    parents[current_term['id']].append(parent_id)
                elif line.startswith("relationship: is_a ") and current_term:
                    parent_id = line.split()[2]
                    parents[current_term['id']].append(parent_id)
                    
        if current_term:  # Capture the last term
            terms[current_term['id']] = current_term

        return terms, parents

    def get_descendants(term_id, parents):
        """Recursively find all descendants of a given term."""
        stack = [term_id]
        descendants = set()

        while stack:
            current = stack.pop()
            if current not in descendants:
                descendants.add(current)
                # Find all children that have this term as a parent
                for child, parent_list in parents.items():
                    if current in parent_list:
                        stack.append(child)

        return descendants

    # Input file path - update this to your file path
    obo_file = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\GO\go-basic.obo"
    
    # Process the OBO file
    terms, parents = parse_obo(obo_file)
    
    # Create a new dictionary with descendants for each term
    result_dict = {}
    for label, go_id in go_terms_dict.items():
        descendants = get_descendants(go_id, parents)
        result_dict[label] = {
            "GO_ID": go_id,
            "Descendants": list(descendants)
        }
    
    return result_dict

# First, modify your filter_by_go_terms function to handle the new structure of d_with_descendants
def filter_by_go_terms(df, go_terms_dict):
    """
    Filters the DataFrame to include only rows that contain at least one of the specified GO terms.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame.
    go_terms_dict (dict): Dictionary with GO terms and their descendants.
    
    Returns:
    pd.DataFrame: The filtered DataFrame.
    """
    # Extract all GO terms from the dictionary
    all_go_terms = []
    for label, term_data in go_terms_dict.items():
        all_go_terms.extend(term_data["Descendants"])
    
    # Filter rows that contain at least one of the specified GO terms
    filtered_df = df[df['GAF_Data'].apply(lambda x: any(go_term in str(x) for go_term in all_go_terms))]
    
    # Check if the filtered DataFrame is empty
    if filtered_df.empty:
        print("WARNING: Filtering resulted in an empty DataFrame. No genes match the GO terms.")
        # Return the original df or a subset to prevent errors
        return df.head(10)  # Return some data to avoid the error
    
    filtered_df = filtered_df.drop_duplicates(subset=['Genes'])
    return filtered_df

def create_gene_clustered_heatmap(df, gene_col, normal_label, normal_cols, hg_label, hg_cols, hl_label, hl_cols, combo_label, combo_cols):
    # Calculate the average for each group
    filtered_df = df.copy()
    filtered_df['Normal_avg'] = filtered_df[normal_cols].mean(axis=1)
    filtered_df['HG_avg'] = filtered_df[hg_cols].mean(axis=1)
    filtered_df['HL_avg'] = filtered_df[hl_cols].mean(axis=1)
    filtered_df['Combo_avg'] = filtered_df[combo_cols].mean(axis=1)

    # Prepare the DataFrame for clustering
    df_cluster = filtered_df.set_index(gene_col)[['Normal_avg', 'HG_avg', 'HL_avg', 'Combo_avg']]
    
    # Calculate the z-score for each value in the DataFrame
    df_zscore = df_cluster.apply(zscore, axis=1)
    
    # Handle NaN values by filling them with 0
    df_zscore = df_zscore.fillna(0)

    num_genes = df_zscore.shape[0]
    
    # Print shape information for debugging
    print(f"Shape of z-score matrix: {df_zscore.shape}")
    print(f"Number of genes: {num_genes}")
    
   
    clustered_map = sns.clustermap(
        df_zscore, 
        col_cluster=False,    # Don't cluster columns
        row_cluster=True,     # Cluster rows (genes)
        cmap='coolwarm',
        figsize=(9, 11),
        xticklabels=[normal_label, hg_label, hl_label, combo_label],
        yticklabels=True,
        # cbar=True,  # Disable the colorbar
        # cbar_kws={
        #     'label': 'Z-Score',
        #     'orientation': 'vertical',  
        # },
        dendrogram_ratio=(0.1, 0),        # Adjust dendrogram size (row, col)
        colors_ratio=(0.03, 0.05)         # Adjust colors ratio
    )
    
    # Move colorbar 
    # clustered_map.cax.set_position([0.92, 0.2, 0.03, 0.5])  # [x, y, width, height]
    
    # Adjust the font size for gene labels
    plt.setp(clustered_map.ax_heatmap.get_yticklabels(), fontsize=10, rotation=0)
    # plt.setp(clustered_map.ax_heatmap.get_yticklabels(), fontsize = max(num_genes/2 , 6), rotation=0)

    clustered_map.ax_cbar.remove() #removed color bar text
    
    # # Adjust dendrogram
    # dendro_box = clustered_map.ax_row_dendrogram.get_position()
    # clustered_map.ax_row_dendrogram.set_position([dendro_box.x0, dendro_box.y0, 0.1, dendro_box.height])
    
    # # Adjust heatmap position to accommodate the adjusted dendrogram
    # heatmap_box = clustered_map.ax_heatmap.get_position()
    # clustered_map.ax_heatmap.set_position([heatmap_box.x0 - 0.02, heatmap_box.y0+0.02, heatmap_box.width - 0.05, heatmap_box.height-0.05])
    
    # clustered_map.fig.subplots_adjust(top=0.96)  # Make room for title

    clustered_map.fig.subplots_adjust(top=0.92)
    # Add title
    # clustered_map.fig.suptitle('Heatmap of DEPs relating to '+str(d), fontsize=16)
    clustered_map.fig.suptitle(('Heatmap of DEPs Relating to Glycolysis'), fontsize=16)

    # Remove any y-axis title
    clustered_map.ax_heatmap.set_ylabel('')
    
    # Rotate x-axis labels for better readability
    plt.setp(clustered_map.ax_heatmap.get_xticklabels(), rotation=0)
    
    # Show the plot
    # plt.tight_layout(rect=[0, 0.05, 0.8, 0.7])  # Adjust layout to make room for colorbar at bottom
    plt.show()

    # print list of DEPs
    for label in clustered_map.ax_heatmap.get_yticklabels():
     print(label.get_text())
    
    return clustered_map


def print_cluster_tables(df, gene_col,  normal_label, normal_cols, hg_label, hg_cols, hl_label, hl_cols, combo_label, combo_cols, num_clusters=4):#desc_col,
    # Compute the averages for each condition
    filtered_df = df.copy()
    filtered_df['Normal_avg'] = filtered_df[normal_cols].mean(axis=1)
    filtered_df['HG_avg'] = filtered_df[hg_cols].mean(axis=1)
    filtered_df['HL_avg'] = filtered_df[hl_cols].mean(axis=1)
    filtered_df['Combo_avg'] = filtered_df[combo_cols].mean(axis=1)
    
    # Prepare for clustering
    df_cluster = filtered_df.set_index(gene_col)[['Normal_avg', 'HG_avg', 'HL_avg', 'Combo_avg']]
    df_zscore = df_cluster.apply(zscore, axis=1).fillna(0)  # Normalize and fill NaN values
    
    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(df_zscore, method='ward')
    cluster_labels = sch.fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    filtered_df['Cluster'] = cluster_labels
    
    # Print tables for each cluster
    for cluster in range(1, num_clusters + 1):
        cluster_data = filtered_df[filtered_df['Cluster'] == cluster]
        selected_columns = [gene_col,  'Normal_avg', 'HG_avg', 'HL_avg', 'Combo_avg'] #desc_col,
        table = cluster_data[selected_columns]
        
        print(f"\nCluster {cluster}: ({len(table)} proteins)")
        print(table.to_string(index=False))
        print("=" * 80)
    
    return filtered_df  # Returning the dataframe with cluster labels for further analysis

def print_cluster_tables2(df, gene_col, desc_col, normal_label, normal_cols, hg_label, hg_cols, hl_label, hl_cols, combo_label, combo_cols, num_clusters=4, output_dir="clusters_output"):
    import scipy.cluster.hierarchy as sch
    from scipy.stats import zscore

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Compute the averages for each condition
    filtered_df = df.copy()
    filtered_df['Normal_avg'] = filtered_df[normal_cols].mean(axis=1)
    filtered_df['HG_avg'] = filtered_df[hg_cols].mean(axis=1)
    filtered_df['HL_avg'] = filtered_df[hl_cols].mean(axis=1)
    filtered_df['Combo_avg'] = filtered_df[combo_cols].mean(axis=1)
    
    # Prepare for clustering
    df_cluster = filtered_df.set_index(gene_col)[['Normal_avg', 'HG_avg', 'HL_avg', 'Combo_avg']]
    df_zscore = df_cluster.apply(zscore, axis=1).fillna(0)
    
    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(df_zscore, method='ward')
    cluster_labels = sch.fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    filtered_df['Cluster'] = cluster_labels
    
    # Save and print tables for each cluster
    for cluster in range(1, num_clusters + 1):
        cluster_data = filtered_df[filtered_df['Cluster'] == cluster]
        selected_columns = [gene_col, desc_col, 'Normal_avg', 'HG_avg', 'HL_avg', 'Combo_avg'] 
        table = cluster_data[selected_columns]

        print(f"\nCluster {cluster}: ({len(table)} proteins)")
        print(table.to_string(index=False))
        print("=" * 80)
        
        # Save to CSV
        csv_path = os.path.join(output_dir, f"cluster_{cluster}.csv")
        table.to_csv(csv_path, index=False)

    # Optionally save the entire dataframe with cluster labels
    full_csv_path = os.path.join(output_dir, "all_clusters.csv")
    filtered_df.to_csv(full_csv_path, index=False)
    
    return filtered_df




if __name__ == "__main__":
    #merge data with GAF
    gaf_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\GO\goa_human.gaf.gz"
    df_path = r"C:\Users\רויטל\Desktop\ISF\Data analysis projects\new_runs cell lysates 06-2025\05-heatmaps\filtered_significant_results.csv"
    merged_df=merge_gaf_with_df(gaf_path, df_path)

    #get related GO terms
    d = {
          'glycolytic process': 'GO:0006096',
          'glucose catabolic process': 'GO:0006007',
          'glucose metabolic process': 'GO:0006006'



        }
    d_with_descendants = get_go_term_descendants(d)

    #filter df by specified GO terms
    filtered_df = filter_by_go_terms(merged_df,d_with_descendants )
    
    # print_cluster_tables(filtered_df, 'Genes','Normal', ['Normal-1','Normal-2','Normal-3'], 'HG', ['HG-1','HG-2','HG-3'], 'HL', ['HL-1', 'HL-2','HL-3'], 'Combo', ['Combo-1','Combo-2','Combo-3'])#,'Description'
    
    create_gene_clustered_heatmap(filtered_df, 'Genes', 'Normal', ['Normal-1','Normal-2','Normal-3'], 'HG', ['HG-1','HG-2','HG-3'], 'HL', ['HL-1', 'HL-2','HL-3'], 'Combo', ['Combo-1','Combo-2','Combo-3'])

    
    