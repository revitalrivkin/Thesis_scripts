import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


def calculate_ssmd_for_row(row, conditions, reference_condition='Normal'):
    """
    Calculate SSMD values for a single row of data.
    
    Args:
        row: pandas Series with columns like <condition>_mean, <condition>_SD, <condition>_N
        conditions: list of condition names
        reference_condition: name of reference condition (default 'Normal')
    
    Returns:
        dict: Dictionary with SSMD values for each condition vs reference
    """
    means = {}
    sds = {}
    ns = {}
    
    for cond in conditions:
        means[cond] = row[f"{cond}_mean"]
        sds[cond] = row[f"{cond}_SD"]
        ns[cond] = row[f"{cond}_N"]
    
    # Step 1: Normalize means
    mu_min = min(means.values())
    mu_max = max(means.values())
    range_val = mu_max - mu_min if mu_max != mu_min else 1
    
    means_norm = {c: (means[c] - mu_min) / range_val for c in conditions}
    
    # Step 2: Normalize variances
    variances = {c: sds[c]**2 for c in conditions}
    variances_norm = {c: variances[c] / (range_val**2) for c in conditions}
    
    # Step 3: Calculate SSMD vs reference
    mu_ref_norm = means_norm[reference_condition]
    var_ref_norm = variances_norm[reference_condition]
    
    ssmd_results = {}
    
    for cond in conditions:
        if cond == reference_condition:
            continue
        mu_cond_norm = means_norm[cond]
        var_cond_norm = variances_norm[cond]
        
        denom = np.sqrt(var_ref_norm + var_cond_norm)
        ssmd = 0 if denom == 0 else (mu_cond_norm - mu_ref_norm) / denom
        
        ssmd_results[f'SSMD_{cond}_vs_{reference_condition}'] = ssmd
    
    return ssmd_results


def process_all_csvs(processed_data_dir, reference_condition='Normal', file_order=None):
    """
    Process all CSV files in processed_data directory and calculate SSMD values.
    
    Args:
        processed_data_dir: Path to directory containing processed CSV files
        reference_condition: Reference condition name (default 'Normal')
        file_order: Optional list of filenames (without extension) in desired order.
                   If None, files are processed alphabetically.
    
    Returns:
        pandas.DataFrame: Combined results with Assay column and SSMD columns
    """
    processed_data_path = Path(processed_data_dir)
    
    if not processed_data_path.exists():
        raise ValueError(f"Directory not found: {processed_data_dir}")
    
    all_results = []
    
    # Get all CSV files
    all_csv_files = list(processed_data_path.glob('*.csv'))
    
    if not all_csv_files:
        print(f"Warning: No CSV files found in {processed_data_dir}")
        return pd.DataFrame()
    
    # If file_order is provided, sort files according to it
    if file_order is not None:
        # Create a mapping of filename (stem) to full path
        file_dict = {f.stem: f for f in all_csv_files}
        csv_files = []
        # Track which files we've added
        added_files = set()
        # Add files in the specified order
        for filename in file_order:
            if filename in file_dict:
                csv_files.append(file_dict[filename])
                added_files.add(filename)
            else:
                print(f"Warning: File '{filename}' from order list not found in directory")
        # Add any remaining files not in the order list
        for f in all_csv_files:
            if f.stem not in added_files:
                csv_files.append(f)
    else:
        # Default: alphabetical order
        csv_files = sorted(all_csv_files)
    
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            assay_name = csv_file.stem
            
            # Detect condition names
            all_cols = df.columns
            conditions = sorted(list({c.split('_')[0] for c in all_cols if '_' in c}))
            
            if reference_condition not in conditions:
                print(f"Warning: Skipping {assay_name} - reference condition '{reference_condition}' not found. Found: {conditions}")
                continue
            
            # Process each row (usually just one row per CSV)
            for idx, row in df.iterrows():
                ssmd_results = calculate_ssmd_for_row(row, conditions, reference_condition)
                
                result_row = {
                    'Assay': assay_name,
                    **ssmd_results
                }
                all_results.append(result_row)
        
        except Exception as e:
            print(f"Error processing {csv_file.name}: {e}")
            continue
    
    if not all_results:
        print("Warning: No valid results generated")
        return pd.DataFrame()
    
    results_df = pd.DataFrame(all_results)
    # Ensure order is preserved by resetting index (DataFrame should already preserve order, but this ensures it)
    results_df = results_df.reset_index(drop=True)
    
    # Debug: Print order to verify
    if file_order is not None:
        print(f"File order requested: {file_order}")
        print(f"Assay order in results: {results_df['Assay'].tolist()}")
    
    return results_df


def create_visualization(results_df, output_path=None, output_filename=None):
    """
    Create horizontal bar chart visualization of SSMD values.
    
    Args:
        results_df: DataFrame with Assay column and SSMD columns
        output_path: Optional directory path to save the plot
        output_filename: Optional custom filename (without extension). If None, uses "SSMD_visualization"
    """
    if results_df.empty:
        print("No data to visualize")
        return
    
    # Set Arial font for all text
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']
    
    # Define fixed bar order: HG, HL, Combo
    bar_order = ['HG', 'HL', 'Combo']
    
    # Color scheme matching reference
    colors = {
        'HG': (0.99, 0.57, 0),      # orange
        'HL': (0.99, 0.14, 0),       # red-orange
        'Combo': (0.57, 0.125, 0.57) # purple
    }
    
    # Find which SSMD columns exist in the data
    ssmd_columns = [col for col in results_df.columns if col.startswith('SSMD_') and col.endswith('_vs_Normal')]
    
    # Filter bar_order to only include conditions that exist in the data
    available_conditions = []
    available_cols = []
    available_colors = []
    
    for cond in bar_order:
        col_name = f'SSMD_{cond}_vs_Normal'
        if col_name in ssmd_columns:
            available_conditions.append(cond)
            available_cols.append(col_name)
            available_colors.append(colors[cond])
    
    # If there are other conditions beyond HG, HL, Combo, add them with default colors
    for col in ssmd_columns:
        cond = col.replace('SSMD_', '').replace('_vs_Normal', '')
        if cond not in available_conditions:
            available_conditions.append(cond)
            available_cols.append(col)
            # Assign a default color (gray-ish)
            available_colors.append((0.5, 0.5, 0.5))
    
    if not available_cols:
        print("No SSMD columns found for visualization")
        return
    
    # Prepare data for plotting
    n_assays = len(results_df)
    bar_height = 0.25
    # Reverse y_positions so first row appears at top (matplotlib barh puts y=0 at bottom)
    y_positions = np.arange(n_assays)[::-1]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6, n_assays * 0.4)))
    
    # Plot each SSMD group in fixed order
    # Note: values are in DataFrame order (first to last), y_positions are reversed (last to first)
    # So first DataFrame row will be at top of plot
    for i, (col, color) in enumerate(zip(available_cols, available_colors)):
        cond_name = col.replace('SSMD_', '').replace('_vs_Normal', '')
        # Convert to list to ensure order is preserved (pandas Series might have index issues)
        values = results_df[col].tolist()
        ax.barh(
            y_positions + (i - 1) * bar_height,
            values,
            height=bar_height,
            color=color,
            label=cond_name
        )
    
    # Calculate symmetric x-axis range
    all_ssmd_values = []
    for col in available_cols:
        all_ssmd_values.extend(results_df[col].dropna().tolist())
    
    if all_ssmd_values:
        max_abs_value = max(abs(min(all_ssmd_values)), abs(max(all_ssmd_values)))
        # Add some padding (10% of the range)
        padding = max_abs_value * 0.1
        x_limit = max_abs_value + padding
        
        # Ensure x_limit is even (if odd, round up to next even number)
        if x_limit % 2 != 0:
            x_limit = int(np.ceil(x_limit / 2)) * 2
        else:
            x_limit = int(x_limit)
        
        ax.set_xlim(-x_limit, x_limit)
        
        # Set x-axis ticks to show only even numbers
        # Generate even numbers from -x_limit to x_limit
        even_ticks = list(range(-x_limit, x_limit + 1, 2))
        ax.set_xticks(even_ticks)
    
    # Set labels and formatting
    ax.set_title("SSMD Scores", fontfamily='Arial', fontsize=17, fontweight='bold')
    ax.set_yticks(y_positions)
    # Labels should match the order: first DataFrame row at top
    # y_positions are reversed [n-1, n-2, ..., 0], so labels should be in DataFrame order
    assay_labels = results_df['Assay'].tolist()
    ax.set_yticklabels(assay_labels, fontfamily='Arial', fontweight='bold')
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("SSMD", fontfamily='Arial', fontweight='bold')
    ax.legend(title="Groups", prop={'family': 'Arial'})  # Legend remains normal weight
    
    # Set tick labels to Arial and bold
    for label in ax.get_xticklabels():
        label.set_fontfamily('Arial')
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontfamily('Arial')
        label.set_fontweight('bold')
    
    plt.tight_layout()
    
    # Save plot if output_path is provided
    if output_path:
        # Determine filename
        if output_filename:
            # Use custom filename
            filename = output_filename if output_filename.endswith('.png') else f"{output_filename}.png"
            save_path = Path(output_path) / filename
        else:
            # Use default filename
            save_path = Path(output_path) / "SSMD_visualization.png"
        
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    # Display plot
    plt.show()


def main():
    """Main function to process CSVs and create visualization."""
    # Get project root (assuming script is in backend folder)
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    processed_data_dir = project_root / 'data' / 'processed_data'
    
    print("="*80)
    print("SSMD VISUALIZATION")
    print("="*80)
    print(f"Processing CSVs from: {processed_data_dir}")
    
    # Process all CSVs
    results_df = process_all_csvs(processed_data_dir, reference_condition='Normal')
    
    if results_df.empty:
        print("No data to visualize. Exiting.")
        return
    
    print(f"\nProcessed {len(results_df)} assays")
    print(f"Found SSMD columns: {[col for col in results_df.columns if col.startswith('SSMD_')]}")
    
    # Create visualization
    create_visualization(results_df, output_path=processed_data_dir, output_filename=None)
    
    print("="*80)


if __name__ == "__main__":
    main()

