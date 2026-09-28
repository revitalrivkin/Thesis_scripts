import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Set up path to import SSMD calculation function from backend
script_dir = Path(__file__).resolve().parent
# Navigate up to "Data analysis projects" directory
data_analysis_projects_dir = script_dir.parents[2]
# Path to SSMD_GUI backend
backend_dir = data_analysis_projects_dir / "SSMD_GUI" / "backend"
sys.path.insert(0, str(backend_dir))

from visualize_ssmd import calculate_ssmd_for_row


def load_dep_data(csv_path):
    """
    Load DEP MS data CSV and validate required columns exist.
    
    Args:
        csv_path: Path to DEP_MS_data.csv
    
    Returns:
        pandas.DataFrame: Loaded data with validated columns
    """
    df = pd.read_csv(csv_path)
    
    # Expected columns
    required_cols = ['Genes'] + [
        f'{cond}-{i}' for cond in ['Normal', 'HG', 'HL', 'Combo'] 
        for i in [1, 2, 3]
    ]
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Coerce replicate columns to numeric
    replicate_cols = [col for col in df.columns if col != 'Genes']
    for col in replicate_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df


def compute_summary_stats_for_row(row, conditions=['Normal', 'HG', 'HL', 'Combo']):
    """
    Compute mean, SD, and N for each condition from replicate columns.
    
    Args:
        row: pandas Series representing one assay row
        conditions: List of condition names
    
    Returns:
        pandas.Series: Summary stats with columns like Normal_mean, Normal_SD, Normal_N, etc.
    """
    summary = {}
    
    for cond in conditions:
        # Get replicate columns for this condition
        replicate_cols = [f'{cond}-{i}' for i in [1, 2, 3]]
        values = []
        
        for col in replicate_cols:
            if col in row.index:
                val = row[col]
                if pd.notna(val):
                    values.append(float(val))
        
        # Compute statistics
        if len(values) == 0:
            mean_val = np.nan
            sd_val = np.nan
            n_val = 0
        elif len(values) == 1:
            mean_val = values[0]
            sd_val = np.nan  # SD undefined with only one value
            n_val = 1
        else:
            mean_val = np.mean(values)
            sd_val = np.std(values, ddof=1)  # Sample standard deviation
            n_val = len(values)
        
        summary[f'{cond}_mean'] = mean_val
        summary[f'{cond}_SD'] = sd_val
        summary[f'{cond}_N'] = n_val
    
    return pd.Series(summary)


def calculate_ssmd_for_dep_assays(df, reference_condition='Normal'):
    """
    Calculate SSMD values for each DEP assay.
    
    Args:
        df: DataFrame with DEP data (columns: Genes, Normal-1, Normal-2, etc.)
        reference_condition: Reference condition name (default 'Normal')
    
    Returns:
        pandas.DataFrame: Results with Assay column and SSMD columns
    """
    conditions = ['Normal', 'HG', 'HL', 'Combo']
    all_results = []
    
    for idx, row in df.iterrows():
        assay_name = row['Genes']
        
        # Compute summary stats for this row
        summary_row = compute_summary_stats_for_row(row, conditions)
        
        # Calculate SSMD using the imported function
        try:
            ssmd_results = calculate_ssmd_for_row(
                summary_row, 
                conditions=conditions, 
                reference_condition=reference_condition
            )
            
            result_row = {
                'Assay': assay_name,
                **ssmd_results
            }
            all_results.append(result_row)
        except Exception as e:
            print(f"Warning: Error calculating SSMD for {assay_name}: {e}")
            continue
    
    if not all_results:
        return pd.DataFrame()
    
    results_df = pd.DataFrame(all_results)
    return results_df


def calculate_mse_from_ssmd(ssmd_df):
    """
    Calculate MSE (mean squared error) for each condition from SSMD values.
    
    Args:
        ssmd_df: DataFrame with SSMD columns (SSMD_HG_vs_Normal, etc.)
    
    Returns:
        pandas.DataFrame: One-row DataFrame with MSE_HG, MSE_HL, MSE_Combo
    """
    mse_results = {}
    
    # SSMD columns of interest
    ssmd_columns = {
        'SSMD_HG_vs_Normal': 'MSE_HG',
        'SSMD_HL_vs_Normal': 'MSE_HL',
        'SSMD_Combo_vs_Normal': 'MSE_Combo'
    }
    
    for ssmd_col, mse_col in ssmd_columns.items():
        if ssmd_col in ssmd_df.columns:
            # Get SSMD values, drop NaNs
            ssmd_values = ssmd_df[ssmd_col].dropna().values
            
            if len(ssmd_values) > 0:
                # MSE = mean(SSMD^2)
                mse = np.mean(ssmd_values ** 2)
                mse_results[mse_col] = mse
            else:
                mse_results[mse_col] = np.nan
        else:
            mse_results[mse_col] = np.nan
    
    # Create one-row DataFrame
    mse_df = pd.DataFrame([mse_results])
    return mse_df


def main():
    """
    Main function to process DEP data and calculate SSMD and MSE.
    """
    # Locate input file relative to script directory
    script_dir = Path(__file__).resolve().parent
    input_csv = script_dir / "DEP_MS_data.csv"

    if not input_csv.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_csv}\n"
            "Make sure 'DEP_MS_data.csv' is in the same folder as 'MSE_of_DEPs.py'."
        )
    
    print("="*80)
    print("DEP SSMD and MSE Calculation")
    print("="*80)
    print(f"Loading data from: {input_csv}")
    
    # Load data
    df = load_dep_data(input_csv)
    print(f"Loaded {len(df)} assays")
    
    # Calculate SSMD for each assay
    print("\nCalculating SSMD values for each assay...")
    ssmd_results_df = calculate_ssmd_for_dep_assays(df, reference_condition='Normal')
    
    if ssmd_results_df.empty:
        print("Error: No SSMD results generated")
        return
    
    print(f"Calculated SSMD for {len(ssmd_results_df)} assays")
    
    # Save SSMD results
    ssmd_output = script_dir / "DEP_SSMD_results.csv"
    ssmd_results_df.to_csv(ssmd_output, index=False)
    print(f"SSMD results saved to: {ssmd_output}")
    
    # Calculate MSE
    print("\nCalculating MSE values...")
    mse_df = calculate_mse_from_ssmd(ssmd_results_df)
    
    # Display MSE results
    print("\n" + "="*80)
    print("MSE RESULTS")
    print("="*80)
    print(mse_df.to_string(index=False))
    print("="*80)
    
    # Save MSE results
    mse_output = script_dir / "DEP_MSE_results.csv"
    mse_df.to_csv(mse_output, index=False)
    print(f"\nMSE results saved to: {mse_output}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
