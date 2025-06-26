import os
import glob
import re
import pandas as pd
import shutil

def _parse_dlg(lines):
    """
    Single-pass scan of .dlg file lines to extract:
      - clustering rows: (cluster_rank, num_in_cluster, run)
      - rmsd rows:      (cluster_rank, binding_energy)
    """
    clust_data = []
    rmsd_data  = []
    in_clust = in_rmsd = False

    for line in lines:
        # section switch
        if 'CLUSTERING HISTOGRAM' in line:
            in_clust, in_rmsd = True, False
            continue
        if 'RMSD TABLE' in line:
            in_clust, in_rmsd = False, True
            continue

        if in_clust:
            # Cluster rank | Cluster size | Run | Cluster RMSD
            # ^: Start of line, \s: whitespace, *: 0 or more, \d+: digit one or more
            # End with |
            if '|' in line and re.match(r'^\s*\d+\s*\|', line):
                parts = [p.strip() for p in line.split('|')]
                try:
                    clust_data.append({
                        'cluster_rank':    int(parts[0]),
                        'run':             int(parts[2]),
                        'num_in_cluster':  int(parts[4]),
                    })
                except (IndexError, ValueError):
                    pass

        elif in_rmsd:
            cols = line.split()
            if len(cols) >= 4 and cols[0].isdigit():
                try:
                    rmsd_data.append({
                        'cluster_rank':   int(cols[0]),
                        'binding_energy': float(cols[3])
                    })
                except ValueError:
                    pass

    return clust_data, rmsd_data


def process_all_dlg(input_dir, csv_dir):
    """
    1) Parses every .dlg in input_dir
    2) Writes one CSV per ligand under csv_dir
    3) Returns dict of DataFrames
    4) Filters out compounds with max binding value > 3 and moves them to a fail_log_dlg folder
    """
    os.makedirs(csv_dir, exist_ok=True)
    
    # Create a folder for failed DLG files
    fail_dir = os.path.join(os.path.dirname(input_dir), "fail_log_dlg")
    os.makedirs(fail_dir, exist_ok=True)
    
    results = {}
    filtered_compounds = []

    dlg_paths = sorted(glob.glob(os.path.join(input_dir, '*.dlg')))
    for dlg in dlg_paths:
        base = os.path.splitext(os.path.basename(dlg))[0]
        lines = open(dlg).read().splitlines()

        clust_data, rmsd_data = _parse_dlg(lines)
        if not clust_data or not rmsd_data:
            continue

        df_clust = pd.DataFrame(clust_data)
        max_size = df_clust['num_in_cluster'].max()
        best = (
            df_clust[df_clust['num_in_cluster'] == max_size]
            .nsmallest(1, 'cluster_rank')
            .iloc[0]
        )
        rk, run = best['cluster_rank'], best['run']

        df_rmsd = pd.DataFrame(rmsd_data)
        df_best = df_rmsd[df_rmsd['cluster_rank'] == rk] \
                     .reset_index(drop=True) \
                     .rename(columns={'binding_energy': f"{base}_{rk}_{run}"})

        # Check if the maximum binding value is > 3
        # Note: binding energies are negative, so we're looking for binding values > 3
        column_name = f"{base}_{rk}_{run}"
        if df_best[column_name].max() > 3:
            filtered_compounds.append(base)
            # Move the DLG file to the fail_log_dlg folder
            fail_dest = os.path.join(fail_dir, os.path.basename(dlg))
            try:
                shutil.move(dlg, fail_dest)
                print(f"Moved {base} to fail_log_dlg (max binding value > 3)")
            except Exception as e:
                print(f"Warning: Could not move DLG file for {base}: {e}")
            continue

        name = f"{base}_{rk}_{run}"
        results[name] = df_best[[f"{base}_{rk}_{run}"]]

        # write individual CSV
        df_best.to_csv(os.path.join(csv_dir, f"{name}.csv"),
                       index=False, encoding='utf-8-sig')

    print(f"Filtered out {len(filtered_compounds)} compounds with binding values > 3")
    if filtered_compounds:
        print(f"Failed DLG files moved to: {fail_dir}")
    return results


def concat_csvs(csv_dir, output_path):
    """
    Reads all CSVs in csv_dir, concatenates them side-by-side,
    and writes the combined DataFrame to output_path.
    """
    paths = sorted(glob.glob(os.path.join(csv_dir, '*.csv')))
    if not paths:
        raise FileNotFoundError(f"No CSVs found in {csv_dir!r}")

    df_list = []
    for p in paths:
        df = pd.read_csv(p)
        # Ensure there's at least two columns
        if df.shape[1] < 2:
            raise ValueError(f"File {os.path.basename(p)!r} has fewer than 2 columns")
        # Select only the 2nd column and preserve its header
        col = df.columns[1]
        # Convert to numeric, coercing errors to NaN
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df_list.append(df[[col]])

    final_df = pd.concat(df_list, axis=1)
    final_df.to_csv(output_path, index=False, encoding='utf-8-sig')
    return final_df


def get_column_stats(df):
    """
    Calculate minimum, maximum, and median values for each column in the DataFrame.
    Returns a new DataFrame with these statistics where each row is a compound
    and columns are adgpu_min, adgpu_max, and adgpu_median.
    
    Args:
        df (pd.DataFrame): Input DataFrame with numeric columns
        
    Returns:
        pd.DataFrame: DataFrame containing statistics for each compound
    """
    # Convert DataFrame to numeric, coercing errors to NaN
    numeric_df = df.apply(pd.to_numeric, errors='coerce')
    
    stats = pd.DataFrame({
        'adgpu_min': numeric_df.min(),
        'adgpu_max': numeric_df.max(),
        'adgpu_median': numeric_df.median()
    })
    return stats
