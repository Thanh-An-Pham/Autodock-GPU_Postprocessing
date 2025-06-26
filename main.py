import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve
import shutil

def calculate_logAUC(y_true, y_pred, alpha=0.001):
    """
    Calculate logAUC score.
    
    Args:
        y_true: True labels (1 for active, 0 for decoy)
        y_pred: Predicted scores
        alpha: Minimum threshold for logAUC calculation
        
    Returns:
        float: logAUC score
    """
    # Sort predictions in descending order
    sorted_idx = np.argsort(-y_pred)
    y_true_sorted = y_true[sorted_idx]
    y_pred_sorted = y_pred[sorted_idx]
    
    # Calculate cumulative sum of actives
    cumsum = np.cumsum(y_true_sorted)
    
    # Calculate logAUC
    logAUC = 0
    for i in range(len(y_pred_sorted)):
        if y_pred_sorted[i] >= alpha:
            logAUC += cumsum[i] / (i + 1)
    
    return logAUC / len(y_pred_sorted)

def calculate_EF1(y_true, y_pred):
    """
    Calculate EF1% (Enrichment Factor at 1%).
    
    Args:
        y_true: True labels (1 for active, 0 for decoy)
        y_pred: Predicted scores
        
    Returns:
        float: EF1% score
    """
    # Sort predictions in descending order
    sorted_idx = np.argsort(-y_pred)
    y_true_sorted = y_true[sorted_idx]
    
    # Calculate number of compounds to consider (1% of total)
    n_compounds = int(len(y_true) * 0.01)
    if n_compounds == 0:
        n_compounds = 1
    
    # Calculate number of actives in top 1%
    n_actives_top = np.sum(y_true_sorted[:n_compounds])
    
    # Calculate total number of actives
    n_actives_total = np.sum(y_true)
    
    # Calculate EF1%
    if n_actives_total == 0:
        return 0
    return (n_actives_top / n_compounds) / (n_actives_total / len(y_true))

def get_labels(stats_df):
    """
    Get labels for compounds based on their names.
    Returns 0 for inactive/decoy compounds, 1 for active compounds.
    """
    return np.array([0 if ('decoy' in idx.lower() or 'inactive' in idx.lower()) else 1 for idx in stats_df.index])

def plot_combined_roc_curves(stats_df, save_path):
    """
    Plot ROC curves for all columns on the same figure.
    
    Args:
        stats_df: DataFrame containing the columns to analyze
        save_path: Path to save the plot
    Returns:
        dict: Dictionary containing G-mean values and cutoff points for each metric
    """
    y_true = get_labels(stats_df)
    
    plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(stats_df.columns)))
    
    g_means_dict = {}
    
    for i, col in enumerate(stats_df.columns):
        y_pred = -stats_df[col].values
        fpr, tpr, thresholds = roc_curve(y_true, y_pred)
        roc_auc = roc_auc_score(y_true, y_pred)
        
        # Calculate G-mean
        g_means = np.sqrt(tpr * (1 - fpr))
        # Find optimal threshold
        optimal_idx = np.argmax(g_means)
        optimal_threshold = thresholds[optimal_idx]
        optimal_g_mean = g_means[optimal_idx]
        
        # Store G-mean values and cutoff
        g_means_dict[col] = {
            'g_mean': optimal_g_mean,
            'threshold': optimal_threshold,
            'tpr': tpr[optimal_idx],
            'fpr': fpr[optimal_idx],
            'cutoff': -optimal_threshold  # Convert back to original scale
        }
        
        # Plot ROC curve
        plt.plot(fpr, tpr, color=colors[i], lw=2, 
                label=f'{col} (AUC = {roc_auc:.3f}, Cutoff = {-optimal_threshold:.2f})')
        # Plot optimal point
        plt.scatter(fpr[optimal_idx], tpr[optimal_idx], 
                   color=colors[i], marker='o')
    
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=16)
    plt.ylabel('True Positive Rate', fontsize=16)
    plt.title('ROC Curves Comparison', fontsize=24, weight='semibold')
    plt.legend(loc="lower right")
    plt.grid(True, alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return g_means_dict

def plot_combined_logAUC_curves(stats_df, save_path):
    """
    Plot logAUC curves for all columns on the same figure.
    
    Args:
        stats_df: DataFrame containing the columns to analyze
        save_path: Path to save the plot
    """
    y_true = get_labels(stats_df)
    
    plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(stats_df.columns)))
    
    for i, col in enumerate(stats_df.columns):
        y_pred = -stats_df[col].values
        sorted_idx = np.argsort(-y_pred)
        y_true_sorted = y_true[sorted_idx]
        cumsum = np.cumsum(y_true_sorted)
        total_actives = np.sum(y_true)
        
        x = np.arange(1, len(y_pred) + 1) / len(y_pred)
        y = cumsum / total_actives
        
        plt.plot(x, y, color=colors[i], lw=2, label=f'{col}')
    
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random')
    plt.xscale('log')
    plt.xlim([0.001, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Log(Fraction of Database)')
    plt.ylabel('Fraction of Actives Found')
    plt.title('logAUC Curves Comparison')
    plt.legend(loc="lower right")
    plt.grid(True, alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_combined_enrichment_curves(stats_df, save_path):
    """
    Plot enrichment curves for all columns on the same figure.
    
    Args:
        stats_df: DataFrame containing the columns to analyze
        save_path: Path to save the plot
    """
    y_true = get_labels(stats_df)
    
    plt.figure(figsize=(10, 8))
    colors = plt.cm.tab10(np.linspace(0, 1, len(stats_df.columns)))
    
    percentages = np.arange(1, 101) / 100
    
    for i, col in enumerate(stats_df.columns):
        y_pred = -stats_df[col].values
        sorted_idx = np.argsort(-y_pred)
        y_true_sorted = y_true[sorted_idx]
        
        enrichment_factors = []
        for p in percentages:
            n_compounds = int(len(y_true) * p)
            if n_compounds == 0:
                n_compounds = 1
            n_actives_top = np.sum(y_true_sorted[:n_compounds])
            n_actives_total = np.sum(y_true)
            if n_actives_total == 0:
                ef = 0
            else:
                ef = (n_actives_top / n_compounds) / (n_actives_total / len(y_true))
            enrichment_factors.append(ef)
        
        plt.plot(percentages * 100, enrichment_factors, color=colors[i], lw=2, label=f'{col}')
    
    plt.axhline(y=1, color='navy', linestyle='--', label='Random')
    plt.xlim([0, 100])
    plt.ylim([0, max([max(plt.gca().lines[i].get_ydata()) for i in range(len(plt.gca().lines))]) * 1.1])
    plt.xlabel('Percentage of Database Screened')
    plt.ylabel('Enrichment Factor')
    plt.title('Enrichment Curves Comparison')
    plt.legend(loc="upper right")
    plt.grid(True, alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def analyze_metrics(final_csv_path, output_dir, is_training=False, train_cutoffs=None, pdbqt_dir=None):
    """
    Analyze metrics and create plots from column stats CSV.
    
    Args:
        final_csv_path: Path to the column stats CSV file
        output_dir: Directory to save the metrics and plots
        is_training: Whether this is the training set
        train_cutoffs: Dictionary of cutoffs from training set (required for test set)
        pdbqt_dir: Directory containing PDBQT files (required for test set)
    """
    # Create output directory if it doesn't exist
    if is_training:
        output_dir = os.path.join(output_dir, 'train_metrics')
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the stats CSV
    stats_df = pd.read_csv(final_csv_path, index_col=0)
    
    if is_training:
        # For training set: calculate metrics and cutoffs
        y_true = get_labels(stats_df)
        
        # Print debug information
        # print("\nLabel distribution:")
        print(f"Total compounds: {len(y_true)}")
        print(f"Active compounds: {np.sum(y_true)}")
        print(f"Inactive compounds: {len(y_true) - np.sum(y_true)}")
        # print("\nExample labels:")
        # for idx, label in zip(stats_df.index[:5], y_true[:5]):
        #     print(f"{idx}: {'Active' if label == 1 else 'Inactive'}")
        
        if len(np.unique(y_true)) < 2:
            print("Warning: No inactive compounds found in the dataset. Metrics may not be meaningful.")
            metrics = []
            for col in stats_df.columns:
                metrics.append({
                    'metric': col,
                    'ROC-AUC': np.nan,
                    'logAUC': np.nan,
                    'EF1%': np.nan,
                    'cutoff': np.nan
                })
            metrics_df = pd.DataFrame(metrics)
        else:
            metrics = []
            # print("\nAnalyzing each metric:")
            #For each scoring metric in the data (adgpu_min, adgpu_max, adgpu_median):
            for col in stats_df.columns:
                # print(f"\nProcessing {col}:")
                y_pred = -stats_df[col].values
                
                # Check for NaN values
                nan_mask = np.isnan(y_pred)
                if np.any(nan_mask):
                    print(f"  Warning: Found {np.sum(nan_mask)} NaN values in {col}")
                    y_pred = y_pred[~nan_mask]
                    y_true_col = y_true[~nan_mask]
                else:
                    y_true_col = y_true
                
                # Check class distribution
                unique_classes = np.unique(y_true_col)
                # print(f"  Classes present: {unique_classes}")
                # print(f"  Active compounds: {np.sum(y_true_col == 1)}")
                # print(f"  Inactive compounds: {np.sum(y_true_col == 0)}")
                
                if len(unique_classes) < 2:
                    print(f"  Warning: Only one class present in {col}, skipping metrics")
                    metrics.append({
                        'metric': col,
                        'ROC-AUC': np.nan,
                        'logAUC': np.nan,
                        'EF1%': np.nan,
                        'cutoff': np.nan
                    })
                    continue
                
                try:
                    # Calculate metrics
                    roc_auc = roc_auc_score(y_true_col, y_pred)
                    log_auc = calculate_logAUC(y_true_col, y_pred)
                    ef1 = calculate_EF1(y_true_col, y_pred)
                    
                    # Calculate optimal cutoff using ROC curve
                    fpr, tpr, thresholds = roc_curve(y_true_col, y_pred)
                    g_means = np.sqrt(tpr * (1 - fpr))
                    optimal_idx = np.argmax(g_means)
                    optimal_threshold = -thresholds[optimal_idx]
                    
                    # print(f"  ROC-AUC: {roc_auc:.3f}")
                    # print(f"  logAUC: {log_auc:.3f}")
                    # print(f"  EF1%: {ef1:.3f}")
                    # print(f"  Cutoff: {optimal_threshold:.3f}")
                    
                    metrics.append({
                        'metric': col,
                        'ROC-AUC': roc_auc,
                        'logAUC': log_auc,
                        'EF1%': ef1,
                        'cutoff': optimal_threshold
                    })
                except Exception as e:
                    print(f"  Error calculating metrics for {col}: {str(e)}")
                    metrics.append({
                        'metric': col,
                        'ROC-AUC': np.nan,
                        'logAUC': np.nan,
                        'EF1%': np.nan,
                        'cutoff': np.nan
                    })
            
            metrics_df = pd.DataFrame(metrics)
            
            # Create plots for training set
            plot_combined_roc_curves(stats_df, os.path.join(output_dir, 'roc_curves.png'))
            plot_combined_logAUC_curves(stats_df, os.path.join(output_dir, 'logAUC_curves.png'))
            plot_combined_enrichment_curves(stats_df, os.path.join(output_dir, 'enrichment_curves.png'))
            
            # Save metrics to CSV
            metrics_df.to_csv(os.path.join(output_dir, 'metrics.csv'), index=False)
            
            return metrics_df
            
    else:
        # For test set: only apply cutoffs to identify active compounds
        if train_cutoffs is None:
            raise ValueError("train_cutoffs must be provided for test set analysis")
        if pdbqt_dir is None:
            raise ValueError("pdbqt_dir must be provided for test set analysis")
        
        # Select the best metric based on highest ROC-AUC value
        best_metric = None
        best_auc = -1
        best_cutoff = None
        
        # Print available metrics and find the one with highest ROC-AUC
        print("\nAvailable metrics from training set:")
        
        # Check if we received a DataFrame or a dictionary
        if isinstance(train_cutoffs, pd.DataFrame):
            # We have a full DataFrame with metrics
            metrics_df = train_cutoffs
            
            # Print available metrics
            for _, row in metrics_df.iterrows():
                metric = row['metric']
                cutoff = row['cutoff']
                auc = row['ROC-AUC']
                print(f"  {metric}: cutoff={cutoff:.3f}, ROC-AUC={auc:.3f}")
            
            # Find the metric with highest ROC-AUC using pandas operations
            if 'ROC-AUC' in metrics_df.columns:
                # Get the index of the row with the maximum ROC-AUC value
                best_idx = metrics_df['ROC-AUC'].idxmax()
                
                # Get the metric, cutoff, and AUC from that row
                best_metric = metrics_df.loc[best_idx, 'metric']
                best_cutoff = metrics_df.loc[best_idx, 'cutoff']
                best_auc = metrics_df.loc[best_idx, 'ROC-AUC']
            else:
                raise ValueError("ROC-AUC column not found in metrics DataFrame")
        else:
            # We have a dictionary mapping metrics to cutoffs
            # Get the corresponding metrics DataFrame from the same location
            train_metrics_path = os.path.join(os.path.dirname(output_dir), "metrics.csv")
            
            if os.path.exists(train_metrics_path):
                metrics_df = pd.read_csv(train_metrics_path)
                
                # Print available metrics
                valid_metrics = []
                for _, row in metrics_df.iterrows():
                    metric = row['metric']
                    if metric in train_cutoffs:  # Make sure this metric has a cutoff
                        cutoff = train_cutoffs[metric]
                        auc = row['ROC-AUC']
                        print(f"  {metric}: cutoff={cutoff:.3f}, ROC-AUC={auc:.3f}")
                        valid_metrics.append(metric)
                
                # Filter metrics to only include those with cutoffs
                filtered_df = metrics_df[metrics_df['metric'].isin(valid_metrics)]
                
                if not filtered_df.empty and 'ROC-AUC' in filtered_df.columns:
                    # Get the index of the row with the maximum ROC-AUC value
                    best_idx = filtered_df['ROC-AUC'].idxmax()
                    
                    # Get the metric, cutoff, and AUC from that row
                    best_metric = metrics_df.loc[best_idx, 'metric']
                    best_cutoff = train_cutoffs[best_metric]  # Get from dictionary
                    best_auc = metrics_df.loc[best_idx, 'ROC-AUC']
                else:
                    raise ValueError("No valid metrics with ROC-AUC values found")
            else:
                raise ValueError(f"Metrics file not found at {train_metrics_path}")
        
        # Confirm we found a best metric
        if best_metric is None:
            raise ValueError("Could not determine best metric based on ROC-AUC values")
            
        print(f"\nSelected best metric: {best_metric} (ROC-AUC={best_auc:.3f}, cutoff={best_cutoff:.3f})")
        
        # Get only the active compounds that meet the cutoff for the best metric
        if best_metric in stats_df.columns:
            # Apply cutoff to identify active compounds
            active_mask = stats_df[best_metric] <= best_cutoff
            active_compounds = stats_df.index[active_mask].tolist()
            
            print(f"Found {len(active_compounds)} active compounds using {best_metric} (cutoff: {best_cutoff:.3f})")
            
            # Save list of active compounds to CSV
            active_df = pd.DataFrame({
                'compound': active_compounds,
                'score': stats_df.loc[active_compounds, best_metric].values
            })
            active_df = active_df.sort_values('score')  # Sort by binding energy (ascending)
            active_df.to_csv(os.path.join(output_dir, 'active_compounds.csv'), index=False)
            
            # No longer extract PDBQT files here - this is now done in process_datasets.py
            
            return active_df
        else:
            print(f"Error: Metric '{best_metric}' not found in test set columns: {stats_df.columns.tolist()}")
            return pd.DataFrame()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze metrics and create plots from column stats CSV')
    parser.add_argument('--final_csv_path', help='Path to the column stats CSV file')
    parser.add_argument('--output_dir', default='metrics_output', 
                       help='Directory to save the metrics and plots (default: metrics_output)')
    parser.add_argument('--is_training', action='store_true',
                       help='Whether this is the training set')
    parser.add_argument('--train_cutoffs', help='Path to training cutoffs CSV file (required for test set)')
    parser.add_argument('--pdbqt_dir', help='Directory containing PDBQT files (required for test set)')
    
    args = parser.parse_args()
    
    train_cutoffs = None
    if not args.is_training:
        if args.train_cutoffs is None:
            raise ValueError("--train_cutoffs must be provided for test set analysis")
        if args.pdbqt_dir is None:
            raise ValueError("--pdbqt_dir must be provided for test set analysis")
        train_cutoffs = pd.read_csv(args.train_cutoffs, index_col=0)['cutoff'].to_dict()
    
    analyze_metrics(args.final_csv_path, args.output_dir, args.is_training, train_cutoffs, args.pdbqt_dir) 