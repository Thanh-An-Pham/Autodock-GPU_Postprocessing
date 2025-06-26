import os
import pandas as pd
import argparse
import tempfile
from autodock_postprocessing import process_all_dlg, concat_csvs, get_column_stats
from main import analyze_metrics
from screening_pdbqt import main as extract_pdbqt, extract_pdbqt_for_actives_list
from convert_pdbqt_to_pdb import convert_pdbqt_to_pdb

def process_dataset(dlg_dir, output_dir, dataset_type):
    """
    Process a dataset (train or test) of DLG files.
    
    Args:
        dlg_dir: Directory containing DLG files
        output_dir: Directory to save output files
        dataset_type: Either "train" or "test"
    """
    # Create output directories
    metrics_dir = os.path.join(output_dir, f"{dataset_type}_metrics")
    os.makedirs(metrics_dir, exist_ok=True)
    
    # Create a temporary directory for intermediate CSV files
    temp_csv_dir = os.path.join(output_dir, f"{dataset_type}_temp_csv")
    os.makedirs(temp_csv_dir, exist_ok=True)
    
    # Process DLG files using the correct functions from autodock_postprocessing
    # 1. Process all DLG files to individual CSVs
    process_all_dlg(dlg_dir, temp_csv_dir)
    
    # 2. Concatenate all CSVs
    combined_csv = os.path.join(metrics_dir, f"{dataset_type}_combined.csv")
    combined_df = concat_csvs(temp_csv_dir, combined_csv)
    
    # 3. Calculate statistics
    stats_df = get_column_stats(combined_df)
    stats_csv = os.path.join(metrics_dir, f"{dataset_type}_stats.csv")
    stats_df.to_csv(stats_csv)
    
    # Analyze metrics
    if dataset_type == "train":
        # For training set, calculate metrics and cutoffs
        print("Calculating metrics and finding optimal cutoffs from training set...")
        metrics_results = analyze_metrics(stats_csv, metrics_dir, is_training=True)
        return metrics_results, stats_csv
    else:
        # For test set, we need the training cutoffs
        # Note: The metrics.csv file is saved within the train_metrics directory by analyze_metrics
        train_metrics_path = os.path.join(output_dir, "train_metrics", "train_metrics", "metrics.csv")
        if not os.path.exists(train_metrics_path):
            # Try alternate path as fallback
            alternate_path = os.path.join(output_dir, "train_metrics", "metrics.csv")
            if os.path.exists(alternate_path):
                train_metrics_path = alternate_path
            else:
                print(f"Error: Could not find training metrics file at {train_metrics_path} or {alternate_path}")
                print(f"Directory contents of {os.path.join(output_dir, 'train_metrics')}:")
                try:
                    print(os.listdir(os.path.join(output_dir, "train_metrics")))
                except:
                    print("Could not list directory contents")
                raise ValueError("Training metrics file not found. Please process training set first.")
        
        print(f"Loading training metrics from: {train_metrics_path}")
        # Load training metrics DataFrame with ROC-AUC values
        train_metrics_df = pd.read_csv(train_metrics_path)
        
        # Process test set using training metrics
        metrics_results = analyze_metrics(
            stats_csv, 
            metrics_dir, 
            is_training=False,
            train_cutoffs=train_metrics_df,  # Pass the full DataFrame with ROC-AUC values
            pdbqt_dir=dlg_dir  # Source directory for PDBQT extraction
        )
        
        # Try to locate active PDBQT directory or PDBQT files
        active_pdbqt_dir = os.path.join(metrics_dir, 'active_pdbqt')
        os.makedirs(active_pdbqt_dir, exist_ok=True)
        
        # Get active compounds identified by metrics analysis
        active_compounds_file = os.path.join(metrics_dir, 'active_compounds.csv')
        if os.path.exists(active_compounds_file):
            # If active compounds list exists, read it directly
            active_df = pd.read_csv(active_compounds_file)
            active_compounds = active_df['compound'].tolist()
        else:
            # If for some reason the file wasn't created, print a warning
            print("Warning: Active compounds list not found. Please check metrics analysis results.")
            return metrics_results, stats_csv
            
        # Direct extraction from DLG files to PDBQT
        extract_success = extract_pdbqt_for_actives_list(active_compounds, dlg_dir, active_pdbqt_dir)
        
        # Convert active PDBQT files to PDB if any were successfully extracted
        if extract_success:
            print(f"Converting active compounds from PDBQT to PDB...")
            convert_pdbqt_to_pdb(active_pdbqt_dir)
            print(f"PDB files saved to: {active_pdbqt_dir}")
        else:
            print("No active PDBQT files were successfully extracted. Check DLG files and compound names.")
        
        return metrics_results, stats_csv

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process docking results from DLG files')
    parser.add_argument('--train_dlg_dir', required=True, help='Directory containing training DLG files')
    parser.add_argument('--test_dlg_dir', required=True, help='Directory containing test DLG files')
    parser.add_argument('--output_dir', default='output', help='Directory to save output files')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Process training set first
    print("Processing training set...")
    train_metrics, train_stats = process_dataset(args.train_dlg_dir, args.output_dir, "train")
    
    # Then process test set
    print("\nProcessing test set...")
    test_metrics, test_stats = process_dataset(args.test_dlg_dir, args.output_dir, "test")
    
    print("\nProcessing complete!")
    print(f"Training metrics saved to: {os.path.join(args.output_dir, 'train_metrics')}")
    print(f"Test active compounds saved to: {os.path.join(args.output_dir, 'test_metrics')}")

if __name__ == "__main__":
    main() 