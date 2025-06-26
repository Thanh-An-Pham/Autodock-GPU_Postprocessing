import re
import sys
import os

def parse_cluster_table(lines):
    """
    Parse the CLUSTERING HISTOGRAM section to find the run number for the cluster with the most members.
    """
    in_cluster = False
    best_run = None
    max_members = -1
    for line in lines:
        if 'CLUSTERING HISTOGRAM' in line:
            in_cluster = True
            continue
        if in_cluster:
            # End of cluster table when RMSD TABLE starts
            if 'RMSD TABLE' in line:
                break
            # Match lines like: "   1 |     -8.63 | 401 |     -7.96 | 987 |"
            m = re.match(r"\s*(\d+)\s*\|\s*([\-\d\.]+)\s*\|\s*(\d+)\s*\|\s*([\-\d\.]+)\s*\|\s*(\d+)", line)
            if m:
                # group(1)=rank, group(2)=lowest binding, group(3)=run, group(4)=mean binding, group(5)=num in cluster
                run = int(m.group(3))
                members = int(m.group(5))
                if members > max_members:
                    max_members = members
                    best_run = run
    return best_run


def extract_run_section(lines, run_number):
    """
    Extract lines from 'DOCKED: USER    Run = {run_number}' until 'DOCKED: ENDMDL'.
    """
    pattern = f"DOCKED: USER    Run = {run_number}"
    section = []
    capturing = False
    
    for line in lines:
        if capturing:
            section.append(line)
            if 'DOCKED: ENDMDL' in line:
                break
        elif pattern in line:
            section.append(line)
            capturing = True
    return section


def main(dlg_path, output_path):
    # Handle both single file and directory
    if os.path.isdir(dlg_path):
        # Process all .dlg files in directory
        for filename in os.listdir(dlg_path):
            if filename.endswith('.dlg'):
                input_file = os.path.join(dlg_path, filename)
                output_file = os.path.join(output_path, f"{os.path.splitext(filename)[0]}_extracted.txt")
                process_single_file(input_file, output_file)
    else:
        # Process single file
        process_single_file(dlg_path, output_path)

def process_single_file(dlg_path, output_path, specified_run=None):
    # Read .dlg file
    with open(dlg_path, 'r') as f:
        lines = f.readlines()

    # Find cluster run with most members, or use specified run
    run_to_extract = specified_run
    if run_to_extract is None:
        run_to_extract = parse_cluster_table(lines)
        
    if run_to_extract is None:
        print(f'Warning: could not find cluster with the most members in {dlg_path}')
        return

    # Extract the docked section for this run
    docked_section = extract_run_section(lines, run_to_extract)
    if not docked_section:
        print(f'Warning: no DOCKED section found for run {run_to_extract} in {dlg_path}')
        return

    # Remove "DOCKED: " prefix from all lines
    cleaned_section = [line.replace("DOCKED: ", "") for line in docked_section]

    # Write to output file with run number and .pdbqt extension
    base_name = os.path.splitext(output_path)[0]
    output_with_run = f"{base_name}_{run_to_extract}.pdbqt"
    with open(output_with_run, 'w') as out:
        out.writelines(cleaned_section)

    # print(f'Extracted DOCKED section for run {run_to_extract} to {output_with_run}')
    return output_with_run

def extract_pdbqt_for_actives_list(active_compounds_list, pdbqt_dir, output_dir):
    """
    Extract PDBQT files for active compounds from a simple list.
    Specifically extracts the correct conformation from DLG files based on the run number in compound name.
    
    Args:
        active_compounds_list: List of active compound names (format: Base_Name_ClusterID_RunID)
        pdbqt_dir: Directory containing DLG files
        output_dir: Directory to save extracted PDBQT files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Extracting PDBQT files for {len(active_compounds_list)} active compounds...")
    
    # Track success and failures
    processed = 0
    successful = 0
    missing = 0
    
    # List all DLG files in the directory
    try:
        all_files = os.listdir(pdbqt_dir)
        dlg_files = [f for f in all_files if f.endswith('.dlg')]
        # print(f"Found {len(dlg_files)} DLG files in {pdbqt_dir}")
    except Exception as e:
        print(f"Error listing directory {pdbqt_dir}: {e}")
        dlg_files = []
    
    # Extract PDBQT files for each active compound using DLG extraction
    for compound in active_compounds_list:
        processed += 1
        
        # Parse the compound name to extract base name and run number
        # Format expected: Base_Name_ClusterID_RunID (e.g., Test_active_91_10_538)
        if '_' in compound:
            parts = compound.split('_')
            if len(parts) >= 4:  # At least 4 parts needed
                # Extract the base name (everything up to the cluster ID) and the run number (last part)
                base_name = '_'.join(parts[:-2])  # "Test_active_91"
                try:
                    run_number = int(parts[-1])   # "538" (as integer)
                    
                    # Look for matching DLG file
                    dlg_file_name = f"{base_name}.dlg"
                    dlg_path = os.path.join(pdbqt_dir, dlg_file_name)
                    
                    if os.path.exists(dlg_path):
                        # print(f"Found DLG file for {compound}: {dlg_file_name}, extracting run {run_number}")
                        
                        # Create temporary output path for extracted PDBQT
                        temp_output = os.path.join(output_dir, f"{compound}_temp")
                        try:
                            # Extract the specified run
                            output_file = process_single_file(dlg_path, temp_output, specified_run=run_number)
                            
                            if output_file and os.path.exists(output_file):
                                # Rename to final name if needed
                                final_file = os.path.join(output_dir, f"{compound}.pdbqt")
                                if output_file != final_file:
                                    os.rename(output_file, final_file)
                                
                                successful += 1
                                # print(f"Successfully extracted run {run_number} from {dlg_file_name}")
                            else:
                                # print(f"Error: Failed to extract run {run_number} from {dlg_file_name}")
                                missing += 1
                        except Exception as e:
                            print(f"Error extracting run {run_number} from {dlg_file_name}: {e}")
                            missing += 1
                    else:
                        print(f"Error: DLG file {dlg_file_name} not found for {compound}")
                        missing += 1
                except ValueError:
                    print(f"Error: Could not parse run number from {compound}")
                    missing += 1
            else:
                print(f"Error: Compound name {compound} does not match expected format (Base_Name_ClusterID_RunID)")
                missing += 1
        else:
            print(f"Error: Compound name {compound} does not contain underscores")
            missing += 1
    
    # print(f"\nPDBQT extraction complete:")
    print(f"  Total compounds processed: {processed}")
    print(f"  Successfully extracted: {successful}")
    print(f"  Failed to extract: {missing}")
        
    return successful > 0


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Extract the DOCKED section for the cluster with the most members from an AutoDock .dlg file.'
    )
    parser.add_argument('--path_dlg_file', help='Path to the .dlg file')
    parser.add_argument('-o', '--out', default='extracted_run.txt', help='Output file for the extracted section')
    args = parser.parse_args()
    main(args.path_dlg_file, args.out)
