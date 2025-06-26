import os
import sys
import pymol
pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet mode

def convert_pdbqt_to_pdb(input_path):
    """
    Convert all .pdbqt files in the given path to .pdb files using PyMOL
    """
    # Check if path exists
    if not os.path.exists(input_path):
        print(f"Error: Path {input_path} does not exist")
        return

    # Process all .pdbqt files in directory
    for filename in os.listdir(input_path):
        if filename.endswith('.pdbqt'):
            input_file = os.path.join(input_path, filename)
            output_file = os.path.join(input_path, filename.replace('.pdbqt', '.pdb'))
            
            # Load the pdbqt file
            pymol.cmd.load(input_file, 'temp')
            
            # Save as pdb
            pymol.cmd.save(output_file, 'temp')
            
            # Delete the temporary object
            pymol.cmd.delete('temp')
            
            # print(f"Converted {filename} to {os.path.basename(output_file)}")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Convert .pdbqt files to .pdb files using PyMOL'
    )
    parser.add_argument('--path', required=True, help='Path containing .pdbqt files')
    args = parser.parse_args()
    
    convert_pdbqt_to_pdb(args.path)
    
    # Quit PyMOL
    pymol.cmd.quit() 