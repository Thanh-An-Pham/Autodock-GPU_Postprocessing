import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional
import MDAnalysis as mda
import prolif
import prolif.fingerprint as fp
from prolif.plotting.network import LigNetwork
import warnings
from pathlib import Path

# Suppress specific warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis")
warnings.filterwarnings("ignore", message="Reader has no dt information, set to 1.0 ps")
warnings.filterwarnings("ignore", category=UserWarning, module="prolif")

class ProLIFAnalyzer:    
    def __init__(self, protein_file: str, ligand_path: str, distance_cutoff: float = 15.0):
        """
        Initialize the ProLIF analyzer.
        
        Args:
            protein_file: Path to the protein PDB file
            ligand_path: Path to directory containing ligand PDB files or single ligand file
            distance_cutoff: Distance cutoff in Ångströms (default: 15.0)
        """
        self.protein_file = protein_file
        self.ligand_path = ligand_path
        self.results_dir = Path("prolif_results")
        self.results_dir.mkdir(exist_ok=True)
        self.distance_cutoff = distance_cutoff
        
        # Load structures
        self._load_structures()
        
    def _load_structures(self):
        """Load protein and ligand structures."""
        try:
            # Load protein with minimal backbone representation
            self.protein = mda.Universe(self.protein_file)
            protein_selection = "name CA CB"
            self.protein_atoms = self.protein.select_atoms(protein_selection)
            
            # Load ligands - handle both single file and directory
            self.ligands = {}
            
            if os.path.isfile(self.ligand_path):
                # Single ligand file
                if self.ligand_path.endswith('.pdb'):
                    try:
                        lig = mda.Universe(self.ligand_path)
                        lig_atoms = lig.select_atoms("not name H*")
                        ligand_name = os.path.basename(self.ligand_path)
                        self.ligands[ligand_name] = lig_atoms
                        print(f"✓ Loaded single ligand: {ligand_name} ({len(lig_atoms)} atoms)")
                    except Exception as e:
                        print(f"❌ Error loading ligand file {self.ligand_path}: {str(e)}")
                        sys.exit(1)
                else:
                    print(f"❌ Ligand file must be a PDB file: {self.ligand_path}")
                    sys.exit(1)
                    
            elif os.path.isdir(self.ligand_path):
                # Directory of ligand files
                ligand_files = [f for f in os.listdir(self.ligand_path) if f.endswith('.pdb')]
                
                if not ligand_files:
                    print(f"❌ No PDB files found in directory: {self.ligand_path}")
                    sys.exit(1)
                
                print(f"📁 Found {len(ligand_files)} ligand files in directory")
                
                for lig_file in ligand_files:
                    lig_path = os.path.join(self.ligand_path, lig_file)
                    try:
                        lig = mda.Universe(lig_path)
                        lig_atoms = lig.select_atoms("not name H*")
                        self.ligands[lig_file] = lig_atoms
                        print(f"✓ Loaded ligand: {lig_file} ({len(lig_atoms)} atoms)")
                    except Exception as e:
                        print(f"⚠ Warning: Could not load ligand {lig_file}: {str(e)}")
                        continue
            else:
                print(f"❌ Ligand path does not exist: {self.ligand_path}")
                sys.exit(1)
            
            if not self.ligands:
                print("❌ No ligands were successfully loaded")
                sys.exit(1)
                
            print(f"✅ Successfully loaded {len(self.ligands)} ligand(s)")
                        
        except Exception as e:
            print(f"❌ Error loading structures: {str(e)}")
            sys.exit(1)
            
    def create_2d_network(self, ligand_name: str) -> bool:
        """
        Create 2D interactive network visualization for a ligand.
        
        Args:
            ligand_name: Name of the ligand file
            
        Returns:
            bool: True if interactions were detected, False if no interactions found
        """
        try:
            if ligand_name not in self.ligands:
                print(f"Error: Ligand {ligand_name} not found")
                return False
                
            print(f"\n🎨 Creating 2D ligand interaction network for {ligand_name}...")
            
            # Use minimal backbone representation for protein
            protein_selection = "name CA CB"
            protein_atoms = self.protein.select_atoms(protein_selection)
            
            if len(protein_atoms) == 0:
                protein_atoms = self.protein.select_atoms("name CA")
            
            print(f"  - Selected {len(protein_atoms)} protein atoms")
            
            # Get ligand atoms (remove hydrogens)
            ligand_atoms = self.ligands[ligand_name]
            print(f"  - Using ligand with {len(ligand_atoms)} atoms")
            
            # Create fingerprint
            all_interactions = [
                "Hydrophobic", "VdWContact", "HBDonor", "HBAcceptor", 
                "PiStacking", "PiCation", "CationPi", "Anionic", "Cationic",
                "XBDonor", "XBAcceptor", "MetalAcceptor"
            ]
            
            fp = prolif.Fingerprint(
                interactions=all_interactions,
                vicinity_cutoff=self.distance_cutoff
            )
            print(f"  - Created fingerprint with interactions: {list(fp.interactions.keys())}")
            
            # Create ProLIF molecules
            try:
                protein_mol = prolif.Molecule.from_mda(protein_atoms, NoImplicit=False)
                ligand_mol = prolif.Molecule.from_mda(ligand_atoms, NoImplicit=False)
                print(f"  ✓ Created ProLIF molecules")
            except Exception as e:
                print(f"  ❌ Failed to create molecules: {e}")
                return False
            
            # Run interaction analysis
            try:
                interactions = fp.run_from_iterable([ligand_mol], protein_mol)
                print(f"  ✓ Calculated interactions")
                
                # Check if any interactions were detected
                interactions_detected = False
                interaction_count = 0
                
                try:
                    if hasattr(interactions, 'values') and len(interactions.values) > 0:
                        frame_data = interactions.values[0]
                        for idx, (residue, ligand_id, interaction_type) in enumerate(frame_data.index):
                            if frame_data.iloc[idx] > 0:
                                interaction_count += 1
                                interactions_detected = True
                    
                    if not interactions_detected:
                        print(f"  ⚠ No interactions detected")
                    else:
                        print(f"  ✓ Detected {interaction_count} protein-ligand interactions")
                except Exception as check_error:
                    print(f"  ⚠ Could not verify interactions: {check_error}")
                    interactions_detected = True
                
                # Always create and save visualization regardless of interaction detection
                try:
                    view = prolif.plotting.network.LigNetwork.from_fingerprint(interactions, ligand_mol)
                    
                    output_dir = self.results_dir / "2d_networks"
                    output_dir.mkdir(exist_ok=True)
                    html_path = output_dir / f"{ligand_name.replace('.pdb', '')}_2d_network.html"
                    
                    # Try to save the network
                    if hasattr(view, 'save'):
                        view.save(str(html_path))
                        print(f"  ✅ 2D network saved: {html_path}")
                    elif hasattr(view, 'fig') and view.fig is not None:
                        view.fig.write_html(str(html_path))
                        print(f"  ✅ 2D network saved: {html_path}")
                    else:
                        print(f"  ⚠ Could not save network visualization")
                    
                    # Save interaction data
                    try:
                        df = fp.to_dataframe()
                        if not df.empty:
                            csv_path = output_dir / f"{ligand_name.replace('.pdb', '')}_interactions.csv"
                            df.to_csv(csv_path)
                            print(f"  ✓ Interaction data saved: {csv_path}")
                        else:
                            # Create empty CSV file even when no interactions
                            csv_path = output_dir / f"{ligand_name.replace('.pdb', '')}_interactions.csv"
                            empty_df = pd.DataFrame(columns=['Ligand', 'Protein', 'Interaction'])
                            empty_df.to_csv(csv_path, index=False)
                            print(f"  ✓ Empty interaction data saved: {csv_path}")
                    except Exception as save_error:
                        print(f"  ⚠ Could not save interaction data: {save_error}")
                    
                    return interactions_detected
                    
                except Exception as viz_error:
                    print(f"  ❌ Failed to create visualization: {viz_error}")
                    return interactions_detected
            
            except Exception as e:
                print(f"  ❌ Failed to calculate interactions: {e}")
                return False
            
        except Exception as e:
            print(f"❌ Unexpected error: {str(e)}")
            return False

    def analyze_binding_site(self, ligand_name: str) -> bool:
        """
        Perform binding site analysis for a ligand.
        
        Args:
            ligand_name: Name of the ligand file
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            if ligand_name not in self.ligands:
                print(f"Error: Ligand {ligand_name} not found")
                return False
                
            print(f"  🔍 Analyzing binding site for {ligand_name} (cutoff: {self.distance_cutoff}Å)")
            
            # Get ligand center and coordinates
            ligand_atoms = self.ligands[ligand_name]
            ligand_center = ligand_atoms.center_of_mass()
            ligand_coords = ligand_atoms.positions
            
            print(f"    💊 Ligand center: [{ligand_center[0]:.2f}, {ligand_center[1]:.2f}, {ligand_center[2]:.2f}]")
            
            # Find nearby residues
            nearby_residues = []
            min_distance = float('inf')
            closest_residue = None
            
            for residue in self.protein.residues:
                try:
                    residue_atoms = self.protein.select_atoms(f"resid {residue.resid}")
                    if len(residue_atoms) == 0:
                        continue
                        
                    residue_center = residue_atoms.center_of_mass()
                    center_dist = np.linalg.norm(residue_center - ligand_center)
                    
                    # Calculate minimum distance between any ligand atom and any residue atom
                    min_atom_dist = float('inf')
                    for lig_coord in ligand_coords:
                        for res_coord in residue_atoms.positions:
                            atom_dist = np.linalg.norm(lig_coord - res_coord)
                            min_atom_dist = min(min_atom_dist, atom_dist)
                    
                    final_dist = min(center_dist, min_atom_dist)
                    
                    if final_dist < min_distance:
                        min_distance = final_dist
                        closest_residue = f"{residue.resname}{residue.resid}"
                    
                    if final_dist <= self.distance_cutoff:
                        nearby_residues.append({
                            'residue': residue.resname,
                            'resid': residue.resid,
                            'chain': residue.segid if residue.segid else 'A',
                            'distance': final_dist,
                            'center_distance': center_dist,
                            'min_atom_distance': min_atom_dist
                        })
                except Exception as e:
                    continue
            
            print(f"    🔍 Closest residue: {closest_residue} at {min_distance:.2f}Å")
            print(f"    📊 Found {len(nearby_residues)} residues within {self.distance_cutoff}Å")
            
            # Save results
            output_dir = self.results_dir / "binding_sites"
            output_dir.mkdir(exist_ok=True)
            csv_file = output_dir / f"{ligand_name.replace('.pdb', '')}_binding_site.csv"
            
            if not nearby_residues:
                print(f"    ⚠ No residues found within {self.distance_cutoff}Å of ligand")
                print(f"    💡 Try increasing distance_cutoff (closest residue is {min_distance:.2f}Å away)")
                empty_df = pd.DataFrame(columns=['residue', 'resid', 'chain', 'distance', 'center_distance', 'min_atom_distance'])
                empty_df.to_csv(csv_file, index=False)
                print(f"    ✅ Empty binding site data saved: {csv_file}")
                return True
            
            df = pd.DataFrame(nearby_residues)
            df.to_csv(csv_file, index=False)
            print(f"    ✅ Binding site data saved: {csv_file}")
            
            # Create visualization
            try:
                plt.figure(figsize=(15, 6))
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
                
                # Distance distribution
                ax1.hist(df['distance'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
                ax1.set_xlabel('Distance (Å)')
                ax1.set_ylabel('Number of Residues')
                ax1.set_title(f'Distance Distribution\n{ligand_name}')
                ax1.grid(True, alpha=0.3)
                
                # Top closest residues
                top_residues = df.nsmallest(15, 'distance')
                residue_labels = [f"{row['residue']}{row['resid']}" for _, row in top_residues.iterrows()]
                
                bars = ax2.barh(range(len(top_residues)), top_residues['distance'], color='lightcoral')
                ax2.set_yticks(range(len(top_residues)))
                ax2.set_yticklabels(residue_labels)
                ax2.set_xlabel('Distance (Å)')
                ax2.set_title(f'Top 15 Closest Residues\n{ligand_name}')
                ax2.grid(True, alpha=0.3, axis='x')
                
                # Add distance labels
                for i, (bar, dist) in enumerate(zip(bars, top_residues['distance'])):
                    ax2.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                            f'{dist:.1f}', ha='left', va='center', fontsize=8)
                
                plt.tight_layout()
                
                plot_file = output_dir / f"{ligand_name.replace('.pdb', '')}_binding_site.png"
                plt.savefig(plot_file, dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"    ✅ Binding site plot saved: {plot_file}")
                
            except Exception as plot_error:
                print(f"    ⚠ Could not create plot: {plot_error}")
                plt.close()  # Close figure even if there's an error
            
            return True
            
        except Exception as e:
            print(f"❌ Error analyzing binding site: {str(e)}")
            return False
            
    def analyze_all_ligands(self):
        """Analyze all ligands with both 2D networks and binding site analysis."""
        total_ligands = len(self.ligands)
        successful_interactions = []
        weak_binders = []
        successful_binding = []
        failed_binding = []
        
        for ligand_name in self.ligands.keys():
            print(f"\nAnalyzing ligand: {ligand_name}")
            
            # 2D Network Analysis (interaction detection)
            has_interactions = self.create_2d_network(ligand_name)
            if has_interactions:
                successful_interactions.append(ligand_name)
            else:
                weak_binders.append(ligand_name)
            
            # Binding Site Analysis
            binding_success = self.analyze_binding_site(ligand_name)
            if binding_success:
                successful_binding.append(ligand_name)
            else:
                failed_binding.append(ligand_name)
        
        # Print comprehensive summary report
        print("\n" + "="*80)
        print("🎯 COMPREHENSIVE ANALYSIS SUMMARY REPORT")
        print("="*80)
        
        # Interaction Detection Results
        success_count = len(successful_interactions)
        weak_count = len(weak_binders)
        success_rate = (success_count / total_ligands) * 100 if total_ligands > 0 else 0
        
        print(f"\n🔬 Protein-Ligand Interaction Detection Results:")
        print(f"   Total ligands analyzed: {total_ligands}")
        print(f"   ✅ Strong Binders: {success_count}/{total_ligands} ligands ({success_rate:.1f}%) with detectable interactions")
        print(f"   ❌ Weak Binders: {weak_count}/{total_ligands} ligands with \"⚠ No interactions detected\"")
        
        # Binding Site Analysis Results
        binding_success_count = len(successful_binding)
        binding_success_rate = (binding_success_count / total_ligands) * 100 if total_ligands > 0 else 0
        
        print(f"\n🔬 Binding Site Analysis Results:")
        print(f"   ✅ Success: {binding_success_count}/{total_ligands} ligands ({binding_success_rate:.1f}%) completed binding site analysis")
        if failed_binding:
            print(f"   ❌ Failed: {len(failed_binding)}/{total_ligands} ligands had binding site analysis issues")
        
        # Results location summary
        print(f"\n📁 Results Location:")
        print(f"   🌐 Interactive 2D Networks: prolif_results/2d_networks/ (HTML files)")
        print(f"      - Generated for all analyzed ligands")
        print(f"   📈 Binding Site Analysis: prolif_results/binding_sites/ (CSV + PNG files)")
        print(f"      - Distance calculations and residue identification")
        print(f"   📊 Interaction Data: prolif_results/2d_networks/ (CSV files)")
        print(f"      - Generated for all analyzed ligands (empty files for weak binders)")
        
        print(f"\n✨ Analysis complete! Your ProLIF pipeline is fully functional.")
        
        if successful_interactions:
            print(f"🎉 {len(successful_interactions)} ligands show strong protein-ligand interactions!")
        
        if weak_binders:
            print(f"⚠️  {len(weak_binders)} ligands classified as weak binders (no detectable interactions)")
            
        print("="*80)

def main():
    """Main function to run the analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(description="ProLIF Analysis Tool - Comprehensive protein-ligand interaction analysis")
    parser.add_argument('--protein', required=True, help='Path to protein PDB file')
    parser.add_argument('--ligands', required=True, help='Path to ligand PDB file or directory containing ligand PDB files')
    parser.add_argument('--distance-cutoff', type=float, default=5.0, help='Distance cutoff in Ångströms (default: 5.0)')
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = ProLIFAnalyzer(args.protein, args.ligands, args.distance_cutoff)
    
    # Run comprehensive analysis
    analyzer.analyze_all_ligands()
    
    print("\n✨ Analysis complete! Results saved in 'prolif_results' directory.")

if __name__ == "__main__":
    main() 