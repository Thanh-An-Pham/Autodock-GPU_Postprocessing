# Autodock-GPU_Postprocessing
Automated large-scale post-processing for AutoDock-GPU Molecular Docking

## Installation
Tested with Python 3.12.5
Complete dependency list for all tools in the requirements.txt

```bash
# Install dependencies
pip install -r requirements.txt

# Optional: Install PyMOL for PDBQT to PDB conversion
conda install -c conda-forge pymol-open-source
```

## How to Run

The project consists of several Python scripts that can be run directly:

### 1. Protein-Ligand Interaction Analysis
```bash
# Basic interaction analysis
python interaction_analysis.py --protein protein.pdb --ligands ligand_directory/

# Single ligand analysis with 2D mode
python interaction_analysis.py --protein protein.pdb --ligands single_ligand.pdb --mode 2d

# Use configuration file
python interaction_analysis.py --config my_config.yaml
```

### 2. AutoDock Pipeline Processing
```bash
# Complete AutoDock pipeline
python main.py \
  --train_dlg_dir train_data/ \
  --test_dlg_dir test_data/ \
  --output_dir results/

# Show all available options
python main.py --help
```

### 3. Individual Tools
```bash
# PDBQT screening
python screening_pdbqt.py --input_dir pdbqt_files/ --output_dir screened/

# Convert PDBQT to PDB
python convert_pdbqt_to_pdb.py --input file.pdbqt --output file.pdb

# Process datasets
python process_datasets.py --train_dlg_dir Data/Train/ --test_dlg_dir Data/Test/
```

## Basic Usage

### ProLIF Analysis
```bash
# Analyze protein-ligand interactions
python interaction_analysis.py --protein protein.pdb --ligands ligand_directory/

# Single ligand analysis
python interaction_analysis.py --protein protein.pdb --ligands single_ligand.pdb --mode 2d

# Use configuration file
python interaction_analysis.py --config my_config.yaml
```

### AutoDock Pipeline
```bash
# Complete pipeline
python main.py \
  --train_dlg_dir train_data/ \
  --test_dlg_dir test_data/ \
  --output_dir results/
```

## 📖 Documentation

### ProLIF Analysis Tool

The `python interaction_analysis.py` command provides comprehensive protein-ligand interaction analysis:

```bash
# Basic analysis (2D networks + binding sites)
python interaction_analysis.py --protein protein.pdb --ligands ligand_dir/

# Different analysis modes
python interaction_analysis.py --mode 2d --protein protein.pdb --ligands ligands/      # Only 2D networks
python interaction_analysis.py --mode binding --protein protein.pdb --ligands ligands/ # Only binding sites
python interaction_analysis.py --mode parallel --protein protein.pdb --ligands ligands/ # Parallel processing
```
**Outputs**:
- `prolif_results/2d_networks/`: Interactive HTML network visualizations
- `prolif_results/binding_sites/`: CSV data and distance plots
- `prolif_results/analysis_summary.json`: Comprehensive analysis report

### AutoDock Processing Pipeline

The `python main.py` command provides complete AutoDock result analysis:

```bash
python main.py --help  # Show all options

# Basic usage
python main.py \
  --train_dlg_dir /path/to/train/dlg/files \
  --test_dlg_dir /path/to/test/dlg/files \
  --output_dir results/
```

**Pipeline Steps**:
1. **DLG Processing**: Parse AutoDock files and extract binding energies
2. **Training Analysis**: Calculate ROC-AUC metrics and find optimal cutoffs
3. **Test Analysis**: Apply cutoffs to identify active compounds
4. **Structure Extraction**: Extract PDBQT conformations for active compounds
5. **Format Conversion**: Convert to PDB format for further analysis

**Output Structure**:
```
results/
├── train_metrics/
│   ├── train_combined.csv
│   ├── train_stats.csv
│   └── train_metrics/
│       ├── metrics.csv
│       ├── roc_curves.png
│       ├── logAUC_curves.png
│       └── enrichment_curves.png
└── test_metrics/
    ├── test_combined.csv
    ├── test_stats.csv
    ├── active_compounds.csv
    └── active_pdbqt/
        ├── compound1.pdb
        ├── compound1.pdbqt
        └── ...
```

## 🔧 Advanced Usage

### Integration Workflow

Combine ProLIF analysis with AutoDock pipeline:

```bash
# 1. Process AutoDock results
python main.py \
  --train_dlg_dir train_dlg/ \
  --test_dlg_dir test_dlg/ \
  --output_dir autodock_results/

# 2. Analyze active compounds with ProLIF
python interaction_analysis.py \
  --protein target_protein.pdb \
  --ligands autodock_results/test_metrics/active_pdbqt/ \
  --output-dir interaction_analysis/

# 3. Results are in:
# - autodock_results/: Machine learning metrics and active compounds
# - interaction_analysis/: Detailed protein-ligand interactions
```
