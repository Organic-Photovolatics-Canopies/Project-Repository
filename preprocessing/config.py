"""
Configuration constants for CEPDB preprocessing pipeline
"""

from pathlib import Path

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "Datasets" / "cepdb_truncated_first_100_querries"
OUTPUT_DIR = PROJECT_ROOT / "preprocessed_data"

# Preprocessing decisions (from recommended options)
QC_METHOD = "calibrated_prelim"  # Use calibrated values only
INCLUDE_EXPLICIT_H = False       # Implicit hydrogens only
SPLIT_METHOD = "scaffold"        # Scaffold split for generalization
SPLIT_RATIOS = {"train": 0.7, "val": 0.15, "test": 0.15}
SCALING_METHOD = "standard"      # StandardScaler
HANDLE_MISSING = "drop"          # Drop molecules without targets
RANDOM_SEED = 42

# Multi-task regression targets
TARGET_COLUMNS = ["pce", "voc", "jsc"]

# Column schemas (from db_setup.sql)
MOLGRAPH_COLS = [
    'id', 'SMILES_str', 'iupac_str', 'inchi_str', 'cas_str',
    'trivial_namestr', 'stoich_str', 'n_el', 'n_heavyatoms',
    'n_bf_sz', 'n_bf_dzp', 'n_bf_tzp', 'mass', 'permission_level'
]

CALCQC_COLS = [
    'id', 'mol_graph_id', 'mol_geom_id', 'calc_id_str',
    'calc_tbz_str', 'calc_archive_subdir_path', 'modelchem_str',
    'e_total', 'e_homo_alpha', 'e_lumo_alpha', 'e_gap_alpha',
    'e_homo_beta', 'e_lumo_beta', 'e_gap_beta', 'e_gap_min',
    'dipmom_total', 's2_val'
]

CALIBQC_COLS = [
    'id', 'mol_graph_id', 'mol_geom_id', 'calc_qcset1_id',
    'calib_type', 'e_homo_alpha', 'e_lumo_alpha', 'e_gap_alpha',
    'e_homo_beta', 'e_lumo_beta', 'e_gap_beta', 'e_gap_min'
]

SCHARBER_COLS = [
    'id', 'mol_graph_id', 'mol_geom_id', 'calib_id',
    'scharber_type', 'pce', 'voc', 'jsc'
]

# Feature groups
QUANTUM_FEATURES = ['e_homo_alpha', 'e_lumo_alpha', 'e_gap_alpha', 'e_gap_min']
MOLECULAR_FEATURES = ['mass', 'n_heavyatoms', 'n_el']
DERIVED_FEATURES = [
    'heteroatom_ratio', 'electron_affinity_proxy',
    'ionization_potential_proxy', 'aromatic_rings',
    'rotatable_bonds', 'conjugation_estimate'
]
