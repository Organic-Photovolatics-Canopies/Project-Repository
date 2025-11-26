# OPV Preprocessing Pipeline - Setup Guide

## Quick Start (Recommended for Next Session)

### Step 1: Create a Fresh Conda Environment

```bash
# Create new environment with Python 3.11
conda create -n opv_preprocessing python=3.11 -y

# Activate the environment
conda activate opv_preprocessing
```

### Step 2: Install All Dependencies (One Command)

```bash
# Install everything via conda for maximum compatibility
conda install -c conda-forge -c pytorch \
  rdkit \
  pytorch \
  torchvision \
  cpuonly \
  pytorch-geometric \
  scikit-learn \
  pandas \
  matplotlib \
  seaborn \
  jupyter \
  -y
```

**Note:** If the above command fails, install in two steps:

```bash
# Step 2a: Core ML packages from PyTorch channel
conda install -c pytorch pytorch torchvision cpuonly -y

# Step 2b: Chemistry and data science packages
conda install -c conda-forge rdkit pytorch-geometric scikit-learn pandas matplotlib seaborn jupyter -y
```

### Step 3: Run the Preprocessing Pipeline

```bash
# Navigate to repository
cd C:\Users\nhamd\dev_space\Capstone_OPV\repo\claude-repo

# Run the pipeline
python -m preprocessing.pipeline
```

### Step 4: Verify Output

Check that the following directory structure was created:

```
preprocessed_data/
├── molecules.csv                    # Master feature table (~186 molecules)
├── graphs/                          # PyG graph objects (mol_0.pt, mol_1.pt, ...)
├── graphormer_inputs/               # Graphormer-encoded graphs
├── scalers/
│   ├── feature_scaler.pkl          # StandardScaler for features
│   └── target_scaler.pkl           # StandardScaler for targets
├── splits/
│   ├── train_indices.npy           # Training set indices (~130 molecules)
│   ├── val_indices.npy             # Validation set indices (~28 molecules)
│   └── test_indices.npy            # Test set indices (~28 molecules)
└── metadata.json                    # Preprocessing configuration
```

---

## What the Pipeline Does

1. **Data Loading**: Loads CEPDB CSV files from `Datasets/cepdb_truncated_first_100_querries/`
2. **Integration**: Merges tables and selects calibrated QC values
3. **Graph Building**: Converts SMILES to PyTorch Geometric graphs
4. **Feature Engineering**: Creates derived molecular features
5. **Target Preparation**: Prepares multi-task targets (PCE, Voc, Jsc)
6. **Validation**: Validates SMILES, removes duplicates, detects outliers
7. **Splitting**: Scaffold-based train/val/test split (70/15/15)
8. **Scaling**: Fits scalers on training data only (prevents data leakage)
9. **Graphormer Encoding**: Adds spatial encoding, centrality, attention bias
10. **Saving**: Saves all preprocessed data and metadata
