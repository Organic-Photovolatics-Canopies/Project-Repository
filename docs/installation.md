# Installation Guide

[← Back to Documentation Home](README.md)

This guide will walk you through setting up the OMID USA development environment for preprocessing molecular data and training GNN models for OPV property prediction.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Step-by-Step Installation](#step-by-step-installation)
- [Verification](#verification)
- [GPU Setup (Optional)](#gpu-setup-optional)
- [Docker Setup (Coming Soon)](#docker-setup-coming-soon)
- [Platform-Specific Notes](#platform-specific-notes)

---

## Prerequisites

Before installing, ensure you have the following:

### Required Software

- **Python 3.11**: The project is developed and tested with Python 3.11
- **Conda**: Anaconda or Miniconda for environment management
  - Download: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/download)
- **Git**: For cloning the repository
  - Download: [Git](https://git-scm.com/downloads)

### System Requirements

- **RAM**: Minimum 8 GB, recommended 16 GB for processing larger datasets
- **Storage**: At least 5 GB free space for dependencies and data
- **CPU**: Multi-core processor recommended for faster preprocessing
- **OS**: Linux, macOS, or Windows (see [platform-specific notes](#platform-specific-notes))

---

## Step-by-Step Installation

### 1. Clone the Repository

```bash
# Clone the repository
git clone https://github.com/your-org/Project-Repository.git

# Navigate to the project directory
cd Project-Repository
```

### 2. Create Conda Environment

Create a dedicated conda environment with Python 3.11:

```bash
# Create the environment
conda create -n opv_preprocessing python=3.11 -y

# Activate the environment
conda activate opv_preprocessing
```

> **Note**: You must activate this environment every time you work with the project.

### 3. Install Core Dependencies

Install PyTorch, PyTorch Geometric, RDKit, and other essential packages:

```bash
# Install from conda-forge and pytorch channels
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

> **For GPU users**: See [GPU Setup](#gpu-setup-optional) for CUDA-enabled installation.

### 4. Install Additional Python Packages

Install any additional dependencies via pip:

```bash
# Install from requirements if available
pip install -r requirements.txt  # If file exists

# Or install individual packages
pip install networkx tqdm
```

### 5. Verify RDKit Installation

RDKit is critical for molecular processing. Verify it's installed correctly:

```bash
python -c "from rdkit import Chem; print('RDKit version:', Chem.rdBase.rdkitVersion)"
```

Expected output: `RDKit version: 2023.x.x` or similar

### 6. Verify PyTorch Geometric

Check PyTorch Geometric installation:

```bash
python -c "import torch_geometric; print('PyG version:', torch_geometric.__version__)"
```

Expected output: `PyG version: 2.x.x` or similar

### 7. Download Dataset

Clone the OPV2D dataset:

```bash
# Clone OPV2D repository
git clone https://github.com/sunyrain/OPV2D.git

# Verify data files
ls OPV2D/data/
```

**Note**: Preprocessing scripts for OPV2D will be added to the repository soon.

---

## Verification

### Test the Installation

Run a quick test to ensure everything is set up correctly:

```bash
# Activate the environment (if not already active)
conda activate opv_preprocessing

# Test imports
python -c "
from rdkit import Chem
import torch
import torch_geometric
from sklearn.preprocessing import StandardScaler
import pandas as pd
print('✅ All core dependencies imported successfully!')
"
```

### Run Preprocessing Pipeline

Test the full preprocessing pipeline on the sample dataset:

```bash
# From the project root directory
python -m preprocessing.pipeline
```

Expected behavior:
- Pipeline initializes and begins processing
- Progress messages appear for each step (1-10)
- Creates `preprocessed_data/` directory with processed outputs
- Displays summary statistics at completion

**Estimated time**: Varies by dataset size (see user guide)

### Check Output Files

Verify the preprocessing created the expected output structure:

```bash
ls -R preprocessed_data/
```

Expected structure:
```
preprocessed_data/
├── molecules.csv          # Master feature table
├── graphs/                # PyG graph objects
│   ├── mol_0.pt
│   ├── mol_1.pt
│   └── ...
├── graphormer_inputs/     # Graphormer-encoded graphs
│   ├── mol_0.pt
│   └── ...
├── scalers/               # Feature and target scalers
│   ├── feature_scaler.pkl
│   └── target_scaler.pkl
├── splits/                # Train/val/test indices
│   ├── train_indices.npy
│   ├── val_indices.npy
│   └── test_indices.npy
└── metadata.json          # Pipeline configuration
```

---

## GPU Setup (Optional)

For faster model training, install PyTorch with CUDA support:

### Check CUDA Availability

```bash
# Check if NVIDIA GPU is available
nvidia-smi
```

### Install PyTorch with CUDA

Replace the CPU-only PyTorch installation:

```bash
# Remove CPU-only PyTorch
conda remove pytorch torchvision cpuonly

# Install CUDA-enabled PyTorch (for CUDA 11.8)
conda install -c pytorch -c nvidia \
  pytorch torchvision pytorch-cuda=11.8 \
  -y

# For CUDA 12.1
conda install -c pytorch -c nvidia \
  pytorch torchvision pytorch-cuda=12.1 \
  -y
```

### Verify GPU Setup

```bash
python -c "
import torch
print('CUDA available:', torch.cuda.is_available())
print('CUDA version:', torch.version.cuda)
print('GPU device:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'None')
"
```

> **Note**: CUDA setup is optional. The preprocessing pipeline works efficiently on CPU.

---

## Docker Setup (Coming Soon)

A Docker containerized version for easy deployment is under development.

```bash
# Future Docker commands (placeholder)
docker pull omidusa/opv-preprocessing:latest
docker run -v $(pwd)/data:/app/data omidusa/opv-preprocessing
```

---

## Platform-Specific Notes

### Linux

The recommended platform. All instructions above should work without modification.

**Ubuntu/Debian specific**:
```bash
# If you encounter build issues
sudo apt-get update
sudo apt-get install build-essential python3-dev
```

### macOS

Works well on both Intel and Apple Silicon (M1/M2/M3) Macs.

**Apple Silicon (M1/M2/M3)**:
```bash
# Use the native ARM64 conda
# Download Miniforge instead of Miniconda
# https://github.com/conda-forge/miniforge

# Installation is otherwise identical
```

### Windows

Supported via Anaconda Prompt or Windows Subsystem for Linux (WSL).

**Option 1: Anaconda Prompt (Native Windows)**
- Use Anaconda Prompt (not PowerShell or CMD)
- All conda commands work the same
- RDKit installation can be slower

**Option 2: WSL (Recommended)**
```bash
# Install WSL2 with Ubuntu
wsl --install

# Then follow Linux instructions inside WSL
```

---

## Troubleshooting Installation

### RDKit Installation Fails

**Issue**: `PackagesNotFoundError: rdkit`

**Solution**:
```bash
# Ensure conda-forge channel is added
conda config --add channels conda-forge
conda config --set channel_priority strict

# Reinstall RDKit
conda install -c conda-forge rdkit -y
```

### PyTorch Geometric Installation Issues

**Issue**: Import errors or missing dependencies

**Solution**:
```bash
# Install PyG with explicit version
conda install pytorch-geometric=2.3.0 -c pytorch -c conda-forge -y

# Or install from PyG directly
pip install torch-geometric
```

### Conda Environment Conflicts

**Issue**: Dependency conflicts during installation

**Solution**:
```bash
# Remove the environment and start fresh
conda deactivate
conda env remove -n opv_preprocessing

# Recreate with explicit Python version
conda create -n opv_preprocessing python=3.11 -y
conda activate opv_preprocessing

# Install packages one by one if needed
conda install -c conda-forge rdkit -y
conda install -c pytorch pytorch cpuonly -y
# ... continue individually
```

### Missing `preprocessed_data/` Directory

**Issue**: Pipeline runs but no output directory created

**Solution**:
```bash
# Check you're in the project root directory
pwd  # Should show: .../Project-Repository

# Manually create the directory if needed
mkdir -p preprocessed_data

# Run pipeline again
python -m preprocessing.pipeline
```

### Memory Errors During Preprocessing

**Issue**: `MemoryError` or system slowdown

**Solution**:
- Close other applications to free RAM
- Process smaller batches by modifying `preprocessing/config.py`
- Use a machine with more RAM
- Consider cloud computing options (AWS, Google Colab)

---

## Environment Variables (Optional)

Set optional environment variables for customization:

```bash
# Set data directory path (if not using default)
export OPV_DATA_PATH="/path/to/your/data"

# Set output directory
export OPV_OUTPUT_PATH="/path/to/output"

# Add to your ~/.bashrc or ~/.zshrc to persist
echo 'export OPV_DATA_PATH="/path/to/your/data"' >> ~/.bashrc
```

---

## Next Steps

✅ **Installation Complete!**

- Proceed to the [User Guide](user-guide.md) for a quick start tutorial
- Read the [User Manual](user-manual.md) for detailed technical information
- Check the [FAQ](faq.md) if you have questions

---

## Getting Help

If you encounter issues not covered here:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Review the [FAQ](faq.md)
3. Search existing issues on GitHub
4. Contact the development team (see [team members](README.md#team-members))

---

[← Back to Documentation Home](README.md)
