# Troubleshooting Guide

[← Back to Documentation Home](README.md)

This guide helps you diagnose and resolve common issues with the OMID USA preprocessing pipeline and model training.

## Table of Contents

- [Installation Issues](#installation-issues)
- [Data Loading Errors](#data-loading-errors)
- [Preprocessing Failures](#preprocessing-failures)
- [Memory Issues](#memory-issues)
- [Model Training Problems](#model-training-problems)
- [Output Interpretation](#output-interpretation)
- [Platform-Specific Issues](#platform-specific-issues)

---

## Installation Issues

### Problem: RDKit Installation Fails

**Error Messages**:
```
PackagesNotFoundError: The following packages are not available from current channels:
  - rdkit
```

**Solutions**:

**Solution 1: Add conda-forge channel**
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install rdkit -y
```

**Solution 2: Explicit channel in install command**
```bash
conda install -c conda-forge rdkit -y
```

**Solution 3: Fresh environment**
```bash
# Remove old environment
conda env remove -n opv_preprocessing

# Create fresh environment with conda-forge default
conda create -n opv_preprocessing -c conda-forge python=3.11 -y
conda activate opv_preprocessing
conda install -c conda-forge rdkit -y
```

**Solution 4: Use Mamba (faster alternative to conda)**
```bash
conda install mamba -c conda-forge
mamba install rdkit -c conda-forge -y
```

**Verification**:
```bash
python -c "from rdkit import Chem; print('RDKit version:', Chem.rdBase.rdkitVersion)"
```

---

### Problem: PyTorch Geometric Installation Issues

**Error Messages**:
```
ImportError: cannot import name 'Data' from 'torch_geometric.data'
```

**Solutions**:

**Solution 1: Install from conda**
```bash
conda install pytorch-geometric -c pytorch -c conda-forge -y
```

**Solution 2: Install from PyG directly**
```bash
pip install torch-geometric
```

**Solution 3: Manual installation with dependencies**
```bash
# Install PyTorch first
conda install pytorch cpuonly -c pytorch -y

# Install PyG and dependencies
pip install torch-geometric
pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
```

**Solution 4: Version-specific installation**
```bash
conda install pytorch-geometric=2.3.0 -c pytorch -y
```

**Verification**:
```bash
python -c "import torch_geometric; print('PyG version:', torch_geometric.__version__)"
```

---

### Problem: Conda Environment Conflicts

**Error Messages**:
```
UnsatisfiableError: The following specifications were found to be incompatible with each other
```

**Solutions**:

**Solution 1: Create minimal environment first**
```bash
# Start with Python only
conda create -n opv_preprocessing python=3.11 -y
conda activate opv_preprocessing

# Install packages one by one
conda install -c conda-forge rdkit -y
conda install -c pytorch pytorch cpuonly -y
conda install -c conda-forge scikit-learn pandas -y
conda install pytorch-geometric -c pytorch -c conda-forge -y
```

**Solution 2: Use explicit package list**
```bash
conda create -n opv_preprocessing python=3.11 rdkit pytorch cpuonly pytorch-geometric scikit-learn pandas matplotlib seaborn jupyter -c conda-forge -c pytorch -y
```

**Solution 3: Relax strict channel priority**
```bash
conda config --set channel_priority flexible
```

**Solution 4: Use environment.yml file**
```yaml
# Create environment.yml
name: opv_preprocessing
channels:
  - conda-forge
  - pytorch
  - defaults
dependencies:
  - python=3.11
  - rdkit
  - pytorch
  - cpuonly
  - pytorch-geometric
  - scikit-learn
  - pandas
  - matplotlib
  - seaborn
  - jupyter
```

```bash
conda env create -f environment.yml
```

---

### Problem: ImportError for installed packages

**Error Messages**:
```
ModuleNotFoundError: No module named 'rdkit'
```

**But conda shows it's installed**:
```bash
conda list | grep rdkit
# rdkit  2023.3.2  py311h...
```

**Solutions**:

**Solution 1: Verify correct environment is active**
```bash
conda env list  # Check active environment (has *)
conda activate opv_preprocessing
```

**Solution 2: Check Python interpreter**
```bash
which python
# Should be: /path/to/conda/envs/opv_preprocessing/bin/python

python --version
# Should be: Python 3.11.x
```

**Solution 3: Reinstall in current environment**
```bash
conda activate opv_preprocessing
conda install --force-reinstall rdkit -c conda-forge
```

**Solution 4: Use full Python path**
```bash
/path/to/conda/envs/opv_preprocessing/bin/python -m preprocessing.pipeline
```

---

## Data Loading Errors

### Problem: FileNotFoundError for dataset

**Error Messages**:
```
FileNotFoundError: [Errno 2] No such file or directory: 'HCEPDB/data_calcqcset1.csv'
```

**Solutions**:

**Solution 1: Check current directory**
```bash
pwd
# Should be: .../Project-Repository

# If not, navigate to project root
cd /path/to/Project-Repository
```

**Solution 2: Verify file exists**
```bash
ls -la HCEPDB/
# Should show: data_calcqcset1.csv
```

**Solution 3: Check file permissions**
```bash
ls -l HCEPDB/data_calcqcset1.csv
# Should be readable: -rw-r--r--

# If not, fix permissions
chmod 644 HCEPDB/data_calcqcset1.csv
```

**Solution 4: Use absolute path**
```python
# Edit preprocessing/config.py
import os
DATA_PATH = os.path.join(os.path.dirname(__file__), '..', 'HCEPDB', 'data_calcqcset1.csv')
DATA_PATH = os.path.abspath(DATA_PATH)
```

---

### Problem: CSV Parsing Errors

**Error Messages**:
```
pandas.errors.ParserError: Error tokenizing data. C error: Expected 20 fields in line 42, saw 21
```

**Solutions**:

**Solution 1: Specify CSV format**
```python
# In data_loader.py
df = pd.read_csv(data_path, sep=',', encoding='utf-8', low_memory=False)
```

**Solution 2: Handle problematic rows**
```python
df = pd.read_csv(data_path, on_bad_lines='skip')  # Skip bad rows
```

**Solution 3: Inspect the problematic line**
```bash
# View line 42
sed -n '42p' HCEPDB/data_calcqcset1.csv
```

**Solution 4: Clean the CSV file**
```python
# Remove problematic characters
import pandas as pd

# Read with error handling
chunks = []
for chunk in pd.read_csv('HCEPDB/data_calcqcset1.csv', chunksize=1000, on_bad_lines='skip'):
    chunks.append(chunk)

df = pd.concat(chunks, ignore_index=True)
df.to_csv('HCEPDB/data_calcqcset1_cleaned.csv', index=False)
```

---

### Problem: Missing Required Columns

**Error Messages**:
```
KeyError: 'smiles' not found in DataFrame
```

**Solutions**:

**Solution 1: Check column names**
```python
import pandas as pd
df = pd.read_csv('HCEPDB/data_calcqcset1.csv')
print(df.columns.tolist())
```

**Solution 2: Rename columns**
```python
# If column is named differently
column_mapping = {
    'Smiles': 'smiles',
    'SMILES': 'smiles',
    'canonical_smiles': 'smiles'
}
df = df.rename(columns=column_mapping)
df.to_csv('HCEPDB/data_calcqcset1.csv', index=False)
```

**Solution 3: Verify required columns exist**
```python
required_columns = ['smiles', 'InChIKey', 'homo', 'lumo', 'pce', 'voc', 'jsc']
missing = [col for col in required_columns if col not in df.columns]
if missing:
    print(f"Missing columns: {missing}")
```

---

## Preprocessing Failures

### Problem: Invalid SMILES Strings

**Error Messages**:
```
Warning: Failed to parse SMILES for molecule 42
RDKit ERROR: Invalid SMILES string
```

**Solutions**:

**Solution 1: Skip invalid molecules** (default behavior)
```python
# Pipeline automatically skips invalid SMILES
# Check logs for skipped molecules
```

**Solution 2: Pre-validate SMILES**
```python
from rdkit import Chem
import pandas as pd

df = pd.read_csv('HCEPDB/data_calcqcset1.csv')

# Add validity column
df['valid_smiles'] = df['smiles'].apply(
    lambda s: Chem.MolFromSmiles(s) is not None if pd.notna(s) else False
)

# Filter to valid only
df_valid = df[df['valid_smiles']].drop(columns=['valid_smiles'])
print(f"Valid molecules: {len(df_valid)} / {len(df)}")

# Save cleaned dataset
df_valid.to_csv(' HCEPDB/data_calcqcset1_valid.csv', index=False)
```

**Solution 3: Sanitize SMILES**
```python
from rdkit import Chem

def sanitize_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol)  # Canonical SMILES
    except:
        pass
    return None

df['smiles'] = df['smiles'].apply(sanitize_smiles)
df = df.dropna(subset=['smiles'])
```

---

### Problem: Graph Construction Failures

**Error Messages**:
```
ValueError: Graph has 0 nodes
RuntimeError: Edge index out of bounds
```

**Solutions**:

**Solution 1: Check molecule size**
```python
# Skip molecules that are too large or too small
from rdkit import Chem

def is_valid_molecule(smiles, min_atoms=3, max_atoms=200):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    num_atoms = mol.GetNumAtoms()
    return min_atoms <= num_atoms <= max_atoms

df['valid'] = df['smiles'].apply(is_valid_molecule)
df = df[df['valid']]
```

**Solution 2: Handle edge cases in graph builder**
```python
# In preprocessing/graph_builder.py
def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None or mol.GetNumAtoms() == 0:
        return None  # Skip invalid molecules
    
    # Continue with graph construction...
```

**Solution 3: Validate edge indices**
```python
import torch

# Check edge_index is valid
num_nodes = graph.num_nodes
edge_index = graph.edge_index

if edge_index.max() >= num_nodes:
    print("ERROR: Edge index out of bounds!")
    # Fix: Recompute edge indices
```

---

### Problem: Feature Engineering Errors

**Error Messages**:
```
ValueError: NaN values in engineered features
AttributeError: 'NoneType' object has no attribute 'GetNumAtoms'
```

**Solutions**:

**Solution 1: Handle None molecules**
```python
# In feature_engineering.py
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {feature: np.nan for feature in feature_names}
    
    # Calculate features...
```

**Solution 2: Fill missing values**
```python
# After feature engineering
df = df.fillna(df.mean())  # Fill with column means

# Or drop rows with missing features
df = df.dropna()
```

**Solution 3: Add error handling**
```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def safe_descriptor(smiles, descriptor_func):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return descriptor_func(mol)
    except Exception as e:
        print(f"Error computing descriptor: {e}")
    return np.nan
```

---

## Memory Issues

### Problem: Out of Memory During Preprocessing

**Error Messages**:
```
MemoryError: Unable to allocate array
Killed
```

**Solutions**:

**Solution 1: Process in batches**
```python
# Modify pipeline to process in chunks
chunk_size = 100
for i in range(0, len(df), chunk_size):
    chunk = df.iloc[i:i+chunk_size]
    process_chunk(chunk)
    # Save intermediate results
```

**Solution 2: Disable Graphormer encoding**
```python
# preprocessing/config.py
USE_GRAPHORMER_ENCODING = False  # Saves significant memory
```

**Solution 3: Use smaller dataset**
```python
# Test with subset first
df = df.head(100)  # Process only first 100 molecules
```

**Solution 4: Increase swap space (Linux)**
```bash
# Check current swap
free -h

# Add swap file
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

**Solution 5: Use cloud computing**
- Google Colab: Free 12GB RAM
- Kaggle: Free 16GB RAM
- AWS/Azure: Pay for larger instances

---

### Problem: CUDA Out of Memory

**Error Messages**:
```
RuntimeError: CUDA out of memory. Tried to allocate 2.00 GiB
```

**Solutions**:

**Solution 1: Reduce batch size**
```python
# In training script
train_loader = DataLoader(train_data, batch_size=8)  # Reduce from 32
```

**Solution 2: Use gradient accumulation**
```python
effective_batch_size = 32
actual_batch_size = 8
accumulation_steps = effective_batch_size // actual_batch_size

for i, batch in enumerate(train_loader):
    loss = model(batch) / accumulation_steps
    loss.backward()
    
    if (i + 1) % accumulation_steps == 0:
        optimizer.step()
        optimizer.zero_grad()
```

**Solution 3: Use smaller model**
```python
model = GCN_OPV_Predictor(
    hidden_dim=64,  # Instead of 128
    num_layers=2    # Instead of 3
)
```

**Solution 4: Clear cache regularly**
```python
import torch

for epoch in range(num_epochs):
    train_one_epoch()
    torch.cuda.empty_cache()  # Clear unused memory
```

**Solution 5: Use mixed precision training**
```python
from torch.cuda.amp import autocast, GradScaler

scaler = GradScaler()

for batch in train_loader:
    optimizer.zero_grad()
    
    with autocast():  # Automatic mixed precision
        output = model(batch)
        loss = criterion(output, batch.y)
    
    scaler.scale(loss).backward()
    scaler.step(optimizer)
    scaler.update()
```

**Solution 6: Use CPU instead**
```python
device = torch.device('cpu')
model = model.to(device)
```

---

## Model Training Problems

### Problem: Model Not Learning (Loss Not Decreasing)

**Symptoms**:
```
Epoch 1: Loss = 10.5
Epoch 10: Loss = 10.4
Epoch 50: Loss = 10.3
```

**Solutions**:

**Solution 1: Check learning rate**
```python
# Try different learning rates
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)   # Higher
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001) # Lower
```

**Solution 2: Verify data scaling**
```python
# Ensure features and targets are normalized
import pickle

with open('preprocessed_data/scalers/feature_scaler.pkl', 'rb') as f:
    feature_scaler = pickle.load(f)

with open('preprocessed_data/scalers/target_scaler.pkl', 'rb') as f:
    target_scaler = pickle.load(f)

# Apply scaling before training
X_scaled = feature_scaler.transform(X)
y_scaled = target_scaler.transform(y)
```

**Solution 3: Check for NaN values**
```python
# After each forward pass
if torch.isnan(loss):
    print("NaN loss detected!")
    print("Model outputs:", output)
    print("Targets:", batch.y)
    break
```

**Solution 4: Use gradient clipping**
```python
# Prevent exploding gradients
torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
```

**Solution 5: Verify model architecture**
```python
# Print model to check layers
print(model)

# Check parameter count
num_params = sum(p.numel() for p in model.parameters())
print(f"Total parameters: {num_params:,}")
```

---

### Problem: Overfitting (Train Loss Low, Val Loss High)

**Symptoms**:
```
Epoch 50:
  Train Loss = 0.2
  Val Loss = 1.5
```

**Solutions**:

**Solution 1: Add dropout**
```python
model = GCN_OPV_Predictor(
    dropout=0.3  # Increase from 0.2
)
```

**Solution 2: Add L2 regularization**
```python
optimizer = torch.optim.Adam(
    model.parameters(),
    lr=0.001,
    weight_decay=1e-4  # L2 penalty
)
```

**Solution 3: Early stopping**
```python
best_val_loss = float('inf')
patience = 10
patience_counter = 0

for epoch in range(max_epochs):
    train_loss = train_one_epoch()
    val_loss = validate()
    
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        patience_counter = 0
        torch.save(model.state_dict(), 'best_model.pt')
    else:
        patience_counter += 1
        if patience_counter >= patience:
            print("Early stopping!")
            break
```

**Solution 4: Reduce model complexity**
```python
model = GCN_OPV_Predictor(
    hidden_dim=64,   # Reduce from 128
    num_layers=2     # Reduce from 3
)
```

**Solution 5: Use more training data**
```python
# Increase training set size
# preprocessing/config.py
TRAIN_RATIO = 0.8  # Instead of 0.7
VAL_RATIO = 0.1
TEST_RATIO = 0.1
```

---

### Problem: Slow Training Speed

**Solutions**:

**Solution 1: Use GPU**
```bash
# Install CUDA-enabled PyTorch
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
```

**Solution 2: Increase batch size**
```python
train_loader = DataLoader(train_data, batch_size=64)  # From 32
```

**Solution 3: Use DataLoader workers**
```python
train_loader = DataLoader(
    train_data,
    batch_size=32,
    num_workers=4,  # Parallel data loading
    pin_memory=True  # Faster GPU transfer
)
```

**Solution 4: Profile code**
```python
import cProfile

cProfile.run('train_model()', 'training_profile.prof')

# Analyze results
import pstats
stats = pstats.Stats('training_profile.prof')
stats.sort_stats('cumulative')
stats.print_stats(20)  # Top 20 slowest functions
```

---

## Output Interpretation

### Problem: Understanding Warning Messages

**Warning**: `Outlier detected in feature 'homo_lumo_gap': value=5.2 (mean=1.8, std=0.3)`

**Meaning**: Statistical outlier (>3 standard deviations from mean)

**Action**: 
- Inspect the molecule: `df[df['homo_lumo_gap'] == 5.2]`
- May be a valid edge case or data error
- Not automatically removed unless `REMOVE_OUTLIERS = True` in config

---

**Warning**: `Skipped 5 molecules with missing targets`

**Meaning**: Some molecules don't have PCE/Voc/Jsc values

**Action**: Normal behavior - these molecules can't be used for supervised learning

---

**Warning**: `Scaffold splitting resulted in unbalanced splits: train=75%, val=15%, test=10%`

**Meaning**: Exact ratios couldn't be achieved due to scaffold constraints

**Action**: Usually acceptable - scaffolds are more important than exact ratios

---

### Problem: Unexpected Output Directory Structure

**Issue**: Missing subdirectories or files

**Solution**: Check preprocessing completed all 10 steps

```bash
# Verify all output files exist
ls preprocessed_data/molecules.csv
ls preprocessed_data/graphs/ | wc -l  # Should equal number of molecules
ls preprocessed_data/splits/
ls preprocessed_data/scalers/
ls preprocessed_data/metadata.json
```

---

## Platform-Specific Issues

### macOS: Permission Denied Errors

```bash
# Fix script permissions
chmod +x scripts/*.sh

# Fix conda permissions (if installed system-wide)
sudo chown -R $(whoami) /path/to/anaconda3
```

### Windows: Path Issues

```python
# Use pathlib for cross-platform paths
from pathlib import Path

DATA_PATH = Path("HCEPDB") / "data_calcqcset1.csv"
OUTPUT_DIR = Path("preprocessed_data")
```

### Linux: Matplotlib Display Issues

```python
# If running on server without display
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
```

---

## Getting Further Help

If these solutions don't resolve your issue:

1. **Check logs**: Review `preprocessing.log` for detailed error messages
2. **Search documentation**: [User Manual](user-manual.md), [FAQ](faq.md)
3. **GitHub Issues**: Search existing issues or create a new one
4. **Contact team**: See [team members](README.md#team-members)

---

[← Back to Documentation Home](README.md)
