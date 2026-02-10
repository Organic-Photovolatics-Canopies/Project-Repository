# User Guide

[← Back to Documentation Home](README.md)

Welcome to the OMID USA User Guide! This guide will help you get started with the OPV materials design pipeline and understand common workflows.

## Table of Contents

- [What is OMID USA?](#what-is-omid-usa)
- [Quick Start (5 Minutes)](#quick-start-5-minutes)
- [Understanding the Workflow](#understanding-the-workflow)
- [Common Workflows](#common-workflows)
  - [Running the Preprocessing Pipeline](#workflow-1-running-the-preprocessing-pipeline)
  - [Exploring Processed Data](#workflow-2-exploring-processed-data)
  - [Training a GNN Model](#workflow-3-training-a-gnn-model)
  - [Using the Web Demo](#workflow-4-using-the-web-demo-coming-soon)
- [Understanding Output](#understanding-output)
- [Tips for Success](#tips-for-success)

---

## What is OMID USA?

**OMID USA** (AI-Powered Design of Organic Photovoltaics Canopies for Agrivoltaics) is a research project that uses artificial intelligence to design better transparent solar panels for agricultural and architectural applications.

### The Problem We're Solving

Traditional solar panels are opaque and block sunlight. This makes them unsuitable for:
- **Greenhouses and vineyards** that need sunlight for crops
- **Building windows** that need natural light
- **Architectural features** where aesthetics matter

**Organic Photovoltaics (OPVs)** can be transparent, flexible, and lightweight—perfect for these applications!

**Our Approach**: We use TD-DFT (Time-Dependent Density Functional Theory) calculations to predict absorption spectra, allowing us to design materials that absorb specific wavelengths while remaining transparent in the visible range.

### Our Solution

We use **Graph Neural Networks** (GNNs) to predict how well different molecular designs will perform as OPV materials. Our system predicts three critical properties:

1. **PCE (Power Conversion Efficiency)**: How efficiently the material converts sunlight to electricity
2. **Voc (Open-Circuit Voltage)**: The maximum voltage the material can produce
3. **Jsc (Short-Circuit Current Density)**: The current generation capability

This accelerates materials discovery—instead of synthesizing and testing thousands of molecules in the lab, we can predict their properties computationally!

---

## Quick Start (5 Minutes)

Get up and running with the preprocessing pipeline in just a few steps.

### Prerequisites

✅ Completed the [Installation Guide](installation.md)

### Run Your First Pipeline

```bash
# 1. Navigate to the project directory
cd Project-Repository

# 2. Activate the conda environment
conda activate opv_preprocessing

# 3. Run the preprocessing pipeline
python -m preprocessing.pipeline
```

### What Happens Next

You'll see progress messages for each of 10 preprocessing steps:

```
[Step 1/10] Loading data...
✓ Loaded molecules from OPV2D dataset

[Step 2/10] Integrating data tables...
✓ Processed molecular data

[Step 3/10] Building graphs from SMILES...
✓ Created molecular graphs

...

[Step 10/10] Saving processed data...
✓ Saved to preprocessed_data/

Pipeline complete! 🎉
```

**Time**: 2-5 minutes depending on your system

### Check Your Results

```bash
# View the output directory
ls preprocessed_data/

# Output:
# molecules.csv  graphs/  graphormer_inputs/  scalers/  splits/  metadata.json
```

🎉 **Success!** You've processed your first molecular dataset.

**Note**: Preprocessing notebooks and updated data files will be added to the repository soon.

---

## Understanding the Workflow

The OMID USA system follows this high-level workflow:

```
┌─────────────────┐
│  Raw Molecular  │
│  Data (SMILES)  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Preprocessing  │  ← You start here
│    Pipeline     │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Graph Dataset  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   GNN Model     │  ← Train Graphormer/GCN
│   Training      │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Predictions:   │
│  PCE, Voc, Jsc  │
└─────────────────┘
```

**Current Release**: Focuses on the preprocessing pipeline and baseline model training.

**Coming Soon**: Web interface for non-technical users to get design recommendations.

---

## Common Workflows

### Workflow 1: Running the Preprocessing Pipeline

**Goal**: Convert raw molecular data into graph representations ready for GNN training.

#### Step-by-Step

**1. Prepare Your Data**

The pipeline expects CSV files with molecular data:

```bash
# Clone OPV2D dataset (recommended)
git clone https://github.com/sunyrain/OPV2D.git

# Or use custom data - place CSV files in Datasets/ directory
mkdir -p Datasets/
# Copy your CSV files there
```

**Note**: Updated preprocessing notebooks for OPV2D will be available soon.

**2. Configure the Pipeline (Optional)**

Edit `preprocessing/config.py` to customize settings:

```python
# Example configurations in config.py
DATA_PATH = "path/to/opv2d/data.csv"  # Data source
OUTPUT_DIR = "preprocessed_data/"      # Output location
TRAIN_RATIO = 0.7                      # Training set size
VAL_RATIO = 0.15                       # Validation set size
TEST_RATIO = 0.15                      # Test set size
```

**3. Run the Pipeline**

```bash
conda activate opv_preprocessing
python -m preprocessing.pipeline
```

**4. Monitor Progress**

The pipeline displays real-time progress:
- ✓ Green checkmarks indicate successful steps
- Warning messages (if any) appear in yellow
- Errors (if any) appear in red

**5. Review Output**

```bash
# Check the molecules dataframe
head preprocessed_data/molecules.csv

# Verify graph files were created
ls preprocessed_data/graphs/ | wc -l  # Should match number of molecules

# View metadata
cat preprocessed_data/metadata.json
```

#### What Gets Created

- **molecules.csv**: Master feature table with all molecular descriptors
- **graphs/**: Individual PyTorch Geometric graph objects (`.pt` files)
- **graphormer_inputs/**: Graphs with Graphormer-specific encoding
- **scalers/**: StandardScaler objects for feature normalization
- **splits/**: Train/validation/test split indices
- **metadata.json**: Configuration and statistics

---

### Workflow 2: Exploring Processed Data

**Goal**: Understand what features were engineered and inspect molecular graphs.

#### Using Python/Jupyter

**1. Load the Processed Data**

```python
import pandas as pd
import torch
from torch_geometric.data import Data

# Load the molecules dataframe
molecules = pd.read_csv('preprocessed_data/molecules.csv')
print(f"Loaded {len(molecules)} molecules")
print(molecules.columns)
```

**2. Inspect Features**

```python
# View summary statistics
print(molecules.describe())

# Check for missing values
print(molecules.isnull().sum())

# View first molecule
print(molecules.iloc[0])
```

**3. Load and Inspect a Graph**

```python
# Load a single molecular graph
graph = torch.load('preprocessed_data/graphs/mol_0.pt')

print(f"Nodes (atoms): {graph.num_nodes}")
print(f"Edges (bonds): {graph.num_edges}")
print(f"Node features shape: {graph.x.shape}")
print(f"Edge features shape: {graph.edge_attr.shape}")
print(f"Target (PCE, Voc, Jsc): {graph.y}")
```

**4. Visualize a Molecule (Optional)**

```python
from rdkit import Chem
from rdkit.Chem import Draw

# Get SMILES string from molecules.csv
smiles = molecules.iloc[0]['smiles']

# Create molecule object
mol = Chem.MolFromSmiles(smiles)

# Draw molecule
img = Draw.MolToImage(mol)
img.show()  # or img.save('molecule.png')
```

**5. Analyze Target Distribution**

```python
import matplotlib.pyplot as plt

# Plot PCE distribution
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.hist(molecules['pce'], bins=30, edgecolor='black')
plt.xlabel('PCE (%)')
plt.ylabel('Count')
plt.title('Power Conversion Efficiency')

plt.subplot(1, 3, 2)
plt.hist(molecules['voc'], bins=30, edgecolor='black')
plt.xlabel('Voc (V)')
plt.ylabel('Count')
plt.title('Open-Circuit Voltage')

plt.subplot(1, 3, 3)
plt.hist(molecules['jsc'], bins=30, edgecolor='black')
plt.xlabel('Jsc (mA/cm²)')
plt.ylabel('Count')
plt.title('Short-Circuit Current Density')

plt.tight_layout()
plt.savefig('target_distribution.png')
plt.show()
```

---

### Workflow 3: Training a GNN Model

**Goal**: Train a Graph Neural Network to predict OPV properties.

> **Note**: This workflow outlines the training process. Full training scripts are under development.

#### Baseline GCN Model

The project includes a baseline Graph Convolutional Network (GCN) implementation:

```python
# Example training workflow (simplified)
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GCNConv, global_mean_pool
import torch.nn.functional as F

# 1. Load preprocessed graphs
graphs = []
for i in range(len(molecules)):
    graph = torch.load(f'preprocessed_data/graphs/mol_{i}.pt')
    graphs.append(graph)

# 2. Load train/val/test splits
train_idx = torch.load('preprocessed_data/splits/train_indices.npy')
val_idx = torch.load('preprocessed_data/splits/val_indices.npy')
test_idx = torch.load('preprocessed_data/splits/test_indices.npy')

train_graphs = [graphs[i] for i in train_idx]
val_graphs = [graphs[i] for i in val_idx]
test_graphs = [graphs[i] for i in test_idx]

# 3. Create data loaders
train_loader = DataLoader(train_graphs, batch_size=32, shuffle=True)
val_loader = DataLoader(val_graphs, batch_size=32)
test_loader = DataLoader(test_graphs, batch_size=32)

# 4. Define model, optimizer, loss
model = YourGCNModel(...)  # See user manual for architecture details
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
criterion = torch.nn.MSELoss()

# 5. Training loop
for epoch in range(100):
    model.train()
    for batch in train_loader:
        optimizer.zero_grad()
        out = model(batch)
        loss = criterion(out, batch.y)
        loss.backward()
        optimizer.step()
    
    # Validation
    model.eval()
    # ... validation code ...

# 6. Evaluate on test set
# ... evaluation code ...
```

#### Graphormer Model (Advanced)

For state-of-the-art performance, use Graphormer-encoded graphs:

```python
# Load Graphormer-encoded graphs
graphormer_graphs = []
for i in range(len(molecules)):
    graph = torch.load(f'preprocessed_data/graphormer_inputs/mol_{i}.pt')
    graphormer_graphs.append(graph)

# Train with Graphormer architecture
# (Full implementation in development)
```

**Current Baseline Performance**:
- R² Score: 0.4335
- RMSE: 0.7259

---

### Workflow 4: Using the Demo

**Goal**: Get OPV design recommendations through an interactive interface.

#### Kaggle Notebook Demo

**Available Now**: [PCE Predictor on Kaggle](https://www.kaggle.com/code/omrjad/pce-predictor)

This interactive Jupyter notebook allows you to:
- Load and explore the OPV dataset
- Run preprocessing steps interactively
- Train GNN models for PCE prediction
- Make predictions on custom molecules
- Visualize results and model performance

**How to Use**:
1. Visit the [Kaggle notebook](https://www.kaggle.com/code/omrjad/pce-predictor)
2. Click "Copy & Edit" to create your own version
3. Run cells sequentially to reproduce results
4. Modify SMILES strings to test your own molecules

#### Web Demo (Coming Soon)

> **Status**: Under development. Full web interface will be deployed to HuggingFace.

#### Planned Features

**1. Simple Design Query**
- Input: Desired transparency level (e.g., 60%)
- Input: Target efficiency (e.g., PCE > 8%)
- Output: Recommended molecular structures with predicted properties

**2. Custom Molecule Prediction**
- Input: SMILES string of your designed molecule
- Output: Predicted PCE, Voc, Jsc values

**3. Batch Processing**
- Upload: CSV file with multiple SMILES structures
- Download: Predictions for all molecules

**4. Interactive Visualization**
- View molecular structures
- Compare predicted properties
- Export results

#### Access

```
📊 Kaggle Notebook (Available): https://www.kaggle.com/code/omrjad/pce-predictor
🌐 Web Demo (Coming Soon): https://huggingface.co/spaces/omid-usa/opv-predictor
```

---

## Understanding Output

### Output Directory Structure

After running the preprocessing pipeline:

```
preprocessed_data/
├── molecules.csv              # Master feature table (all molecules)
├── graphs/                    # PyTorch Geometric graphs
│   ├── mol_0.pt              # Graph for molecule 0
│   ├── mol_1.pt              # Graph for molecule 1
│   └── ...
├── graphormer_inputs/         # Graphormer-encoded graphs
│   ├── mol_0.pt
│   └── ...
├── scalers/                   # Normalization scalers
│   ├── feature_scaler.pkl    # For input features
│   └── target_scaler.pkl     # For PCE/Voc/Jsc targets
├── splits/                    # Dataset splits
│   ├── train_indices.npy     # Training set indices
│   ├── val_indices.npy       # Validation set indices
│   └── test_indices.npy      # Test set indices
└── metadata.json              # Pipeline configuration & stats
```

### Key Files Explained

#### `molecules.csv`

Contains all processed molecules with features:

| Column | Description | Type |
|--------|-------------|------|
| `smiles` | SMILES molecular representation | string |
| `atomic_number` | Element composition encoding | float |
| `degree` | Atomic connectivity | float |
| `homo_energy` | Highest Occupied Molecular Orbital energy | float |
| `lumo_energy` | Lowest Unoccupied Molecular Orbital energy | float |
| `homo_lumo_gap` | Electronic bandgap | float |
| `heteroatom_ratio` | Ratio of non-C/H atoms | float |
| `aromatic_density` | Aromaticity measure | float |
| `pce` | Target: Power Conversion Efficiency (%) | float |
| `voc` | Target: Open-Circuit Voltage (V) | float |
| `jsc` | Target: Short-Circuit Current Density (mA/cm²) | float |

#### Graph Objects (`graphs/mol_*.pt`)

PyTorch Geometric `Data` objects with:

```python
Data(
    x=[num_nodes, 7],           # Node features (atoms)
    edge_index=[2, num_edges],  # Edge connectivity
    edge_attr=[num_edges, 3],   # Edge features (bonds)
    y=[3],                      # Targets: [PCE, Voc, Jsc]
    smiles='C1=CC=C...'         # Original SMILES
)
```

**Node Features** (7 per atom):
1. Atomic number
2. Degree (number of bonds)
3. Formal charge
4. Hybridization (sp, sp2, sp3)
5. Aromaticity (True/False)
6. Hydrogen count
7. Ring membership

**Edge Features** (3 per bond):
1. Bond type (single, double, triple, aromatic)
2. Conjugation (True/False)
3. Ring membership

#### `metadata.json`

Pipeline configuration and statistics:

```json
{
  "num_molecules": 186,
  "num_graphs": 186,
  "train_size": 130,
  "val_size": 28,
  "test_size": 28,
  "feature_dim": 7,
  "edge_feature_dim": 3,
  "target_dim": 3,
  "split_method": "scaffold",
  "preprocessing_date": "2026-02-10"
}
```

---

## Tips for Success

### 🔧 Preprocessing Tips

**1. Start Small**: Use the sample dataset first to verify everything works.

**2. Monitor Memory**: Large datasets may require more RAM. Processing 1000+ molecules can use 8-16 GB.

**3. Check Data Quality**: Review `molecules.csv` for:
   - Missing values (`NaN`)
   - Outliers in targets
   - Invalid SMILES (filtered automatically)

**4. Customize Splits**: Edit `preprocessing/config.py` to change train/val/test ratios.

**5. Save Intermediate Results**: The pipeline saves after each major step, so you can resume if interrupted.

### 🧪 Model Training Tips

**1. Use Scaffold Splitting**: Already implemented—ensures models generalize to new molecular structures.

**2. Normalize Features**: Scalers are provided—always use them to transform features and targets.

**3. Start with Baseline**: Train the simpler GCN before attempting Graphormer.

**4. Monitor Overfitting**: Compare train vs. validation performance regularly.

**5. Use GPU if Available**: Install CUDA-enabled PyTorch for 5-10x speedup (see [Installation Guide](installation.md#gpu-setup-optional)).

### 📊 Data Exploration Tips

**1. Visualize Distributions**: Check target distributions for imbalance or outliers.

**2. Inspect Graphs**: Load and inspect a few graphs to understand the structure.

**3. Feature Correlations**: Analyze correlations between features and targets to understand predictive relationships.

**4. Use Jupyter Notebooks**: Interactive exploration is easier with notebooks.

### ❗ Common Pitfalls

**Pitfall 1**: Running pipeline outside the project root directory
- **Solution**: Always `cd Project-Repository` first

**Pitfall 2**: Forgetting to activate conda environment
- **Solution**: `conda activate opv_preprocessing` before running anything

**Pitfall 3**: Modifying graphs after splitting
- **Solution**: Always split first, then transform separately for train/val/test

**Pitfall 4**: Using unscaled features for training
- **Solution**: Load and apply the scalers from `scalers/`

---

## What's Next?

✅ You're now familiar with the basic workflows!

**Next Steps**:
- **Deep dive**: Read the [User Manual](user-manual.md) for technical details
- **Troubleshooting**: Check [Troubleshooting Guide](troubleshooting.md) if issues arise
- **Questions**: Review the [FAQ](faq.md)
- **Data details**: Read about [Data Sources](data-sources.md)

**For Developers**:
- Explore the `preprocessing/` module source code
- Review the [Design Diagrams](../Design%20Diagrams/) for system architecture
- Check the [Literature Review](../Research/litReview.tex) for background

---

## Need Help?

- 📖 [User Manual](user-manual.md) - Detailed technical reference
- ❓ [FAQ](faq.md) - Common questions
- 🔧 [Troubleshooting](troubleshooting.md) - Problem solving
- 👥 [Team Contacts](README.md#team-members) - Reach out to developers

---

[← Back to Documentation Home](README.md)
