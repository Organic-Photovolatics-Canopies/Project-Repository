# User Manual

[← Back to Documentation Home](README.md)

This is the comprehensive technical reference for the OMID USA system. This manual provides detailed information about the architecture, preprocessing pipeline, data schemas, and model implementations.

## Table of Contents

- [System Architecture](#system-architecture)
- [Preprocessing Pipeline](#preprocessing-pipeline)
- [Configuration Options](#configuration-options)
- [Data Schema](#data-schema)
- [Model Architecture](#model-architecture)
- [API Reference](#api-reference-coming-soon)
- [Advanced Usage](#advanced-usage)

---

## System Architecture

### Overview

OMID USA follows a modular pipeline architecture designed for scalability and maintainability. The system is organized into three main layers:

```
┌─────────────────────────────────────────┐
│         User Interface Layer            │
│  (Web App - Coming Soon / CLI)          │
└─────────────────┬───────────────────────┘
                  │
┌─────────────────▼───────────────────────┐
│      Machine Learning Layer             │
│  • GNN Models (Graphormer, GCN)         │
│  • Training & Inference                 │
│  • Multi-task Prediction                │
└─────────────────┬───────────────────────┘
                  │
┌─────────────────▼───────────────────────┐
│       Data Processing Layer             │
│  • Preprocessing Pipeline               │
│  • Graph Building                       │
│  • Feature Engineering                  │
└─────────────────────────────────────────┘
```

### Design Diagrams

The project includes three levels of design diagrams (see [Design Diagrams/](../Design%20Diagrams/)):

**D0 - Context Diagram**: High-level system context showing user interactions
**D1 - Container Diagram**: Major components and their communication
**D2 - Component Diagram**: Detailed component relationships within the ML system

### Technology Stack

| Layer | Technologies |
|-------|-------------|
| **Data Processing** | Python, pandas, NumPy, RDKit |
| **Graph Representation** | PyTorch Geometric, NetworkX |
| **Machine Learning** | PyTorch, scikit-learn |
| **Chemistry** | RDKit, SMILES, molecular descriptors |
| **Data Storage** | CSV, PyTorch tensors (.pt), pickle |
| **Quantum Chemistry** | DFT calculations from OPV2D |

---

## Preprocessing Pipeline

The preprocessing pipeline transforms raw quantum chemistry data into graph representations suitable for GNN training. It consists of 10 sequential steps.

### Pipeline Architecture

```python
# File: preprocessing/pipeline.py

class PreprocessingPipeline:
    """
    End-to-end pipeline for OPV molecular data preprocessing.
    
    Steps:
        1. Data Loading
        2. Data Integration
        3. Graph Building
        4. Feature Engineering
        5. Target Preparation
        6. Data Validation
        7. Scaffold Splitting
        8. Feature Scaling
        9. Graphormer Encoding
        10. Saving
    """
```

### Step-by-Step Breakdown

#### Step 1: Data Loading

**Module**: `preprocessing/data_loader.py`

Loads CSV files from the OPV2D dataset:

```python
def load_data(data_path: str) -> pd.DataFrame:
    """
    Load molecular data from CSV file.
    
    Args:
        data_path: Path to CSV file
        
    Returns:
        DataFrame with columns: smiles, homo, lumo, gap, ...
    """
```

**Input**: CSV files with quantum chemistry calculations and device data
**Output**: pandas DataFrame with molecular data
**Key Columns**: `smiles`, `homo`, `lumo`, `gap`, `pce`, `voc`, `jsc`

#### Step 2: Data Integration

**Module**: `preprocessing/data_integration.py`

Integrates and cleans molecular data:

```python
def integrate_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Process and integrate molecular property data.
    
    Returns:
        Integrated DataFrame with all features and targets.
    """
```

**Processing**: Handles missing values, validates data quality
**Output**: Single DataFrame with all electronic features and OPV targets

#### Step 3: Graph Building

**Module**: `preprocessing/graph_builder.py`

Converts SMILES strings to PyTorch Geometric graph objects:

```python
def smiles_to_graph(smiles: str) -> Data:
    """
    Convert SMILES string to PyG graph representation.
    
    Args:
        smiles: SMILES molecular representation
        
    Returns:
        PyG Data object with:
            - x: Node features [num_atoms, 7]
            - edge_index: Edge connectivity [2, num_bonds]
            - edge_attr: Edge features [num_bonds, 3]
    """
```

**Node Features** (7 dimensions per atom):

| Feature | Description | Values |
|---------|-------------|---------|
| `atomic_number` | Element type | 1-118 (H=1, C=6, N=7, ...) |
| `degree` | Number of bonds | 0-6 |
| `formal_charge` | Ionic charge | -2 to +2 |
| `hybridization` | Orbital hybridization | 0=s, 1=sp, 2=sp2, 3=sp3, 4=other |
| `is_aromatic` | Aromatic atom | 0 or 1 |
| `num_hydrogens` | Attached H atoms | 0-4 |
| `is_in_ring` | Part of ring structure | 0 or 1 |

**Edge Features** (3 dimensions per bond):

| Feature | Description | Values |
|---------|-------------|---------|
| `bond_type` | Chemical bond type | 1=single, 2=double, 3=triple, 1.5=aromatic |
| `is_conjugated` | Conjugated bond | 0 or 1 |
| `is_in_ring` | Part of ring | 0 or 1 |

**Chemistry Background**:
- SMILES (Simplified Molecular Input Line Entry System) is a text notation for molecular structures
- Example: `C1=CC=CC=C1` represents benzene (6-carbon aromatic ring)
- Atoms become graph nodes, bonds become edges

**Optical Properties Note**:
- Absorption spectra will be generated using TD-DFT calculations
- These predict which wavelengths each molecule absorbs
- Critical for transparency-optimized OPV design

#### Step 4: Feature Engineering

**Module**: `preprocessing/feature_engineering.py`

Creates 11 derived molecular descriptors:

```python
def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create engineered features from quantum chemistry properties.
    
    Features created:
        1. homo_lumo_gap: Electronic bandgap
        2. lumo_squared: Quadratic LUMO term
        3. homo_times_lumo: Interaction term
        4. dipole_magnitude: Total dipole moment
        5. heteroatom_ratio: Non-C/H atom ratio
        6. aromatic_atom_count: Aromaticity measure
        7. molecular_weight: MW from SMILES
        8. heavy_atom_count: Non-H atoms
        9. ring_count: Number of rings
        10. rotatable_bonds: Flexibility measure
        11. aromatic_density: Aromaticity per atom
    """
```

**Feature Categories**:

1. **Electronic Properties**: Derived from HOMO/LUMO energies
   - Important for charge transport in OPVs
   - `homo_lumo_gap = lumo - homo`

2. **Structural Properties**: Derived from SMILES via RDKit
   - Molecular size, flexibility, complexity
   - Ring systems affect optical properties

3. **Chemical Composition**: Heteroatom content
   - N, S, Se atoms are crucial for OPV performance
   - `heteroatom_ratio = (num_heteroatoms) / (total_atoms)`

#### Step 5: Target Preparation

**Module**: `preprocessing/target_preparation.py`

Prepares the three regression targets:

```python
def prepare_targets(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract and prepare PCE, Voc, Jsc targets.
    
    Targets:
        - PCE: Power Conversion Efficiency (%)
        - Voc: Open-Circuit Voltage (V)
        - Jsc: Short-Circuit Current Density (mA/cm²)
        
    Filters:
        - Removes molecules with missing targets
        - Validates target ranges
    """
```

**Target Property Explanations**:

**PCE (Power Conversion Efficiency)**:
- Measures overall performance: PCE = (Voc × Jsc × FF) / P_in
- FF = Fill Factor (device characteristic)
- P_in = Incident light power
- **Typical range**: 3-15% for OPVs (20-25% for silicon)
- **Target**: Maximize while maintaining transparency

**Voc (Open-Circuit Voltage)**:
- Maximum voltage when no current flows
- Related to HOMO-LUMO gap: Voc ≈ (E_LUMO - E_HOMO) - losses
- **Typical range**: 0.5-1.2 V
- **Target**: Higher is better

**Jsc (Short-Circuit Current Density)**:
- Current when voltage is zero
- Depends on light absorption and charge extraction
- **Typical range**: 5-20 mA/cm²
- **Target**: Higher is better, but conflicts with transparency

#### Step 6: Data Validation

**Module**: `preprocessing/data_validation.py`

Validates and cleans the dataset:

```python
def validate_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform data quality checks.
    
    Validations:
        - SMILES validity (RDKit parsing)
        - Duplicate detection (InChIKey)
        - Outlier detection (3-sigma rule)
        - Missing value handling
        - Feature range checks
    """
```

**Validation Checks**:

1. **SMILES Validity**: Can RDKit parse the SMILES?
   - Invalid SMILES are removed (~1-2% typically)

2. **Duplicate Detection**: Multiple entries for the same molecule
   - Keep first occurrence

3. **Outlier Detection**: Statistical outliers in features/targets
   - Flag but don't automatically remove (may be valid edge cases)

4. **Missing Values**: Incomplete quantum chemistry calculations
   - Remove molecules with missing critical features

5. **Physical Constraints**: Chemically unrealistic values
   - e.g., negative HOMO-LUMO gap

#### Step 7: Scaffold Splitting

**Module**: `preprocessing/data_splitting.py`

Splits data by molecular scaffold for better generalization:

```python
def scaffold_split(
    df: pd.DataFrame,
    train_ratio: float = 0.7,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Split dataset by Bemis-Murcko scaffolds.
    
    Ensures train/val/test sets have structurally different molecules,
    preventing data leakage and testing true generalization.
    
    Returns:
        train_indices, val_indices, test_indices
    """
```

**Why Scaffold Splitting?**

Traditional random splitting can lead to data leakage when molecules in train and test sets are structurally similar. Scaffold splitting groups molecules by their core structure (Bemis-Murcko scaffold) and assigns entire groups to train/val/test.

**Example**:
```
Molecule A: C1=CC=C(C=C1)C2=CC=CC=C2  (biphenyl scaffold)
Molecule B: C1=CC=C(C=C1)C2=CC=C(C=C2)O  (substituted biphenyl)
→ Both go in the same split (train OR test, not both)
```

**Split Ratios**:
- Training: 70% (~130 molecules)
- Validation: 15% (~28 molecules)
- Test: 15% (~28 molecules)

#### Step 8: Feature Scaling

**Module**: `preprocessing/feature_scaling.py`

Normalizes features and targets using StandardScaler:

```python
def scale_features(
    train_features: np.ndarray,
    val_features: np.ndarray,
    test_features: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, StandardScaler]:
    """
    Scale features using StandardScaler fitted on training data only.
    
    Formula: x_scaled = (x - mean) / std
    
    Critical: Fit on training data, transform all sets.
    """
```

**Scaling Strategy**:

1. **Fit on training set only**: Learn mean and std from training data
2. **Transform all sets**: Apply same transformation to val and test
3. **Separate scalers**: Different scalers for features and targets
4. **Save scalers**: Store for inference on new molecules

**Why Scale?**
- Neural networks train better with normalized inputs (mean=0, std=1)
- Prevents features with large magnitudes from dominating
- Helps gradient descent converge faster

#### Step 9: Graphormer Encoding

**Module**: `preprocessing/graphormer_encoding.py`

Adds Graphormer-specific encodings:

```python
def encode_graphormer(graph: Data) -> Data:
    """
    Add spatial encoding, centrality, and attention bias for Graphormer.
    
    Additions:
        - Spatial encoding: Shortest path distances between all node pairs
        - Centrality encoding: Node importance (degree centrality)
        - Edge encoding: Spatial distance along edges
        - Attention bias: Graph distance-based attention weights
    """
```

**Graphormer Enhancements**:

Graphormer is a transformer-based GNN that requires additional encodings:

1. **Spatial Encoding**: Shortest path distances
   - Helps model understand global graph structure
   - Matrix of size [num_nodes, num_nodes]

2. **Centrality Encoding**: Node importance
   - Degree centrality, eigenvector centrality
   - Identifies key atoms in the molecule

3. **Edge Encodings**: Path-based edge features
   - Captures multi-hop relationships

4. **Attention Bias**: Pre-computed attention weights
   - Guides transformer attention mechanism

**Reference**: "Do Transformers Really Perform Bad for Graph Representation?" (NeurIPS 2021)

#### Step 10: Saving

**Module**: `preprocessing/pipeline.py`

Saves all processed data:

```python
def save_preprocessed_data(output_dir: str, **kwargs):
    """
    Save all pipeline outputs.
    
    Saved files:
        - molecules.csv: Master feature table
        - graphs/*.pt: Individual PyG graphs
        - graphormer_inputs/*.pt: Graphormer-encoded graphs
        - scalers/*.pkl: Feature and target scalers
        - splits/*.npy: Train/val/test indices
        - metadata.json: Configuration and statistics
    """
```

---

## Configuration Options

### Main Configuration File

**File**: `preprocessing/config.py`

```python
# Data paths
DATA_PATH = "path/to/opv2d/data.csv"
OUTPUT_DIR = "preprocessed_data/"

# Dataset splitting
TRAIN_RATIO = 0.7
VAL_RATIO = 0.15
TEST_RATIO = 0.15
SPLIT_METHOD = "scaffold"  # Options: "scaffold", "random"

# Feature engineering
USE_ENGINEERED_FEATURES = True
NORMALIZE_FEATURES = True

# Graph building
MAX_ATOMS = 200  # Maximum atoms per molecule
INCLUDE_EDGE_FEATURES = True

# Validation
REMOVE_OUTLIERS = False  # Set True to automatically remove outliers
OUTLIER_STD_THRESHOLD = 3.0

# Graphormer encoding
USE_GRAPHORMER_ENCODING = True
MAX_SPATIAL_DISTANCE = 20  # Maximum path length for spatial encoding

# Logging
VERBOSE = True
LOG_FILE = "preprocessing.log"
```

### Customization Examples

**Example 1: Use Random Splitting Instead of Scaffold**

```python
# config.py
SPLIT_METHOD = "random"
```

**Example 2: Increase Training Set Size**

```python
# config.py
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1
TEST_RATIO = 0.1
```

**Example 3: Disable Graphormer Encoding (for faster preprocessing)**

```python
# config.py
USE_GRAPHORMER_ENCODING = False
```

**Example 4: Process Custom Dataset**

```python
# config.py
DATA_PATH = "path/to/your/custom_data.csv"
OUTPUT_DIR = "custom_preprocessed/"
```

---

## Data Schema

### Input Data Schema (OPV2D)

Expected CSV columns:

| Column | Type | Description | Required |
|--------|------|-------------|----------|
| `smiles` | string | SMILES molecular representation | Yes |
| `InChIKey` or `id` | string | Unique molecular identifier | Yes |
| `homo` | float | HOMO energy (eV) | Yes |
| `lumo` | float | LUMO energy (eV) | Yes |
| `gap` | float | HOMO-LUMO gap (eV) | Yes |
| `pce` | float | Power Conversion Efficiency (%) | Yes |
| `voc` | float | Open-Circuit Voltage (V) | Yes |
| `jsc` | float | Short-Circuit Current Density (mA/cm²) | Yes |
| `optical_gap` | float | Optical bandgap (eV) | No |
| `absorption` | float | Absorption properties | No |

### Output Data Schema

#### `molecules.csv`

Master feature table with all processed features:

| Column | Type | Description |
|--------|------|-------------|
| `molecule_id` | int | Sequential molecule index (0, 1, 2, ...) |
| `smiles` | string | Original SMILES |
| `InChIKey` | string | Molecular identifier |
| `homo` | float | HOMO energy (eV) |
| `lumo` | float | LUMO energy (eV) |
| `homo_lumo_gap` | float | Engineered: LUMO - HOMO |
| `lumo_squared` | float | Engineered: LUMO² |
| `homo_times_lumo` | float | Engineered: HOMO × LUMO |
| `dipole_magnitude` | float | Engineered: \|\|dipole\|\| |
| `heteroatom_ratio` | float | Engineered: heteroatoms / total atoms |
| `aromatic_atom_count` | float | Engineered: number of aromatic atoms |
| `molecular_weight` | float | Engineered: molecular weight (g/mol) |
| `heavy_atom_count` | float | Engineered: non-H atoms |
| `ring_count` | float | Engineered: number of rings |
| `rotatable_bonds` | float | Engineered: rotatable bond count |
| `aromatic_density` | float | Engineered: aromatic atoms / total atoms |
| `pce` | float | Target: PCE (%) |
| `voc` | float | Target: Voc (V) |
| `jsc` | float | Target: Jsc (mA/cm²) |
| `split` | string | Dataset split: "train", "val", or "test" |

#### PyG Graph Objects (`graphs/mol_*.pt`)

```python
# Load example
import torch
graph = torch.load('preprocessed_data/graphs/mol_0.pt')

# Structure
Data(
    x: Tensor[num_nodes, 7],         # Node features
    edge_index: Tensor[2, num_edges], # Edge connectivity (COO format)
    edge_attr: Tensor[num_edges, 3],  # Edge features
    y: Tensor[3],                     # Targets: [PCE, Voc, Jsc]
    smiles: str,                      # Original SMILES
    num_nodes: int,                   # Number of atoms
    num_edges: int                    # Number of bonds
)
```

#### Metadata (`metadata.json`)

```json
{
  "preprocessing_version": "1.0",
  "preprocessing_date": "2026-02-10",
  "dataset_source": "OPV2D",
  "num_molecules": 15000,
  "num_valid_molecules": 14850,
  "train_size": 10395,
  "val_size": 2228,
  "test_size": 2227,
  "split_method": "scaffold",
  "feature_dim": 7,
  "edge_feature_dim": 3,
  "target_dim": 3,
  "num_engineered_features": 11,
  "config": {
    "DATA_PATH": "OPV2D/data/opv_data.csv",
    "TRAIN_RATIO": 0.7,
    "VAL_RATIO": 0.15,
    "TEST_RATIO": 0.15
  }
}
```

---

## Model Architecture

### Baseline: Graph Convolutional Network (GCN)

**Current Performance**:
- R² Score: 0.4335
- RMSE: 0.7259
- Training time: ~5-10 minutes on CPU

**Architecture**:

```python
import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool

class GCN_OPV_Predictor(nn.Module):
    """
    Baseline GCN for multi-task OPV property prediction.
    """
    def __init__(
        self,
        node_feature_dim: int = 7,
        hidden_dim: int = 128,
        num_layers: int = 3,
        dropout: float = 0.2,
        num_targets: int = 3
    ):
        super().__init__()
        
        # Graph convolutional layers
        self.conv1 = GCNConv(node_feature_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, hidden_dim)
        
        # Batch normalization
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)
        self.bn3 = nn.BatchNorm1d(hidden_dim)
        
        # Dropout
        self.dropout = nn.Dropout(dropout)
        
        # Prediction head
        self.fc1 = nn.Linear(hidden_dim, hidden_dim // 2)
        self.fc2 = nn.Linear(hidden_dim // 2, num_targets)
        
    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        
        # Graph convolutions with ReLU and batch norm
        x = self.conv1(x, edge_index)
        x = self.bn1(x)
        x = torch.relu(x)
        x = self.dropout(x)
        
        x = self.conv2(x, edge_index)
        x = self.bn2(x)
        x = torch.relu(x)
        x = self.dropout(x)
        
        x = self.conv3(x, edge_index)
        x = self.bn3(x)
        x = torch.relu(x)
        
        # Global pooling (graph-level representation)
        x = global_mean_pool(x, batch)
        
        # Prediction head
        x = torch.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)  # Output: [batch_size, 3] (PCE, Voc, Jsc)
        
        return x
```

**Training Configuration**:

```python
# Hyperparameters
model = GCN_OPV_Predictor(
    node_feature_dim=7,
    hidden_dim=128,
    num_layers=3,
    dropout=0.2
)

optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
criterion = nn.MSELoss()  # Mean Squared Error for regression
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, mode='min', factor=0.5, patience=10
)

# Training loop: 100-200 epochs
```

### Advanced: Graphormer

**Status**: Under development

**Architecture**: Transformer-based Graph Neural Network

**Key Features**:
- Self-attention mechanism for capturing long-range dependencies
- Spatial encoding for global graph structure understanding
- Centrality encoding for node importance
- Superior performance on molecular property prediction benchmarks

**Planned Architecture**:

```python
# Pseudocode - full implementation in progress
class Graphormer_OPV_Predictor(nn.Module):
    """
    Graphormer model for OPV property prediction.
    """
    def __init__(
        self,
        node_feature_dim: int = 7,
        num_heads: int = 8,
        num_layers: int = 12,
        hidden_dim: int = 256,
        ffn_dim: int = 1024,
        dropout: float = 0.1
    ):
        super().__init__()
        
        # Spatial encoding embedding
        self.spatial_encoder = SpatialEncoder(max_distance=20)
        
        # Centrality encoding embedding
        self.centrality_encoder = CentralityEncoder()
        
        # Transformer encoder layers
        self.transformer_layers = nn.ModuleList([
            GraphormerLayer(
                hidden_dim=hidden_dim,
                num_heads=num_heads,
                ffn_dim=ffn_dim,
                dropout=dropout
            )
            for _ in range(num_layers)
        ])
        
        # Prediction head
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 3)  # PCE, Voc, Jsc
        )
```

---

## API Reference (Coming Soon)

Placeholder for future web API documentation.

### Planned Endpoints

**POST `/predict`**
```json
{
  "smiles": "C1=CC=C(C=C1)C2=CC=CC=C2",
  "return_features": false
}
```

Response:
```json
{
  "predictions": {
    "pce": 8.24,
    "voc": 0.92,
    "jsc": 12.5
  },
  "molecule_id": "uuid-here",
  "timestamp": "2026-02-10T12:34:56Z"
}
```

---

## Example Implementation

### Kaggle PCE Predictor Notebook

**Interactive Demo**: [https://www.kaggle.com/code/omrjad/pce-predictor](https://www.kaggle.com/code/omrjad/pce-predictor)

This notebook demonstrates a complete end-to-end workflow for PCE prediction:
- Data loading and preprocessing
- Graph neural network implementation
- Model training and validation
- Performance evaluation
- Interactive prediction interface

**Key Features**:
- Self-contained Jupyter notebook
- Reproducible results
- Commented code for learning
- GPU acceleration on Kaggle
- Example predictions with custom molecules

**Usage**:
```python
# The notebook includes sections for:
# 1. Environment setup and imports
# 2. Data loading from OPV2D
# 3. SMILES to graph conversion
# 4. GNN model definition
# 5. Training loop with early stopping
# 6. Evaluation metrics (R², RMSE, MAE)
# 7. Prediction function for new molecules
```

**To Run**:
1. Visit the Kaggle notebook link
2. Click "Copy & Edit" to fork to your account
3. Enable GPU accelerator (Settings → Accelerator → GPU)
4. Click "Run All" or execute cells sequentially

---

## Advanced Usage

### Custom Feature Engineering

Add your own molecular descriptors:

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def custom_features(smiles: str) -> dict:
    """Add custom molecular descriptors."""
    mol = Chem.MolFromSmiles(smiles)
    
    return {
        'logP': Descriptors.MolLogP(mol),
        'TPSA': Descriptors.TPSA(mol),
        'num_valence_electrons': Descriptors.NumValenceElectrons(mol)
    }

# Integrate into preprocessing pipeline
# (modify preprocessing/feature_engineering.py)
```

### Batch Inference

Predict properties for multiple molecules:

```python
import torch
from torch_geometric.loader import DataLoader

# Load model
model = torch.load('models/best_model.pt')
model.eval()

# Load test graphs
test_graphs = [torch.load(f'preprocessed_data/graphs/mol_{i}.pt') 
               for i in test_indices]

# Create data loader
loader = DataLoader(test_graphs, batch_size=32)

# Predict
predictions = []
with torch.no_grad():
    for batch in loader:
        pred = model(batch)
        predictions.append(pred.cpu().numpy())

predictions = np.concatenate(predictions, axis=0)
# Shape: [num_molecules, 3] (PCE, Voc, Jsc)
```

### Inverse Scaling Targets

Convert normalized predictions back to original scale:

```python
import pickle

# Load target scaler
with open('preprocessed_data/scalers/target_scaler.pkl', 'rb') as f:
    target_scaler = pickle.load(f)

# Inverse transform
predictions_original_scale = target_scaler.inverse_transform(predictions_normalized)

print(f"PCE: {predictions_original_scale[0, 0]:.2f}%")
print(f"Voc: {predictions_original_scale[0, 1]:.3f} V")
print(f"Jsc: {predictions_original_scale[0, 2]:.2f} mA/cm²")
```

---

## Related Documentation

- [User Guide](user-guide.md) - Practical workflows and tutorials
- [Installation](installation.md) - Environment setup
- [FAQ](faq.md) - Common questions
- [Troubleshooting](troubleshooting.md) - Problem solving
- [Data Sources](data-sources.md) - Dataset information

---

[← Back to Documentation Home](README.md)
