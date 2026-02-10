# Frequently Asked Questions (FAQ)

[← Back to Documentation Home](README.md)

Find answers to common questions about OMID USA, organic photovoltaics, and the preprocessing pipeline.

## Table of Contents

- [General Questions](#general-questions)
- [Technical Questions](#technical-questions)
- [Domain-Specific Questions](#domain-specific-questions)
- [Troubleshooting Questions](#troubleshooting-questions)

---

## General Questions

### What are Organic Photovoltaics (OPVs)?

**Organic Photovoltaics** are solar cells made from carbon-based materials (organic compounds) rather than traditional silicon. They offer unique advantages:

- **Transparency**: Can be semi-transparent, allowing light through
- **Flexibility**: Can be printed on flexible substrates
- **Lightweight**: Much lighter than silicon panels
- **Low-cost manufacturing**: Solution-processable (printing, coating)
- **Tunable properties**: Chemical structure determines optical/electronic properties

**Disadvantages** compared to silicon:
- Lower efficiency (8-15% vs. 20-25% for silicon)
- Shorter lifespan (5-10 years vs. 25+ years for silicon)
- Less stable under environmental stress

### Why transparent solar panels for agriculture?

**Agrivoltaics** combines agriculture with photovoltaics. Transparent OPV panels can:

1. **Generate electricity** while allowing sunlight through for crops
2. **Protect crops** from excessive heat and UV damage
3. **Reduce water usage** by providing shade and reducing evaporation
4. **Dual land use**: Same land produces both food and energy

**Applications**:
- Vineyard canopies
- Greenhouse roofs
- High-value crop protection (berries, leafy greens)

### What properties does the model predict?

The model predicts three critical OPV performance metrics:

**1. PCE (Power Conversion Efficiency)**
- **What it is**: Percentage of sunlight converted to electricity
- **Formula**: PCE = (Voc × Jsc × FF) / P_incident × 100%
- **Typical range**: 3-15% for OPVs
- **Why it matters**: Primary measure of solar cell performance

**2. Voc (Open-Circuit Voltage)**
- **What it is**: Maximum voltage when no current flows
- **Typical range**: 0.5-1.2 volts
- **Why it matters**: Higher voltage = better energy yield

**3. Jsc (Short-Circuit Current Density)**
- **What it is**: Current per area when voltage is zero
- **Typical range**: 5-20 mA/cm²
- **Why it matters**: Higher current = more power generation

### How accurate are the predictions?

**Current Baseline GCN**:
- R² Score: **0.4335** (explains 43% of variance)
- RMSE: **0.7259**

This means:
- Predictions are moderately correlated with actual values
- Useful for comparing molecules and screening candidates
- Not yet accurate enough for quantitative design without experimental validation

**Expected with Graphormer**: R² > 0.7 (based on literature benchmarks)

**Limitations**:
- Training data from computational predictions (Scharber model), not experiments
- Real-world performance includes device engineering factors not captured
- Best used for relative comparisons, not absolute performance prediction

### Can I train on my own molecular dataset?

**Yes!** The pipeline is designed to be flexible.

**Requirements for your dataset**:

1. **CSV format** with required columns:
   - `smiles`: SMILES molecular representation
   - `InChIKey`: Unique identifier (can be generated)
   - `homo`: HOMO energy (eV)
   - `lumo`: LUMO energy (eV)
   - `pce`: Power Conversion Efficiency (%)
   - `voc`: Open-Circuit Voltage (V)
   - `jsc`: Short-Circuit Current Density (mA/cm²)

2. **Valid SMILES strings**: Must be parseable by RDKit

3. **Complete quantum chemistry data**: HOMO/LUMO energies typically from DFT calculations

**Steps**:
```python
# 1. Edit preprocessing/config.py
DATA_PATH = "path/to/your/custom_dataset.csv"
OUTPUT_DIR = "custom_output/"

# 2. Run pipeline
python -m preprocessing.pipeline

# 3. Train model on your data
```

### How long does preprocessing take?

**Depends on dataset size**:

| Dataset Size | Estimated Time (CPU) | Memory Usage |
|--------------|---------------------|--------------|
| 100 molecules | 1-2 minutes | 2-4 GB |
| 500 molecules | 5-10 minutes | 4-8 GB |
| 1,000 molecules | 15-30 minutes | 8-16 GB |
| 10,000 molecules | 3-5 hours | 16-32 GB |

**Bottlenecks**:
- SMILES → Graph conversion (RDKit)
- Feature engineering (molecular descriptors)
- Graphormer encoding (spatial distances)

**Speed tips**:
- Use multi-core CPU (parallelization in future updates)
- Disable Graphormer encoding if not needed (`USE_GRAPHORMER_ENCODING = False`)
- Use GPU for model training (not preprocessing)

### What's the difference between HCEPDB and OPV2D datasets?

**HCEPDB (Harvard Clean Energy Project Database)**:
- **Our primary dataset**
- 2.3 million molecules (we use a subset)
- DFT calculations at PBE/6-31G* level
- Includes Scharber model predictions for OPV properties
- Reference: [https://www.cepdb.net](https://www.cepdb.net)

**OPV2D (Referenced dataset)**:
- External reference dataset for comparison
- Focus on 2D organic semiconductors
- Different molecular diversity
- Reference: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)

We primarily use HCEPDB but reference OPV2D for benchmarking and literature comparison.

---

## Technical Questions

### Why does preprocessing take so long?

**Main computational costs**:

1. **RDKit SMILES parsing** (20-30% of time)
   - Converting text to molecular objects
   - Computationally intensive for complex molecules

2. **Graph construction** (20-30% of time)
   - Computing node and edge features
   - Building adjacency matrices

3. **Feature engineering** (10-20% of time)
   - Calculating molecular descriptors
   - Ring detection, aromaticity analysis

4. **Graphormer encoding** (30-40% of time)
   - All-pairs shortest paths: O(n³) complexity
   - Most expensive step for large molecules

**Solutions**:
- Disable Graphormer encoding if using baseline GCN
- Process in batches (save intermediate results)
- Use a machine with more CPU cores

### What if I have limited RAM/CPU?

**For Limited RAM (<8 GB)**:

1. Process smaller datasets (use a subset)
2. Reduce batch sizes during training
3. Disable Graphormer encoding
4. Use cloud computing (free options: Google Colab, Kaggle)

**For Limited CPU**:

1. Expect longer preprocessing times
2. Consider cloud computing with better CPUs
3. Let preprocessing run overnight
4. Use pre-processed datasets if available

**Cloud Options**:
- **Google Colab**: Free GPU/CPU, 12 hours runtime
- **Kaggle Notebooks**: Free GPU, 30 hours/week
- **AWS Free Tier**: Limited but sufficient for testing

### How do I update the dataset?

**Option 1: Replace the data file**

```bash
# Backup old data
mv HCEPDB/data_calcqcset1.csv HCEPDB/data_calcqcset1_backup.csv

# Copy new data file
cp /path/to/new_data.csv HCEPDB/data_calcqcset1.csv

# Rerun preprocessing
python -m preprocessing.pipeline
```

**Option 2: Use a different data directory**

```python
# Edit preprocessing/config.py
DATA_PATH = "path/to/new_data.csv"
OUTPUT_DIR = "new_preprocessed_data/"
```

**Option 3: Append new molecules**

```python
import pandas as pd

# Load existing data
old_data = pd.read_csv('HCEPDB/data_calcqcset1.csv')
new_data = pd.read_csv('new_molecules.csv')

# Combine
combined = pd.concat([old_data, new_data], ignore_index=True)

# Remove duplicates by InChIKey
combined = combined.drop_duplicates(subset='InChIKey', keep='first')

# Save
combined.to_csv('HCEPDB/data_combined.csv', index=False)
```

### What's scaffold splitting vs. random splitting?

**Random Splitting**:
- Randomly assigns molecules to train/val/test
- **Problem**: Similar molecules may appear in both train and test
- **Result**: Overly optimistic performance (data leakage)
- **Use when**: Predicting properties of similar molecules

**Scaffold Splitting** (Default):
- Groups molecules by Bemis-Murcko scaffold (core structure)
- Ensures train and test sets have structurally different molecules
- **Advantage**: Tests true generalization to novel structures
- **Result**: More realistic performance estimates
- **Use when**: Designing new molecules (our use case)

**Example**:
```
Random: 
  Train: Molecule A (biphenyl-based)
  Test:  Molecule B (slight variant of biphenyl) ← Data leakage!

Scaffold:
  Train: All biphenyl-based molecules
  Test:  All thiophene-based molecules ← True generalization test
```

**How to change**:
```python
# preprocessing/config.py
SPLIT_METHOD = "random"  # or "scaffold"
```

### How are SMILES strings converted to graphs?

**Step-by-step process**:

1. **Parse SMILES with RDKit**
   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles('C1=CC=CC=C1')  # Benzene
   ```

2. **Extract atoms (graph nodes)**
   ```python
   for atom in mol.GetAtoms():
       features = [
           atom.GetAtomicNum(),      # Element
           atom.GetDegree(),          # Number of bonds
           atom.GetFormalCharge(),    # Charge
           atom.GetHybridization(),   # sp, sp2, sp3
           atom.GetIsAromatic(),      # Aromatic?
           atom.GetTotalNumHs(),      # H atoms
           atom.IsInRing()            # In ring?
       ]
   ```

3. **Extract bonds (graph edges)**
   ```python
   for bond in mol.GetBonds():
       start_atom = bond.GetBeginAtomIdx()
       end_atom = bond.GetEndAtomIdx()
       features = [
           bond.GetBondTypeAsDouble(),  # Single, double, etc.
           bond.GetIsConjugated(),       # Conjugated?
           bond.IsInRing()               # In ring?
       ]
   ```

4. **Build PyTorch Geometric graph**
   ```python
   from torch_geometric.data import Data
   import torch
   
   graph = Data(
       x=torch.tensor(node_features),      # [num_atoms, 7]
       edge_index=torch.tensor(edge_list), # [2, num_bonds]
       edge_attr=torch.tensor(edge_features) # [num_bonds, 3]
   )
   ```

**Why graphs?**
- Molecules have natural graph structure (atoms = nodes, bonds = edges)
- Graph Neural Networks can learn from this structure
- Captures 3D spatial relationships and chemical interactions

### What happens if SMILES parsing fails?

**Causes of parsing failures**:
- Invalid SMILES syntax
- Uncommon elements not supported
- Stereochemistry not specified
- Corrupt data

**Pipeline behavior**:
1. Attempts to parse SMILES with RDKit
2. If parsing fails:
   - Logs a warning with the molecule ID
   - Skips that molecule
   - Continues with remaining molecules
3. Reports summary of skipped molecules

**Example output**:
```
Warning: Failed to parse SMILES for molecule 42: "C1=CC=C(invalid)..."
Processed 185/186 molecules (1 skipped)
```

**Solution**:
- Review skipped molecules in the logs
- Manually correct invalid SMILES
- Or remove those entries from the dataset

---

## Domain-Specific Questions

### What HOMO-LUMO gap is ideal for OPVs?

**HOMO-LUMO gap** (electronic bandgap) determines:
- **Light absorption**: Smaller gap → absorbs more light → higher Jsc
- **Voltage output**: Larger gap → higher Voc

**Ideal range for OPVs**: 1.5-2.0 eV

**Trade-offs**:
- **Gap too small** (<1.3 eV):
  - Absorbs more light (good for Jsc)
  - But low Voc → poor overall PCE
  - And poor chemical stability

- **Gap too large** (>2.5 eV):
  - High Voc (good)
  - But absorbs less sunlight → low Jsc
  - More transparent but less efficient

**For transparent OPVs**: Slightly larger gaps (1.8-2.2 eV) sacrifice some efficiency for transparency

**Related to Voc**:
```
Voc ≈ (|E_LUMO| - |E_HOMO|) - 0.3V
     = HOMO-LUMO_gap - 0.3V
```
The 0.3V loss comes from charge separation and extraction losses.

### How does transparency affect efficiency?

**Fundamental trade-off**:
- More transparency → less light absorbed → lower efficiency
- Higher efficiency → more light absorbed → less transparency

**Transparency mechanisms**:
1. **Selective absorption**: Absorb UV/near-IR, transmit visible light
2. **Thin films**: Reduce material thickness
3. **Larger bandgap**: Only absorb higher-energy photons

**Typical values**:
- **10-20% transparent**: PCE = 12-15% (semi-transparent)
- **30-50% transparent**: PCE = 8-10% (visibly transparent)
- **>60% transparent**: PCE = 3-5% (highly transparent)

**Our use case** (agrivoltaics):
- Target ~30-50% transparency (allow enough light for crops)
- Accept PCE = 8-10% (still economically viable)

### What's the typical PCE for organic vs. silicon panels?

**Silicon Solar Cells**:
- Lab record: ~26.7%
- Commercial: 20-22%
- Mature technology (60+ years development)
- Opaque, rigid, heavy

**Organic Photovoltaics (OPVs)**:
- Lab record: ~18-19%
- Commercial: Not widely commercialized yet
- Promising but still improving
- Can be transparent, flexible, lightweight

**Other Technologies**:
- Perovskite: 25-26% (lab), stability issues
- Thin-film CdTe: 22% (commercial)
- Thin-film CIGS: 20% (commercial)
- Dye-sensitized: 12-13%

**Why use OPVs despite lower efficiency?**
- Specific applications where transparency/flexibility matters
- Lower manufacturing cost potential
- Lighter weight (building integration)
- Tunable optical properties

### What molecules perform best as OPV materials?

**Key structural features**:

**1. Conjugated π-system**
- Alternating single/double bonds
- Enables charge delocalization
- Example: Thiophene, benzene rings

**2. Donor-acceptor architecture**
- Electron-rich donor units
- Electron-poor acceptor units
- Creates intramolecular charge transfer

**3. Heteroatoms**
- Nitrogen, sulfur, selenium
- Tune electronic properties
- Improve charge transport

**4. Planarity**
- Flat molecular structure
- Better π-π stacking
- Enhanced charge mobility

**Example high-performing scaffold**:
- Benzodithiophene (donor) + benzothiadiazole (acceptor)
- Forms conjugated polymers
- PCE > 10% in lab devices

**Molecules to avoid**:
- Too many sp³ carbons (breaks conjugation)
- Bulky side groups (hinders packing)
- Unstable functional groups (poor lifetime)

### Why use Graph Neural Networks instead of traditional ML?

**Traditional ML** (e.g., Random Forest, SVM):
- Requires hand-crafted molecular descriptors
- Fixed-length feature vectors
- Loses structural information
- Limited expressiveness

**Graph Neural Networks**:
- Learn representations directly from molecular structure
- Capture atom-atom interactions
- Variable-size molecules (flexible input)
- State-of-the-art performance on molecular tasks

**Performance comparison** (molecular property prediction):
- Random Forest: R² ≈ 0.3-0.5
- Feed-forward NN: R² ≈ 0.4-0.6
- GCN/GNN: R² ≈ 0.6-0.8
- Graphormer: R² ≈ 0.7-0.9

**Why GNNs work better**:
- Learn chemical structure-property relationships
- Message passing between atoms
- Capture multi-hop interactions
- Inductive bias matches molecular structure

---

## Troubleshooting Questions

### Why am I getting RDKit import errors?

**Common causes**:

**1. RDKit not installed**
```bash
# Solution: Install from conda-forge
conda install -c conda-forge rdkit -y
```

**2. Wrong conda channel**
```bash
# Solution: Add conda-forge channel
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install rdkit -y
```

**3. Python version incompatibility**
```bash
# Solution: Use Python 3.11
conda create -n opv_preprocessing python=3.11 -y
conda activate opv_preprocessing
conda install -c conda-forge rdkit -y
```

**4. Conflicting packages**
```bash
# Solution: Fresh environment
conda env remove -n opv_preprocessing
conda create -n opv_preprocessing python=3.11 -y
conda activate opv_preprocessing
conda install -c conda-forge rdkit -y
```

### "No valid molecules found" error?

**Cause**: All SMILES in the dataset failed to parse

**Solutions**:

**1. Check data file exists**
```bash
ls HCEPDB/data_calcqcset1.csv
```

**2. Verify CSV format**
```python
import pandas as pd
df = pd.read_csv('HCEPDB/data_calcqcset1.csv')
print(df.head())
print(df.columns)  # Should include 'smiles' column
```

**3. Check SMILES column name**
```python
# If column is named differently, rename it
df = df.rename(columns={'Smiles': 'smiles'})
df.to_csv('HCEPDB/data_calcqcset1.csv', index=False)
```

**4. Validate a few SMILES manually**
```python
from rdkit import Chem
smiles = df.iloc[0]['smiles']
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print(f"Invalid SMILES: {smiles}")
```

### CUDA out of memory errors?

**Solutions**:

**1. Reduce batch size**
```python
# In training script
train_loader = DataLoader(train_data, batch_size=16)  # Reduce from 32
```

**2. Use gradient accumulation**
```python
# Accumulate gradients over multiple batches
accumulation_steps = 4
for i, batch in enumerate(train_loader):
    loss = model(batch) / accumulation_steps
    loss.backward()
    
    if (i + 1) % accumulation_steps == 0:
        optimizer.step()
        optimizer.zero_grad()
```

**3. Use CPU for preprocessing**
```python
# preprocessing/config.py
DEVICE = 'cpu'  # Don't use GPU for preprocessing
```

**4. Clear GPU cache**
```python
import torch
torch.cuda.empty_cache()
```

**5. Use smaller model**
```python
# Reduce hidden dimensions
model = GCN_OPV_Predictor(hidden_dim=64)  # Instead of 128
```

### Missing dataset files?

**If `HCEPDB/data_calcqcset1.csv` doesn't exist**:

**Option 1**: Download from repository
```bash
# Check if it's in git history
git pull
```

**Option 2**: Use your own data
```python
# preprocessing/config.py
DATA_PATH = "path/to/your/data.csv"
```

**Option 3**: Download from HCEPDB website
```bash
# Visit https://www.cepdb.net
# Download dataset
# Place in HCEPDB/ directory
```

For more troubleshooting, see the [Troubleshooting Guide](troubleshooting.md).

---

## Still Have Questions?

- 📘 [User Manual](user-manual.md) - Comprehensive technical documentation
- 🔧 [Troubleshooting Guide](troubleshooting.md) - Detailed problem solving
- 💬 Contact the [development team](README.md#team-members)
- 🐛 [Report an issue on GitHub](https://github.com/)

---

[← Back to Documentation Home](README.md)
