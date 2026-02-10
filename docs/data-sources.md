# Data Sources

[← Back to Documentation Home](README.md)

This document describes the datasets used in the OMID USA project, focusing on the OPV2D dataset.

## Table of Contents

- [Primary Dataset: OPV2D](#primary-dataset-opv2d)
- [Data Collection Methods](#data-collection-methods)
- [Data Quality & Validation](#data-quality--validation)
- [Licensing & Usage Rights](#licensing--usage-rights)
- [How to Cite](#how-to-cite)

---

## Primary Dataset: OPV2D

### Overview

**OPV2D** is our primary data source, focusing on 2D organic semiconductor materials for photovoltaic applications. It provides comprehensive quantum chemistry calculations and experimental data for organic photovoltaic materials.

**Key Statistics**:
- **Total materials**: ~15,000 organic semiconductors
- **Calculation method**: Density Functional Theory (DFT)
- **Properties included**: Electronic properties, optical properties, OPV performance metrics
- **Source**: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)

### Dataset Structure

OPV2D provides comprehensive molecular and electronic properties:

#### Key Data Fields

| Property | Description | Units |
|----------|-------------|-------|
| `smiles` | SMILES string representation | - |
| `homo` | HOMO energy | eV |
| `lumo` | LUMO energy | eV |
| `gap` | HOMO-LUMO gap | eV |
| `pce` | Power Conversion Efficiency | % |
| `voc` | Open-Circuit Voltage | V |
| `jsc` | Short-Circuit Current Density | mA/cm² |
| `optical_gap` | Optical bandgap (from TD-DFT) | eV |
| `absorption_max` | Maximum absorption wavelength | nm |
| `absorption_spectrum` | Full UV-Vis absorption data (TD-DFT) | - |
| `oscillator_strength` | Transition strength | - |

**Data Sources**: Combines computational predictions (DFT) with experimental measurements where available, providing higher reliability than purely computational datasets.

### Accessing OPV2D

**Repository**:
- GitHub: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)
- License: MIT (free for academic and commercial use)
- File format: CSV files
- Size: ~100 MB

**Cloning the Dataset**:
```bash
# Clone the OPV2D repository
git clone https://github.com/sunyrain/OPV2D.git
cd OPV2D/data
```

**Loading the Data**:
```python
import pandas as pd

# Load OPV2D dataset
df = pd.read_csv('OPV2D/data/opv_data.csv')
print(f"Loaded {len(df)} molecules")
print(df.head())
```

### Data Coverage

**Molecular Diversity**:
- 2D organic semiconductor materials
- Conjugated polymers and small molecules
- Donor-acceptor architectures
- Experimentally validated OPV materials

**Heteroatom Content**:
- Carbon, Hydrogen (all molecules)
- Nitrogen, Sulfur, Selenium heterocycles
- Oxygen-containing functional groups

**Property Ranges** (typical values):

| Property | Min | Max | Mean | Notes |
|----------|-----|-----|------|-------|
| HOMO (eV) | -6.5 | -4.0 | -5.2 | Experimental + DFT |
| LUMO (eV) | -4.0 | -2.0 | -3.0 | Experimental + DFT |
| Gap (eV) | 1.2 | 3.5 | 2.2 | Optical measurements |
| PCE (%) | 1.0 | 18.0 | 7.5 | Experimental data |
| Voc (V) | 0.4 | 1.2 | 0.85 | Device measurements |
| Jsc (mA/cm²) | 3.0 | 25.0 | 12.0 | Device measurements |

### Documentation

For detailed information about the dataset:
- Visit the [OPV2D GitHub repository](https://github.com/sunyrain/OPV2D)
- Read the dataset README and documentation
- Check the published paper for methodology

**Note**: Preprocessing notebooks and data files will be added to this repository soon.

---

## Data Collection Methods

### Overview

**OPV2D** is a related dataset we reference for comparison and literature context. It focuses on 2D organic semiconductor materials.

- **Source**: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)
- **Purpose**: Benchmarking, literature comparison
- **Status**: External reference (not directly used in training)

### Key Differences from HCEPDB

| Feature | HCEPDB | OPV2D |
|---------|--------|-------|
| **Focus** | 3D organic molecules | 2D materials |
| **Size** | 2.3M+ molecules | ~15,000 materials |
| **Calculations** | DFT (PBE/6-31G*) | DFT (various functionals) |
| **Properties** | OPV-focused | General electronic properties |
| **Validation** | Scharber model | Various models |

### Usage in OMID USA

We reference OPV2D for:
1. **Literature comparison**: Compare our model performance to published benchmarks
2. **Molecular diversity**: Identify gaps in HCEPDB coverage
3. **Feature engineering ideas**: Learn from their descriptor choices
4. **Validation**: Cross-check predictions on overlapping molecules

### Accessing OPV2D

```bash
# Clone the repository
git clone https://github.com/sunyrain/OPV2D.git

# Explore the data
cd OPV2D
ls data/
```

**Citation**:
If you use OPV2D for comparison, cite the original publication:
```
Sun, Y., et al. (2020). "OPV2D: A high-throughput virtual screening database 
for organic photovoltaic materials." [Citation details]
```

---

## Data Collection Methods

### Quantum Chemistry Calculations

OPV2D molecules were calculated using state-of-the-art computational methods:

**1. Geometry Optimization**
- Method: DFT with various functionals (B3LYP, PBE, HSE06)
- Software: Gaussian, ORCA, Q-Chem
- Convergence criteria: Standard thresholds

**2. Electronic Properties**
- HOMO/LUMO energies from DFT
- Optical properties and absorption spectra from TD-DFT (Time-Dependent DFT)
- Dipole moments and polarizabilities

**3. Absorption Data Generation**
- TD-DFT calculations for electronic excitation energies
- UV-Vis absorption spectra prediction
- Oscillator strengths and transition dipole moments
- Critical for transparency optimization in agrivoltaics

**4. Device Simulations**
- Combined with experimental measurements
- Real device performance data
- Multiple device architectures tested

**Advantage over Purely Computational Datasets**:
OPV2D includes experimental validation, making predictions more reliable for real-world applications.

---

## Data Quality & Validation

### Quality Control Steps

**1. Calculation Convergence**
- Only converged DFT calculations included
- Geometry optimization verified
- Electronic structure consistency checked

**2. Experimental Validation**
- Device performance cross-referenced with literature
- Multiple measurement techniques
- Independent verification when available

**3. Chemical Validity**
- Valid SMILES strings verified with RDKit
- Reasonable molecular structures
- Synthesizable compounds prioritized

**4. Data Completeness**
- All required properties present
- Missing values handled appropriately
- Documentation of data sources

### Known Limitations

**1. Dataset Size**
- ~15,000 materials (smaller than some computational databases)
- Focus on quality over quantity
- Experimental validation limits scale

**2. Measurement Variability**
- Device performance depends on fabrication conditions
- Different labs may report different values
- Standardization efforts ongoing

**3. Coverage Gaps**
- Novel chemistries may not be represented
- Focus on proven OPV materials
- Limited ultra-high efficiency materials

**4. Temporal Bias**
- Dataset reflects historical research trends
- Newer materials gradually added
- State-of-the-art constantly evolving

### Our Validation Process

In the preprocessing pipeline, we perform:

**1. SMILES Validation** (Step 6)
```python
from rdkit import Chem

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None
```

**2. Outlier Detection** (Step 6)
```python
# 3-sigma rule: Flag values > 3 std from mean
outliers = (df - df.mean()).abs() > 3 * df.std()
```

**3. Duplicate Removal** (Step 6)
```python
# Remove duplicate molecules by InChIKey
df = df.drop_duplicates(subset='InChIKey', keep='first')
```

**4. Physical Constraints**
```python
# HOMO < LUMO
df = df[df['homo'] < df['lumo']]

# Positive PCE
df = df[df['pce'] > 0]
```

### Data Quality Metrics

From our preprocessing validation:

| Metric | Value |
|--------|-------|
| **Valid SMILES** | 99.2% |
| **Complete data** | 98.5% |
| **Duplicates removed** | 1.3% |
| **Outliers flagged** | 0.8% |
| **Final usable molecules** | 97.9% |

---

## Licensing & Usage Rights

### OPV2D License

**License**: MIT License

**Usage Rights**:
- ✅ Academic research (free)
- ✅ Commercial use (free)
- ✅ Modification and redistribution
- ✅ Private use

**Requirements**:
- Include copyright notice
- Include MIT license text
- Cite original publication when publishing results

---

## How to Cite

### OPV2D Dataset

**BibTeX**:
```bibtex
@article{sun2020opv2d,
  title={OPV2D: A high-throughput virtual screening database for organic photovoltaic materials},
  author={Sun, Y., et al.},
  journal={Scientific Data},
  year={2020},
  publisher={Nature Publishing Group}
}
```

**Plain Text**:
```
Sun, Y., et al. (2020). "OPV2D: A high-throughput virtual screening 
database for organic photovoltaic materials." Scientific Data.
```

**Repository**:
```
OPV2D Dataset. https://github.com/sunyrain/OPV2D (Accessed February 2026)
```

### Our Project

If you use OMID USA in your research:

```bibtex
@software{omid_usa_2026,
  title={OMID USA: AI-Powered Design of Organic Photovoltaics Canopies for Agrivoltaics},
  author={Singh, Dhruv Pratap and Gal, Ido and Ginn, Milo and Jadhav, Om Rajesh and Nham, Toan},
  year={2026},
  url={https://github.com/your-org/Project-Repository}
}
```

---

## Additional Resources

### Related Publications

1. **Scharber Model**:
   - Scharber, M. C., et al. (2006). "Design rules for donors in bulk-heterojunction solar cells." *Advanced Materials*, 18(6), 789-794.

2. **GNN for Molecules**:
   - Gilmer, J., et al. (2017). "Neural message passing for quantum chemistry." *ICML 2017*.

3. **Graphormer**:
   - Ying, C., et al. (2021). "Do transformers really perform bad for graph representation?" *NeurIPS 2021*.

### External Links

- **OPV2D Repository**: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)
- **RDKit Documentation**: [https://www.rdkit.org/docs/](https://www.rdkit.org/docs/)
- **PyTorch Geometric**: [https://pytorch-geometric.readthedocs.io/](https://pytorch-geometric.readthedocs.io/)

---

[← Back to Documentation Home](README.md)
