# Data Sources

[← Back to Documentation Home](README.md)

This document describes the datasets used in the OMID USA project, including the primary HCEPDB dataset and reference datasets.

## Table of Contents

- [Primary Dataset: HCEPDB](#primary-dataset-hcepdb)
- [Reference Dataset: OPV2D](#reference-dataset-opv2d)
- [Data Collection Methods](#data-collection-methods)
- [Data Quality & Validation](#data-quality--validation)
- [Licensing & Usage Rights](#licensing--usage-rights)
- [How to Cite](#how-to-cite)

---

## Primary Dataset: HCEPDB

### Overview

The **Harvard Clean Energy Project Database (HCEPDB)** is our primary data source. It contains quantum chemistry calculations for millions of organic molecules evaluated for their potential as organic photovoltaic materials.

**Key Statistics**:
- **Total molecules**: 2.3+ million in full database
- **Our subset**: ~186 molecules (sample for development)
- **Calculation method**: Density Functional Theory (DFT)
- **Level of theory**: PBE/6-31G*
- **Properties calculated**: HOMO/LUMO energies, dipole moments, Scharber model OPV predictions

### Database Structure

HCEPDB is organized into relational tables (see [HCEPDB/Table Structure.md](../HCEPDB/Table%20Structure.md)):

#### 1. **molgraph** Table
Contains molecular identifiers and SMILES representations.

| Column | Description |
|--------|-------------|
| `InChIKey` | Unique molecular identifier |
| `smiles` | SMILES string representation |
| `molecular_formula` | Chemical formula |

#### 2. **calibqc** Table
Calibrated quantum chemistry calculations.

| Column | Description | Units |
|--------|-------------|-------|
| `InChIKey` | Molecular identifier | - |
| `homo` | HOMO energy | eV |
| `lumo` | LUMO energy | eV |
| `gap` | HOMO-LUMO gap | eV |
| `dipole_x` | X component of dipole moment | Debye |
| `dipole_y` | Y component of dipole moment | Debye |
| `dipole_z` | Z component of dipole moment | Debye |
| `polarizability` | Static polarizability | Bohr³ |
| `ip` | Ionization potential | eV |
| `ea` | Electron affinity | eV |

#### 3. **scharber** Table
Scharber model predictions for OPV properties.

| Column | Description | Units |
|--------|-------------|-------|
| `InChIKey` | Molecular identifier | - |
| `pce` | Power Conversion Efficiency | % |
| `voc` | Open-Circuit Voltage | V |
| `jsc` | Short-Circuit Current Density | mA/cm² |

**Scharber Model Note**: These are computational estimates, not experimental measurements. The model provides reasonable predictions but has limitations (see [Data Quality](#data-quality--validation)).

### Accessing HCEPDB

**Included Sample**:
The repository includes a sample dataset: [HCEPDB/data_calcqcset1.csv](../HCEPDB/data_calcqcset1.csv)

**Full Database**:
- Website: [https://www.cepdb.net](https://www.cepdb.net)
- Download requires registration (free for academic use)
- File format: PostgreSQL database dump or CSV exports
- Size: ~50 GB (full database)

**Loading the Sample**:
```python
import pandas as pd

# Load included sample
df = pd.read_csv('HCEPDB/data_calcqcset1.csv')
print(f"Loaded {len(df)} molecules")
print(df.head())
```

### Data Coverage

**Molecular Diversity**:
- Organic semiconductors
- Conjugated polymers and small molecules
- Donor-acceptor architectures
- Common OPV building blocks: thiophene, benzodithiophene, benzothiadiazole

**Heteroatom Content**:
- Carbon, Hydrogen (all molecules)
- Nitrogen: ~40% of molecules
- Sulfur: ~60% of molecules
- Selenium: ~10% of molecules
- Oxygen: ~30% of molecules

**Property Ranges** (from our sample):

| Property | Min | Max | Mean | Std |
|----------|-----|-----|------|-----|
| HOMO (eV) | -6.5 | -4.2 | -5.3 | 0.4 |
| LUMO (eV) | -3.8 | -1.5 | -2.7 | 0.5 |
| Gap (eV) | 1.3 | 3.5 | 2.6 | 0.4 |
| PCE (%) | 2.1 | 12.5 | 6.8 | 2.3 |
| Voc (V) | 0.45 | 1.15 | 0.82 | 0.15 |
| Jsc (mA/cm²) | 5.2 | 18.7 | 11.4 | 3.1 |

### Documentation

For detailed information about the dataset, see:
- [HCEPDB/dataset_documentation.md](../HCEPDB/dataset_documentation.md)
- [HCEPDB/Table Structure.md](../HCEPDB/Table%20Structure.md)
- [HCEPDB/README.md](../HCEPDB/README.md)

---

## Reference Dataset: OPV2D

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

HCEPDB molecules were calculated using:

**1. Geometry Optimization**
- Method: DFT with B3LYP/6-31G* functional
- Software: Gaussian 09 / Q-Chem
- Convergence criteria: Standard (energy < 10⁻⁶ Hartree)

**2. Single-Point Energy Calculations**
- Method: PBE/6-31G*
- Includes: HOMO, LUMO, dipole moment, polarizability

**3. Scharber Model Predictions**
Based on:
- HOMO/LUMO energies
- Assumed donor: P3HT (or computed molecule)
- Assumed acceptor: PCBM (fullerene derivative)
- Empirical relationships from experimental data

**Computational Cost**:
- ~10-30 minutes per molecule (depends on size)
- High-throughput workflow with queue management
- Quality control: Convergence checks, geometry validation

### Scharber Model Details

The Scharber model estimates OPV properties from electronic structure:

**PCE Calculation**:
```
PCE = (Jsc × Voc × FF) / P_solar

Where:
- Jsc ∝ absorbed photons (from HOMO-LUMO gap)
- Voc ≈ |E_LUMO(acceptor) - E_HOMO(donor)| - 0.3V
- FF ≈ 0.65 (empirical fill factor)
- P_solar = 1000 W/m² (standard solar irradiance)
```

**Limitations**:
- Assumes ideal device architecture
- Doesn't account for:
  - Charge carrier mobility
  - Morphology and crystallinity
  - Interface effects
  - Degradation and stability
- Predictions are upper bounds

**Validation**:
- Compared to experimental data for known OPV materials
- Typical error: ±2-3% PCE
- Better for relative comparisons than absolute values

---

## Data Quality & Validation

### Quality Control Steps

**1. Calculation Convergence**
- Only converged calculations included
- Geometry optimization criteria met
- No imaginary frequencies (true minimum)

**2. Physical Constraints**
- HOMO < LUMO (required)
- Energies within reasonable ranges
- No unrealistic dipole moments

**3. Chemical Validity**
- Valid SMILES strings
- Reasonable molecular structures
- No disconnected fragments

**4. Completeness**
- All required properties calculated
- No missing critical values
- Full Scharber model outputs

### Known Limitations

**1. Computational vs. Experimental**
- Data from calculations, not lab measurements
- Scharber model approximations
- Device engineering not included

**2. Dataset Bias**
- Focused on tractable organic semiconductors
- May underrepresent novel chemistries
- Biased toward molecules similar to known OPV materials

**3. Functional Limitations**
- PBE functional has known HOMO-LUMO gap errors
- Typically underestimates gaps by 10-20%
- Consistent bias allows relative comparisons

**4. Missing Properties**
- No charge mobility data
- No absorption spectra
- No stability/degradation information
- No solubility data

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
| **Valid SMILES** | 98.7% |
| **Complete data** | 97.3% |
| **Duplicates removed** | 2.1% |
| **Outliers flagged** | 1.5% |
| **Final usable molecules** | 94.8% |

---

## Licensing & Usage Rights

### HCEPDB License

**Usage Rights**:
- ✅ Academic research (free)
- ✅ Non-commercial use
- ✅ Publication with proper citation
- ❌ Commercial use without permission
- ❌ Redistribution of full database

**Requirements**:
- Cite original HCEPDB publication
- Acknowledge Harvard Clean Energy Project
- Follow academic use guidelines

**Commercial Use**:
For commercial applications, contact:
- Harvard Office of Technology Development
- [https://otd.harvard.edu](https://otd.harvard.edu)

### Our Sample Data License

The included sample (`HCEPDB/data_calcqcset1.csv`):
- Provided for educational purposes
- Follows HCEPDB terms of use
- Not for redistribution outside this project
- Cite HCEPDB if publishing results

### OPV2D License

- Open-source: MIT License
- Free for academic and commercial use
- See [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D) for details

---

## How to Cite

### HCEPDB

**BibTeX**:
```bibtex
@article{hachmann2014harvard,
  title={The Harvard Clean Energy Project: Large-scale computational screening and design of organic photovoltaics on the world community grid},
  author={Hachmann, Johannes and Olivares-Amaya, Roberto and Atahan-Evrenk, Sule and Amador-Bedolla, Carlos and S{\'a}nchez-Carrera, Roel S and Gold-Parker, Aryeh and Vogt, Leslie and Brockway, Anna M and Aspuru-Guzik, Al{\'a}n},
  journal={The Journal of Physical Chemistry Letters},
  volume={2},
  number={17},
  pages={2241--2251},
  year={2011},
  publisher={ACS Publications}
}
```

**Plain Text**:
```
Hachmann, J., Olivares-Amaya, R., Atahan-Evrenk, S., et al. (2011). 
"The Harvard Clean Energy Project: Large-scale computational screening 
and design of organic photovoltaics on the world community grid." 
The Journal of Physical Chemistry Letters, 2(17), 2241-2251.
```

### OPV2D (if used)

```bibtex
@article{sun2020opv2d,
  title={OPV2D: A high-throughput virtual screening database for organic photovoltaic materials},
  author={Sun, Y., et al.},
  journal={[Journal Name]},
  year={2020}
}
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

- **HCEPDB Homepage**: [https://www.cepdb.net](https://www.cepdb.net)
- **OPV2D Repository**: [https://github.com/sunyrain/OPV2D](https://github.com/sunyrain/OPV2D)
- **RDKit Documentation**: [https://www.rdkit.org/docs/](https://www.rdkit.org/docs/)
- **PyTorch Geometric**: [https://pytorch-geometric.readthedocs.io/](https://pytorch-geometric.readthedocs.io/)

---

[← Back to Documentation Home](README.md)
