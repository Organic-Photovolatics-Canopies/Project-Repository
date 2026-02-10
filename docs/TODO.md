# TODO: Documentation Updates Needed

This file tracks documentation updates that need to be made as the project develops.

## Priority Updates

### 1. OPV2D Preprocessing Notebook
**Status**: Pending
**What's needed**:
- Jupyter notebook showing how to preprocess OPV2D dataset
- Step-by-step data loading from OPV2D repository
- Feature engineering specific to OPV2D data format
- Example outputs and visualizations

**Location**: Should be added to:
- `preprocessing/opv2d_preprocessing.ipynb` (main notebook)
- Referenced in `docs/user-guide.md` (Workflow 1)
- Add link in `docs/README.md`

### 2. TD-DFT Absorption Data Generation
**Status**: In Progress
**What's needed**:
- Run TD-DFT calculations on OPV2D molecules
- Generate UV-Vis absorption spectra
- Calculate oscillator strengths and transition dipole moments
- Add absorption wavelength features to preprocessing pipeline
- Document TD-DFT methodology in user manual

**Applications**:
- Transparency optimization for agrivoltaics
- Selective absorption targeting (UV/near-IR)
- Design molecules with specific transparency ranges

**Technical Notes**:
- Method: TD-DFT with range-separated hybrid functionals (CAM-B3LYP, ωB97XD)
- Solvent effects: PCM (Polarizable Continuum Model)
- Basis set: 6-31G(d) or larger for accuracy
- Number of excited states: 10-20 lowest singlet states

### 3. Documentation Images/Screenshots
**Status**: Pending
**What's needed**:

**For `docs/images/`**:
- `pipeline-output-example.png` - Terminal screenshot of successful pipeline run
- `output-directory-structure.png` - File tree showing preprocessed_data/ contents  
- `molecular-graph-visualization.png` - Example molecule as graph
- `training-loss-curve.png` - Loss vs epoch plot
- `property-distributions.png` - Histograms of PCE, Voc, Jsc
- `web-demo-interface.png` - Screenshot of HuggingFace/Kaggle demo (when available)

**Usage**:
- Update `docs/images/README.md` with actual images
- Reference in:
  - `docs/installation.md` - successful setup screenshots
  - `docs/user-guide.md` - workflow visuals
  - `docs/user-manual.md` - architecture diagrams
  - `docs/data-sources.md` - data distribution plots

### 3. Configuration File Updates
**Status**: Pending
**What's needed**:
- Update `preprocessing/config.py` with correct OPV2D paths
- Example: `DATA_PATH = "OPV2D/data/opv_data.csv"`
- Update column name mappings if OPV2D uses different names

### 4. Data Integration Updates
**Status**: Pending
**What's needed**:
- Update `preprocessing/data_integration.py` for OPV2D schema
- Remove HCEPDB-specific table merging (calibqc, scharber, molgraph)
- Add OPV2D-specific data processing
- Handle experimental vs computational data appropriately

### 5. README Updates in Other Directories
**Status**: Pending
**Check and update**:
- `preprocessing/SETUP_GUIDE.md` - Update data source references
- `preprocessing/README.md` - If it exists, update dataset info
- Root `README.md` - Remove any remaining HCEPDB links (check team member contributions section)

## Completed Updates ✅

- [x] Primary documentation files updated (docs/README.md, data-sources.md, etc.)
- [x] Changed primary dataset from HCEPDB to OPV2D
- [x] Updated installation guide
- [x] Updated user guide workflows
- [x] Updated user manual technical sections
- [x] Updated FAQ
- [x] Updated troubleshooting guide
- [x] Root README.md documentation section updated

## How to Use This File

When you complete an update:
1. Move item from "Priority Updates" to "Completed Updates"
2. Add checkmark [x] and date
3. Commit changes with descriptive message

Example:
```bash
git add docs/ preprocessing/
git commit -m "Add OPV2D preprocessing notebook and screenshots"
```

## Notes for Team

- **Milo** (Front-End): Web demo screenshots needed when ready
- **Ido** (Data Processing): OPV2D preprocessing notebook priority
- **Om** (GNN Development): Training curve visualizations
- **Dhruv** (LLM): API documentation when integrated
- **Toan** (Model Training): Model performance plots

---

**Last Updated**: February 10, 2026
**Next Review**: When OPV2D preprocessing is implemented
