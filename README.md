# AI-Powered OPV Design

## Fall Design Report

1. Team Members
   - **Advisor:** Dr. Mohsen Rezayat
   - Milo Ginn
   - Om Jadhav
   - Dhruv Singh
   - Ido Gal
   - Toan Nham
2. [Project Description](./Project-Description.md)
3. User Stories and Design Diagrams
   - [User Stories](./User%20Stories.md)
   - [Design Diagrams](./Design%20Diagrams/README.md)
4. Project Tasks and Timeline
   - [Task List](./Task%20List.md)
   - [Timeline](./Milestones,%20Timeline,%20and%20Effort%20Matrix.pdf)
5. [ABET Concerns Essay](./Project%20Constraints%20Essay.pdf)
6. [PPT Presentation](./Fall%20Design%20Presentation.pptx)
7. Self-Assessment Essays
   - [Milo Ginn](./Self%20Assessment%20Essays/Milo%20Ginn.md)
8. Professional Biographies
   - [Dhruv Pratap Singh](./Member%20Bios/Dhruv%20Pratap%20Signh.md)
   - [Milo Ginn](./Member%20Bios/Milo%20Ginn.md)
   - [Om Rajesh Jadhav](./Member%20Bios/Om%20Rajesh%20Jadhav.md)
   - [Ido Gal](./Member%20Bios/Ido%20Gal.md)
   - [Toan Nham](./Member%20Bios/Toan%20Nham.md)
9. Budget
   - This project has not incurred any expenses so far.
10. Appendix
    - [Harvard Clean Energy Project Research](./HCEPDB/README.md)

This section provides references, citations, links to code repositories, meeting
notes, and evidence of work for each team member.

**Team Meetings:**

- The CS team met twice weekly for a total of 2.5 hours per week, focusing on
  technical development, data analysis, and code review.
- The full project team also met weekly to coordinate interdisciplinary tasks
  and project milestones.

**Evidence of Effort (45+ hours per member):**

Each member will add their own contributions individually below.

- **Ido Gal:**
  - Work in progress is documented in the `progress/` folder, including:
    - Quantum chemistry dataset analysis
      ([Quantum_Chemistry_Analysis_Report.md](progress/Quantum_Chemistry_Analysis_Report.md))
    - Data cleaning and pipeline documentation
      ([dataset_documentation.md](progress/dataset_documentation.md),
      [pipeline.sql](progress/pipeline.sql))
    - Data files and analysis images (e.g., `data_calcqcset1.csv`,
      `comprehensive_analysis.png`, `database_structure_analysis.png`)
  - Researched the Harvard Clean Energy dataset, performed data cleaning, and
    prepared the data for machine learning model training. This involved:
    - Downloading, exploring, and understanding the dataset structure
    - Cleaning quantum chemistry tables using statistical filtering (see
      `dataset_documentation.md`)
    - Writing SQL and Python scripts to process and validate the data (see
      `pipeline.sql`)
    - Documenting the process and results in progress reports
- **Dhruv:**
  - Work in progress is documented in the `pyg-exploration/` folder, including:
    - PyTorch Geometric implementation and baseline analysis
      (pyg_tutorial_complete.py, outputs.pdf)
    - Baseline model visualizations and performance metrics
      (baseline_predictions.png, feature_importance.png,
      transparency_distribution.png)
    - Data files and preprocessed datasets (e.g.,
      transparent_opv_candidates.csv, compound_quality_metrics.csv,
      compound_smiles_sample.csv)
  - Researched PyTorch Geometric for molecular property prediction, performed
    baseline analysis, and prepared the graph neural network implementation for
    the OPV transparency prediction task. This involved:
    - Researching graph neural network architectures for molecular property
      prediction and understanding GCN layers, message passing, and molecular
      graph representations (see pyg_tutorial_complete.py)
    - Setting up development environment with PyTorch, PyTorch Geometric, and
      RDKit chemistry library
    - Implementing baseline model achieving R² = 0.4335, RMSE = 0.7259 using
      spectral features
    - Writing 300+ lines of code for SMILES-to-graph conversion, GCN
      architecture, and training loops (see pyg_tutorial_complete.py)
    - Documenting the process and results in technical reports

- **Om:**
  - Work in progress is documented in the `Research/` folder. Including:

    - A comprehensive literature review on graph neural network (GNN) pretraining, transformer-based molecular generation, reinforcement learning for de novo design, and methodological pipelines for OPV material discovery (`litReview.tex`).
    - A structured analysis of Qiu et al.’s methodology, along with recommendations for adapting the pipeline to transparent OPV optimization (see `litReview.tex`).
    - Skeleton code and an executable test script to benchmark Hugging Face GNN models (e.g., Graphormer)(`hf_gnn_test.py`).
  - Researched modern GNN architectures, self-supervised molecular representation learning, cross-attention donor–acceptor fusion models, and transformer-based molecular generation pipelines. This involved:

    - Reviewing ~15 papers on GNN pretraining, spectral/electronic descriptors, and transformer+RL molecular design, focusing on how these techniques apply to OPV PCE and transparency.
    - Writing a 4-page IEEE-style literature review analyzing strengths, limitations, and methodological implications of current OPV ML pipelines (see `litReview.tex`).
    - Identifying gaps in existing OPV datasets, including the scarcity of AVT/spectral data, the heterogeneity of PCE measurements, and the need for DFT-calibrated optical descriptors.
    - Extracting actionable recommendations for transparency optimization, including:

      - integrating spectral/optical property pretraining targets,
      - adopting geometry-aware GNNs,
      - multi-objective prediction of PCE and AVT,
      - incorporating synthetic accessibility and uncertainty-aware scoring.
  - Developed a runnable baseline environment for evaluating Hugging Face GNNs on molecular graph inputs. This included:
    - Writing skeleton code to load, configure, and run Graphormer-style models from Hugging Face using PyTorch (`hf_gnn_test.py`).
    

**References and Citations:**

- PyTorch Geometric Documentation: https://pytorch-geometric.readthedocs.io
- Baseline Model Analysis: See pyg-exploration/outputs.pdf
- PyG Implementation Pipeline: See pyg-exploration/pyg_tutorial_complete.py
- Code Repository:
  [Project GitHub](https://github.com/Organic-Photovoltaics-Canopies/Project-Repository)

**References and Citations:**

- PyTorch Geometric Documentation:
  [https://pytorch-geometric.readthedocs.io](https://pytorch-geometric.readthedocs.io)
- Baseline Model Analysis: See `pyg-exploration/outputs.pdf`
- PyG Implementation Pipeline: See `pyg-exploration/pyg_tutorial_complete.py`
- Code Repository:
  [Project GitHub](https://github.com/Organic-Photovoltaics-Canopies/Project-Repository)

**References and Citations:**

- Harvard Clean Energy Project Dataset:
  [https://www.census.gov/data/datasets/clean-energy.html](https://www.census.gov/data/datasets/clean-energy.html)
- Quantum Chemistry Data Analysis: See
  `progress/Quantum_Chemistry_Analysis_Report.md`
- Data Cleaning Pipeline: See `progress/dataset_documentation.md` and
  `progress/pipeline.sql`
- Code Repository:
  [Project GitHub](https://github.com/Organic-Photovolatics-Canopies/Project-Repository)
