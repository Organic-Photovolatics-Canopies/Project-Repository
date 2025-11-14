# Dataset Construction Pipeline for GNN Transparency Prediction

This document explains how the final machine‑learning dataset
(`v_gnn_dataset`) used for the GNN was produced. The process transforms
multiple raw quantum‑chemical and Scharber-model tables into a cleaned,
physically consistent, ML‑ready dataset.

------------------------------------------------------------------------

## 1. Source Data Tables

Data originates from four core tables:

-   **`data_molgraph`** -- molecular identity (SMILES, atom count, mass)
-   **`data_molgeom`** -- molecular geometries\
-   **`data_calcqcset1`** -- DFT electronic structure (HOMO/LUMO
    energies, band gaps, dipole moment)
-   **`data_scharber`** -- photovoltaic metrics (PCE, Jsc)

These contain multiple experiment entries per molecule.

------------------------------------------------------------------------

## 2. Cleaning Quantum Calculations

Repeated DFT runs per molecule can vary. We clean by **statistical
filtering based on variance**.

### 2.1 Variance thresholds

    variance(e_gap_min)   ≤ 0.03 eV
    variance(e_gap_alpha) ≤ 0.03 eV
    variance(e_gap_beta)  ≤ 0.03 eV
    variance(dipole)      ≤ 0.15 Debye²

Chosen based on typical DFT noise across conformers and SCF stability.

### 2.2 Creating `data_calcqcset1_clean`

Steps:

1.  Group calculations by `mol_graph_id`.
2.  Compute `AVG()` and `VAR_POP()`.
3.  Keep molecules with stable values.
4.  Store averaged fields (e.g., `e_gap_min_avg`, `dipmom_total_avg`).

------------------------------------------------------------------------

## 3. Cleaning Calibrated QC Data

`data_calibqcset1` often contains multiple calibration runs.

We repeat the variance-based filtering:

-   Group by `mol_graph_id` and `calib_type`.
-   Apply same thresholds.
-   Produce the table `data_calibqcset1_clean`.

------------------------------------------------------------------------

## 4. Aggregating Scharber Metrics

Some molecules have several Scharber entries.

We compute:

-   `jsc_avg`
-   `pce_avg`

These appear in `v_scharber_agg`.

------------------------------------------------------------------------

## 5. Building Absorption Features

A unified view `v_absorption_features` combines:

-   Cleaned DFT averages (`data_calcqcset1_clean`)
-   Aggregated Scharber metrics (`v_scharber_agg`)
-   Molecule metadata (`data_molgraph`)

### Derived key features

  Feature                      Meaning
  ---------------------------- -----------------------------------------
  `gap_asymmetry`              
  `absorption_strength_norm`   normalized Jsc (absorption)
  `onset_norm`                 normalized bandgap (transparency proxy)

### Filters applied

-   Molecule must appear in Scharber data.
-   `0 ≤ pce_avg ≤ 40`
-   `0 ≤ jsc_avg ≤ 30`
-   `gap_asymmetry < 0.2 eV` to exclude SCF-broken molecules.

------------------------------------------------------------------------

## 6. Creating the Final Label

We classify molecules as **transparent ≥ 50%** using the
physics‑informed proxy `onset_norm`:

    y_transparent = 1 if onset_norm >= 0.5 else 0

This threshold cleanly separates high‑gap (transparent) and low‑gap
(absorbing) systems.

------------------------------------------------------------------------

## 7. Final Machine‑Learning Dataset

The final view `v_gnn_dataset` includes:

-   SMILES (used to build graph structure)
-   Cleaned electronic structure features
-   Scharber-based absorption features
-   Normalized features (onset, absorption)
-   **Binary transparency label**

### Perfect for GNNs:

-   **Graph features** → from SMILES\
-   **Global numeric features** → band gaps, dipole, Jsc, PCE\
-   **Target** → binary transparency class

------------------------------------------------------------------------

## 8. Pipeline Summary

    data_calcqcset1        → cleaned → data_calcqcset1_clean
    data_calibqcset1       → cleaned → data_calibqcset1_clean
    data_scharber          → aggregated → v_scharber_agg
    molgraph + clean DFT + Scharber → v_absorption_features
                                             ↓
               transparency label (onset_norm ≥ 0.5)
                                             ↓
                              v_gnn_dataset

This is the dataset used for training the transparency prediction GNN.
