"""
Integrate relational tables into unified molecule-level dataset
"""

import pandas as pd
from .config import QC_METHOD


def integrate_data(raw_data):
    """
    Merge molgraph, calibqc, and scharber tables

    Implements:
    - Step 1: Data integration
    - Step 2: QC method selection (calibrated values only)

    Args:
        raw_data: dict from data_loader.load_raw_data()

    Returns:
        pd.DataFrame: Unified dataset with one row per molecule
    """
    molgraph = raw_data['molgraph']
    calibqc = raw_data['calibqc']
    scharber = raw_data['scharber']

    # Step 2: Filter for calibrated preliminary values
    if QC_METHOD == "calibrated_prelim":
        calibqc_filtered = calibqc[calibqc['calib_type'] == 'prelim'].copy()
        print(f"Filtered to {len(calibqc_filtered)} calibrated 'prelim' values")
    else:
        raise NotImplementedError(f"QC method '{QC_METHOD}' not implemented")

    # Step 1: Merge tables
    # molgraph ← calibqc (left join to keep all molecules)
    df = molgraph.merge(
        calibqc_filtered,
        left_on='id',
        right_on='mol_graph_id',
        how='left',
        suffixes=('', '_calib')
    )

    # Handle multiple Scharber types by averaging
    # Group by mol_graph_id and average pce, voc, jsc
    scharber_avg = scharber.groupby('mol_graph_id').agg({
        'pce': 'mean',
        'voc': 'mean',
        'jsc': 'mean'
    }).reset_index()

    # Merge Scharber predictions
    df = df.merge(
        scharber_avg,
        left_on='id',
        right_on='mol_graph_id',
        how='left',
        suffixes=('', '_scharber')
    )

    # Flag molecules with targets
    df['has_target'] = ~df['pce'].isna()

    print(f"✓ Integrated dataset: {len(df)} molecules")
    print(f"  - With targets: {df['has_target'].sum()}")
    print(f"  - Without targets: {(~df['has_target']).sum()}")

    return df
