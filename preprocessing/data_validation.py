"""
Data validation and cleaning
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from scipy import stats
from .config import QUANTUM_FEATURES


def validate_smiles(df, smiles_col='SMILES_str'):
    """
    Check for invalid SMILES strings

    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column

    Returns:
        pd.DataFrame: DataFrame with 'valid_smiles' column added
    """
    print("Validating SMILES...")

    df['valid_smiles'] = df[smiles_col].apply(
        lambda s: Chem.MolFromSmiles(s) is not None
    )

    invalid = df[~df['valid_smiles']]
    if len(invalid) > 0:
        print(f"  WARNING: {len(invalid)} invalid SMILES found")
        print(invalid[['id', smiles_col]])
    else:
        print(f"✓ All {len(df)} SMILES valid")

    return df


def detect_outliers(df, features, threshold=3):
    """
    Detect outliers using z-score method

    Args:
        df: DataFrame
        features: List of feature column names
        threshold: Z-score threshold

    Returns:
        dict: Outlier counts per feature
    """
    print(f"Detecting outliers (z-score > {threshold})...")

    outlier_counts = {}
    for feature in features:
        if feature in df.columns:
            z_scores = np.abs(stats.zscore(df[feature].dropna()))
            outlier_counts[feature] = (z_scores > threshold).sum()

    total_outliers = sum(outlier_counts.values())
    if total_outliers > 0:
        print(f"  Found outliers in {len(outlier_counts)} features:")
        for feat, count in outlier_counts.items():
            if count > 0:
                print(f"    - {feat}: {count} outliers")
    else:
        print(f"✓ No significant outliers detected")

    return outlier_counts


def remove_duplicates(df, smiles_col='SMILES_str'):
    """
    Remove duplicate molecules by canonical SMILES

    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column

    Returns:
        pd.DataFrame: DataFrame with duplicates removed
    """
    print("Checking for duplicates...")

    def canonicalize_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)

    df['canonical_smiles'] = df[smiles_col].apply(canonicalize_smiles)

    initial_count = len(df)
    df = df.drop_duplicates(subset=['canonical_smiles'], keep='first')

    removed = initial_count - len(df)
    if removed > 0:
        print(f"  Removed {removed} duplicate molecules")
    else:
        print(f"✓ No duplicates found")

    return df


def validate_and_clean(df, graphs):
    """
    Run all validation checks

    Args:
        df: DataFrame with molecular data
        graphs: List of PyG graph objects

    Returns:
        tuple: (cleaned_df, cleaned_graphs)
    """
    print("\n" + "="*60)
    print("DATA VALIDATION & CLEANING")
    print("="*60)

    # Validate SMILES
    df = validate_smiles(df)

    # Remove invalid SMILES
    valid_mask = df['valid_smiles']
    df = df[valid_mask].reset_index(drop=True)
    graphs = [g for i, g in enumerate(graphs) if valid_mask.iloc[i]]

    # Detect outliers (informational only)
    detect_outliers(df, QUANTUM_FEATURES)

    # Remove duplicates
    df_before_dedup = df.copy()
    df = remove_duplicates(df)

    # Update graphs to match deduplicated df
    if len(df) < len(df_before_dedup):
        # Find indices that were kept
        kept_indices = df.index.tolist()
        graphs = [graphs[i] for i in kept_indices]

    df = df.reset_index(drop=True)

    print(f"\n✓ Validation complete: {len(df)} clean molecules")

    return df, graphs
