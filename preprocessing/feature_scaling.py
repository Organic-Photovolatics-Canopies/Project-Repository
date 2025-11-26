"""
Feature and target scaling
"""

import pandas as pd
from sklearn.preprocessing import StandardScaler
from .config import QUANTUM_FEATURES, MOLECULAR_FEATURES, DERIVED_FEATURES, TARGET_COLUMNS


def scale_features(df, train_idx):
    """
    Scale features and targets using StandardScaler

    CRITICAL: Fit only on training set to prevent data leakage

    Args:
        df: Full dataframe
        train_idx: Training set indices

    Returns:
        tuple: (df_scaled, scalers_dict)
    """
    print("Scaling features and targets...")

    # Define all feature columns
    all_features = QUANTUM_FEATURES + MOLECULAR_FEATURES + DERIVED_FEATURES

    # Filter to actually present columns
    all_features = [f for f in all_features if f in df.columns]

    df_scaled = df.copy()

    # Scale features
    feature_scaler = StandardScaler()
    feature_scaler.fit(df.loc[train_idx, all_features])
    df_scaled[all_features] = feature_scaler.transform(df[all_features])

    print(f"✓ Scaled {len(all_features)} features")

    # Scale targets
    target_scaler = StandardScaler()
    # Only fit on molecules with targets in training set
    train_with_targets = df.loc[train_idx].dropna(subset=TARGET_COLUMNS)
    target_scaler.fit(train_with_targets[TARGET_COLUMNS])

    # Transform all targets
    non_null_targets = ~df[TARGET_COLUMNS].isna().any(axis=1)
    df_scaled.loc[non_null_targets, TARGET_COLUMNS] = target_scaler.transform(
        df.loc[non_null_targets, TARGET_COLUMNS]
    )

    print(f"✓ Scaled {len(TARGET_COLUMNS)} targets")

    scalers = {
        'feature': feature_scaler,
        'target': target_scaler
    }

    return df_scaled, scalers
