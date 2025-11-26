"""
Prepare multi-task regression targets
"""

import pandas as pd
from .config import TARGET_COLUMNS, HANDLE_MISSING


def prepare_targets(df, graphs=None):
    """
    Prepare PCE, Voc, Jsc targets

    Implements recommended option:
    - Multi-task regression (all three targets)
    - Drop molecules without targets

    Args:
        df: DataFrame with target columns
        graphs: Optional list of graph objects to filter alongside df

    Returns:
        pd.DataFrame or tuple: Filtered df, or (df, graphs) if graphs provided
    """
    print("Preparing multi-task targets...")

    initial_count = len(df)

    # Check for missing targets
    missing_mask = df[TARGET_COLUMNS].isna().any(axis=1)
    missing_count = missing_mask.sum()

    print(f"  Missing targets: {missing_count} / {initial_count} molecules")

    if HANDLE_MISSING == "drop":
        # Keep track of which indices to keep
        keep_mask = ~missing_mask
        df = df[keep_mask].reset_index(drop=True)

        # Filter graphs if provided
        if graphs is not None:
            graphs = [g for i, g in enumerate(graphs) if keep_mask.iloc[i]]

        print(f"✓ Dropped molecules without targets: {len(df)} remaining")
    else:
        raise NotImplementedError(f"Missing handling '{HANDLE_MISSING}' not implemented")

    # Validate target ranges
    for target in TARGET_COLUMNS:
        print(f"  {target}: [{df[target].min():.3f}, {df[target].max():.3f}]")

    # Return both df and graphs if graphs were provided
    if graphs is not None:
        return df, graphs
    else:
        return df
