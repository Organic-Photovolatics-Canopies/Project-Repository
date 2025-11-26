"""
Load CEPDB CSV files with column headers
"""

import pandas as pd
from pathlib import Path
from .config import (
    DATA_DIR, MOLGRAPH_COLS, CALCQC_COLS, CALIBQC_COLS, SCHARBER_COLS
)


def load_raw_data(data_dir=DATA_DIR):
    """
    Load all CEPDB CSV files with proper column names

    Args:
        data_dir: Path to CEPDB data directory

    Returns:
        dict: {'molgraph': df, 'calcqc': df, 'calibqc': df, 'scharber': df}
    """
    data_dir = Path(data_dir)

    print("Loading raw CSV files...")

    # Load with schema
    molgraph = pd.read_csv(
        data_dir / 'data_molgraph.csv',
        header=None,
        names=MOLGRAPH_COLS
    )

    calcqc = pd.read_csv(
        data_dir / 'data_calcqcset1.csv',
        header=None,
        names=CALCQC_COLS
    )

    calibqc = pd.read_csv(
        data_dir / 'data_calibqcset1.csv',
        header=None,
        names=CALIBQC_COLS
    )

    scharber = pd.read_csv(
        data_dir / 'data_scharber.csv',
        header=None,
        names=SCHARBER_COLS
    )

    print(f"✓ Loaded {len(molgraph)} molecules")
    print(f"✓ Loaded {len(calcqc)} QC calculations")
    print(f"✓ Loaded {len(calibqc)} calibrated QC values")
    print(f"✓ Loaded {len(scharber)} Scharber predictions")

    return {
        'molgraph': molgraph,
        'calcqc': calcqc,
        'calibqc': calibqc,
        'scharber': scharber
    }
