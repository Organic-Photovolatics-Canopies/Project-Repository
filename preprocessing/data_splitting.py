"""
Data splitting using scaffold-based strategy
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from .config import SPLIT_RATIOS, RANDOM_SEED


def scaffold_split(df, smiles_col='SMILES_str', seed=RANDOM_SEED):
    """
    Scaffold-based train/val/test split

    Ensures molecules with same core structure are in same split
    Tests true generalization to new scaffolds

    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column
        seed: Random seed for reproducibility

    Returns:
        tuple: (train_idx, val_idx, test_idx)
    """
    print(f"Performing scaffold split (train={SPLIT_RATIOS['train']}, "
          f"val={SPLIT_RATIOS['val']}, test={SPLIT_RATIOS['test']})...")

    np.random.seed(seed)

    # Generate scaffolds
    scaffolds = {}
    for idx, smiles in enumerate(df[smiles_col]):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Assign orphan scaffold
            scaffold = f"invalid_{idx}"
        else:
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)

        if scaffold not in scaffolds:
            scaffolds[scaffold] = []
        scaffolds[scaffold].append(idx)

    # Sort scaffolds by size (largest first)
    scaffold_sets = list(scaffolds.values())
    scaffold_sets.sort(key=len, reverse=True)

    # Greedily assign scaffolds to splits
    train_idx, val_idx, test_idx = [], [], []
    train_size, val_size, test_size = 0, 0, 0
    total_size = len(df)

    train_target = SPLIT_RATIOS['train'] * total_size
    val_target = SPLIT_RATIOS['val'] * total_size

    for scaffold_set in scaffold_sets:
        if train_size < train_target:
            train_idx.extend(scaffold_set)
            train_size += len(scaffold_set)
        elif val_size < val_target:
            val_idx.extend(scaffold_set)
            val_size += len(scaffold_set)
        else:
            test_idx.extend(scaffold_set)
            test_size += len(scaffold_set)

    print(f"✓ Split complete:")
    print(f"  - Train: {len(train_idx)} molecules ({len(train_idx)/total_size*100:.1f}%)")
    print(f"  - Val:   {len(val_idx)} molecules ({len(val_idx)/total_size*100:.1f}%)")
    print(f"  - Test:  {len(test_idx)} molecules ({len(test_idx)/total_size*100:.1f}%)")
    print(f"  - Unique scaffolds: {len(scaffolds)}")

    return train_idx, val_idx, test_idx
