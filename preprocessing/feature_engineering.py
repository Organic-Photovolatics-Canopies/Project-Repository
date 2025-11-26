"""
Engineer additional features from molecular and quantum data
"""

import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def parse_stoichiometry(stoich_str):
    """
    Extract element counts from molecular formula

    Args:
        stoich_str: Molecular formula (e.g., "C21H11NOSSe")

    Returns:
        dict: Element counts (e.g., {'C': 21, 'H': 11, 'N': 1, ...})
    """
    pattern = r'([A-Z][a-z]?)(\d*)'
    elements = {}
    for match in re.finditer(pattern, stoich_str):
        elem = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        elements[elem] = count
    return elements


def engineer_features(df):
    """
    Add derived features to dataframe

    New columns:
        - C_count, H_count, heteroatom_count, heteroatom_ratio
        - electron_affinity_proxy, ionization_potential_proxy
        - aromatic_rings, rotatable_bonds
        - hbd, hba (H-bond donors/acceptors)
        - conjugation_estimate

    Args:
        df: DataFrame with molecular and quantum features

    Returns:
        pd.DataFrame: DataFrame with additional columns
    """
    print("Engineering derived features...")

    df = df.copy()

    # Convert all numeric columns to proper types (in case they're loaded as strings)
    numeric_cols = ['mass', 'n_heavyatoms', 'n_el', 'e_homo_alpha', 'e_lumo_alpha',
                    'e_gap_alpha', 'e_gap_min', 'pce', 'voc', 'jsc']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Stoichiometry features
    df['C_count'] = df['stoich_str'].apply(lambda x: parse_stoichiometry(x).get('C', 0))
    df['H_count'] = df['stoich_str'].apply(lambda x: parse_stoichiometry(x).get('H', 0))
    df['heteroatom_count'] = df['stoich_str'].apply(
        lambda x: sum(v for k, v in parse_stoichiometry(x).items() if k not in ['C', 'H'])
    )
    df['heteroatom_ratio'] = df['heteroatom_count'] / df['n_heavyatoms'].replace(0, 1)

    # Electronic properties
    df['electron_affinity_proxy'] = -df['e_lumo_alpha']
    df['ionization_potential_proxy'] = -df['e_homo_alpha']

    # RDKit descriptors
    def safe_rdkit_calc(smiles, func):
        """Safely calculate RDKit descriptor"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return func(mol)
        except:
            return None

    df['aromatic_rings'] = df['SMILES_str'].apply(
        lambda s: safe_rdkit_calc(s, Descriptors.NumAromaticRings)
    )
    df['rotatable_bonds'] = df['SMILES_str'].apply(
        lambda s: safe_rdkit_calc(s, Descriptors.NumRotatableBonds)
    )
    df['hbd'] = df['SMILES_str'].apply(
        lambda s: safe_rdkit_calc(s, Descriptors.NumHDonors)
    )
    df['hba'] = df['SMILES_str'].apply(
        lambda s: safe_rdkit_calc(s, Descriptors.NumHAcceptors)
    )

    # Rough conjugation estimate
    df['conjugation_estimate'] = df['SMILES_str'].apply(
        lambda s: s.count('=') + s.count(':')
    )

    print(f"✓ Added 11 derived features")

    return df
