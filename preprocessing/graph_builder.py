"""
Generate molecular graphs from SMILES strings
"""

import torch
from torch_geometric.data import Data
from rdkit import Chem
from .config import INCLUDE_EXPLICIT_H


def smiles_to_graph(smiles, include_explicit_h=INCLUDE_EXPLICIT_H):
    """
    Convert SMILES to PyTorch Geometric graph

    Node features (7-dim):
        - Atomic number
        - Total degree
        - Formal charge
        - Hybridization (as int)
        - Aromaticity (0/1)
        - Total H count
        - In ring (0/1)

    Edge features (3-dim):
        - Bond type (1=single, 2=double, 3=triple, 1.5=aromatic)
        - Conjugation (0/1)
        - In ring (0/1)

    Args:
        smiles: SMILES string
        include_explicit_h: Whether to include explicit hydrogen atoms

    Returns:
        torch_geometric.data.Data or None if invalid SMILES
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    if include_explicit_h:
        mol = Chem.AddHs(mol)

    # Extract node features
    atom_features = []
    for atom in mol.GetAtoms():
        features = [
            atom.GetAtomicNum(),
            atom.GetTotalDegree(),
            atom.GetFormalCharge(),
            int(atom.GetHybridization()),  # SP=2, SP2=3, SP3=4
            int(atom.GetIsAromatic()),
            atom.GetTotalNumHs(),
            int(atom.IsInRing())
        ]
        atom_features.append(features)

    # Extract edge features
    edge_index = []
    edge_features = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        bond_type = bond.GetBondTypeAsDouble()
        is_conjugated = int(bond.GetIsConjugated())
        in_ring = int(bond.IsInRing())

        # Undirected graph: add both directions
        edge_index.extend([[i, j], [j, i]])
        edge_features.extend([[bond_type, is_conjugated, in_ring]] * 2)

    # Convert to tensors
    x = torch.tensor(atom_features, dtype=torch.float)
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous() if edge_index else torch.empty((2, 0), dtype=torch.long)
    edge_attr = torch.tensor(edge_features, dtype=torch.float) if edge_features else torch.empty((0, 3), dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)


def build_all_graphs(df, smiles_col='SMILES_str'):
    """
    Build graphs for all molecules in dataframe

    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column

    Returns:
        list: PyG Data objects (None for invalid SMILES)
    """
    print(f"Building molecular graphs from {len(df)} SMILES...")

    graphs = []
    invalid_count = 0

    for idx, smiles in enumerate(df[smiles_col]):
        graph = smiles_to_graph(smiles)
        if graph is None:
            invalid_count += 1
            print(f"  WARNING: Invalid SMILES at index {idx}: {smiles}")
        graphs.append(graph)

    valid_count = len(graphs) - invalid_count
    print(f"✓ Built {valid_count} valid graphs ({invalid_count} invalid)")

    return graphs
