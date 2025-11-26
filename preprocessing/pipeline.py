"""
Main preprocessing pipeline
Orchestrates all preprocessing steps in correct order
"""

import json
import pickle
import numpy as np
import torch
from pathlib import Path

from .config import OUTPUT_DIR, QUANTUM_FEATURES, MOLECULAR_FEATURES, DERIVED_FEATURES, TARGET_COLUMNS
from .config import QC_METHOD, INCLUDE_EXPLICIT_H, SPLIT_METHOD, SPLIT_RATIOS, SCALING_METHOD, HANDLE_MISSING, RANDOM_SEED
from .data_loader import load_raw_data
from .data_integration import integrate_data
from .graph_builder import build_all_graphs
from .feature_engineering import engineer_features
from .target_preparation import prepare_targets
from .data_validation import validate_and_clean
from .data_splitting import scaffold_split
from .feature_scaling import scale_features
from .graphormer_encoding import batch_encode_for_graphormer


def save_preprocessed_data(df, graphs, graphormer_graphs, scalers, splits, output_dir=OUTPUT_DIR):
    """
    Save all preprocessed data with comprehensive metadata

    Args:
        df: Processed DataFrame
        graphs: List of PyG graph objects
        graphormer_graphs: List of Graphormer-encoded graphs
        scalers: Dict with 'feature' and 'target' scalers
        splits: Dict with 'train', 'val', 'test' indices
        output_dir: Output directory path
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    print(f"\nSaving preprocessed data to {output_dir}...")

    # 1. Save master CSV
    df.to_csv(output_path / 'molecules.csv', index=False)
    print(f"✓ Saved molecules.csv ({len(df)} molecules)")

    # 2. Save PyG graphs
    graphs_dir = output_path / 'graphs'
    graphs_dir.mkdir(exist_ok=True)
    for idx, graph in enumerate(graphs):
        if graph is not None:
            torch.save(graph, graphs_dir / f'mol_{idx}.pt')
    print(f"✓ Saved {len(graphs)} PyG graph objects")

    # 3. Save Graphormer-encoded graphs
    graphormer_dir = output_path / 'graphormer_inputs'
    graphormer_dir.mkdir(exist_ok=True)
    for idx, encoded_graph in enumerate(graphormer_graphs):
        if encoded_graph is not None:
            torch.save(encoded_graph, graphormer_dir / f'mol_{idx}.pt')
    print(f"✓ Saved {len(graphormer_graphs)} Graphormer-encoded graphs")

    # 4. Save scalers
    scalers_dir = output_path / 'scalers'
    scalers_dir.mkdir(exist_ok=True)
    with open(scalers_dir / 'feature_scaler.pkl', 'wb') as f:
        pickle.dump(scalers['feature'], f)
    with open(scalers_dir / 'target_scaler.pkl', 'wb') as f:
        pickle.dump(scalers['target'], f)
    print(f"✓ Saved scalers")

    # 5. Save splits
    splits_dir = output_path / 'splits'
    splits_dir.mkdir(exist_ok=True)
    np.save(splits_dir / 'train_indices.npy', splits['train'])
    np.save(splits_dir / 'val_indices.npy', splits['val'])
    np.save(splits_dir / 'test_indices.npy', splits['test'])
    print(f"✓ Saved splits")

    # 6. Save metadata
    metadata = {
        'preprocessing_config': {
            'qc_method': QC_METHOD,
            'include_explicit_h': INCLUDE_EXPLICIT_H,
            'split_method': SPLIT_METHOD,
            'split_ratios': SPLIT_RATIOS,
            'scaling_method': SCALING_METHOD,
            'handle_missing': HANDLE_MISSING,
            'random_seed': RANDOM_SEED,
        },
        'dataset_statistics': {
            'n_molecules': len(df),
            'n_with_targets': int(df['has_target'].sum()),
            'feature_columns': QUANTUM_FEATURES + MOLECULAR_FEATURES + DERIVED_FEATURES,
            'target_columns': TARGET_COLUMNS,
        },
        'split_statistics': {
            'train_size': len(splits['train']),
            'val_size': len(splits['val']),
            'test_size': len(splits['test']),
        },
        'graphormer_info': {
            'max_nodes': int(max(g['num_nodes'] for g in graphormer_graphs if g is not None)),
            'avg_nodes': float(np.mean([g['num_nodes'] for g in graphormer_graphs if g is not None])),
        }
    }

    with open(output_path / 'metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"✓ Saved metadata.json")

    print(f"\n{'='*60}")
    print("✅ ALL PREPROCESSED DATA SAVED")
    print(f"{'='*60}\n")


def run_preprocessing_pipeline():
    """
    Execute complete preprocessing pipeline
    """
    print("\n" + "="*80)
    print("CEPDB PREPROCESSING PIPELINE FOR GRAPHORMER")
    print("="*80 + "\n")

    # Step 1-2: Load and integrate data
    print("[STEP 1-2] Loading and integrating data...")
    raw_data = load_raw_data()
    df = integrate_data(raw_data)

    # Step 3: Build molecular graphs
    print(f"\n[STEP 3] Building molecular graphs...")
    graphs = build_all_graphs(df)

    # Step 4: Feature engineering
    print(f"\n[STEP 4] Engineering features...")
    df = engineer_features(df)

    # Step 5: Prepare targets
    print(f"\n[STEP 5] Preparing targets...")
    df, graphs = prepare_targets(df, graphs)

    # Step 6: Validate and clean
    print(f"\n[STEP 6] Validating and cleaning...")
    df, graphs = validate_and_clean(df, graphs)

    # Step 7: Split data (BEFORE scaling!)
    print(f"\n[STEP 7] Splitting data...")
    train_idx, val_idx, test_idx = scaffold_split(df)
    splits = {'train': train_idx, 'val': val_idx, 'test': test_idx}

    # Step 8: Scale features
    print(f"\n[STEP 8] Scaling features...")
    df, scalers = scale_features(df, train_idx)

    # GRAPHORMER-SPECIFIC: Encode graphs
    print(f"\n[STEP 9 - GRAPHORMER] Encoding for Graphormer...")
    graphormer_graphs = batch_encode_for_graphormer(graphs)

    # Step 10: Save everything
    print(f"\n[STEP 10] Saving preprocessed data...")
    save_preprocessed_data(df, graphs, graphormer_graphs, scalers, splits)

    print("\n" + "="*80)
    print("✅ PREPROCESSING COMPLETE!")
    print("="*80)
    print(f"\nOutput location: {OUTPUT_DIR}")
    print(f"Total molecules: {len(df)}")
    print(f"  Train: {len(train_idx)}")
    print(f"  Val:   {len(val_idx)}")
    print(f"  Test:  {len(test_idx)}")

    return df, graphs, graphormer_graphs, scalers, splits


if __name__ == '__main__':
    run_preprocessing_pipeline()
