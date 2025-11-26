"""
Graphormer-specific graph encoding

Converts PyTorch Geometric graphs to Graphormer input format:
- Spatial encoding: Shortest path distances between all node pairs
- Centrality encoding: In-degree and out-degree for each node
- Edge encoding: Features along shortest paths
- Attention bias: Precomputed bias matrices
"""

import torch
import numpy as np
from scipy.sparse.csgraph import shortest_path
from torch_geometric.utils import to_dense_adj, degree


def compute_spatial_pos(edge_index, num_nodes):
    """
    Compute shortest path distances between all node pairs

    Graphormer uses this for spatial encoding in attention mechanism

    Args:
        edge_index: Edge connectivity [2, num_edges]
        num_nodes: Number of nodes in graph

    Returns:
        torch.LongTensor: [num_nodes, num_nodes] distance matrix
    """
    # Convert edge_index to dense adjacency matrix
    adj = to_dense_adj(edge_index, max_num_nodes=num_nodes)[0]

    # Compute shortest paths using Floyd-Warshall
    adj_np = adj.cpu().numpy()

    # scipy expects 0 for no edge, we need to set non-edges to inf
    adj_np[adj_np == 0] = np.inf
    np.fill_diagonal(adj_np, 0)  # Distance to self is 0

    distances = shortest_path(adj_np, directed=False, method='FW')

    # Replace inf with max_path_len (common practice in Graphormer)
    max_path_len = num_nodes  # or a fixed value like 512
    distances[np.isinf(distances)] = max_path_len

    return torch.from_numpy(distances).long()


def compute_degree_features(edge_index, num_nodes):
    """
    Compute in-degree and out-degree for centrality encoding

    For undirected graphs, in-degree == out-degree

    Args:
        edge_index: Edge connectivity [2, num_edges]
        num_nodes: Number of nodes in graph

    Returns:
        tuple: (in_degree, out_degree) as LongTensors [num_nodes]
    """
    # Compute degree
    deg = degree(edge_index[0], num_nodes=num_nodes, dtype=torch.long)

    # For undirected graphs
    in_degree = deg
    out_degree = deg

    return in_degree, out_degree


def compute_attn_bias(edge_index, edge_attr, num_nodes, spatial_pos):
    """
    Compute attention bias matrix

    Encodes edge features and spatial distances into attention scores
    This is a simplified version - full Graphormer uses more complex encoding

    Args:
        edge_index: Edge connectivity [2, num_edges]
        edge_attr: Edge features [num_edges, edge_feature_dim]
        num_nodes: Number of nodes
        spatial_pos: Shortest path distances [num_nodes, num_nodes]

    Returns:
        torch.FloatTensor: [num_nodes, num_nodes] attention bias
    """
    # Initialize bias matrix
    attn_bias = torch.zeros(num_nodes, num_nodes)

    # Add bias for direct edges (from edge_attr)
    for idx in range(edge_index.size(1)):
        i, j = edge_index[0, idx], edge_index[1, idx]
        if edge_attr is not None and edge_attr.size(0) > 0:
            # Use bond type as bias (first feature)
            attn_bias[i, j] = edge_attr[idx, 0]

    # Add spatial bias (negative distance encourages attention to nearby nodes)
    spatial_bias = -spatial_pos.float() * 0.1  # Scale factor
    attn_bias = attn_bias + spatial_bias

    return attn_bias


def compute_edge_input(edge_index, edge_attr, num_nodes, spatial_pos, max_dist=5):
    """
    Compute edge input encoding

    For Graphormer: encode edge features along shortest paths

    Args:
        edge_index: Edge connectivity [2, num_edges]
        edge_attr: Edge features [num_edges, edge_feature_dim]
        num_nodes: Number of nodes
        spatial_pos: Shortest path distances [num_nodes, num_nodes]
        max_dist: Maximum distance to consider

    Returns:
        torch.FloatTensor: [num_nodes, num_nodes, edge_feature_dim]
    """
    if edge_attr is None or edge_attr.size(0) == 0:
        # No edge features - use distance only
        edge_feature_dim = 1
        edge_input = spatial_pos.unsqueeze(-1).float()
    else:
        edge_feature_dim = edge_attr.size(1)
        edge_input = torch.zeros(num_nodes, num_nodes, edge_feature_dim)

        # Fill in direct edge features
        for idx in range(edge_index.size(1)):
            i, j = edge_index[0, idx], edge_index[1, idx]
            edge_input[i, j] = edge_attr[idx]

    return edge_input


def encode_for_graphormer(graph_data):
    """
    Convert PyG Data object to Graphormer input format

    Args:
        graph_data: torch_geometric.data.Data object

    Returns:
        dict: {
            'input_nodes': Node features [num_nodes, node_feature_dim],
            'input_edges': Edge input [num_nodes, num_nodes, edge_feature_dim],
            'attn_bias': Attention bias [num_nodes, num_nodes],
            'spatial_pos': Spatial encoding [num_nodes, num_nodes],
            'in_degree': In-degree [num_nodes],
            'out_degree': Out-degree [num_nodes],
            'num_nodes': int
        }
    """
    edge_index = graph_data.edge_index
    edge_attr = graph_data.edge_attr
    x = graph_data.x
    num_nodes = x.size(0)

    # Check Graphormer limitation
    if num_nodes > 100:
        print(f"WARNING: Graph has {num_nodes} nodes (>100). "
              f"Graphormer may struggle with large graphs.")

    # Compute Graphormer-specific encodings
    spatial_pos = compute_spatial_pos(edge_index, num_nodes)
    in_degree, out_degree = compute_degree_features(edge_index, num_nodes)
    attn_bias = compute_attn_bias(edge_index, edge_attr, num_nodes, spatial_pos)
    edge_input = compute_edge_input(edge_index, edge_attr, num_nodes, spatial_pos)

    return {
        'input_nodes': x,
        'input_edges': edge_input,
        'attn_bias': attn_bias,
        'spatial_pos': spatial_pos,
        'in_degree': in_degree,
        'out_degree': out_degree,
        'num_nodes': num_nodes
    }


def batch_encode_for_graphormer(graphs):
    """
    Encode a list of graphs for Graphormer

    Args:
        graphs: List of PyG Data objects

    Returns:
        list: List of encoded graph dictionaries
    """
    print("Encoding graphs for Graphormer...")

    encoded_graphs = []
    large_graph_count = 0

    for idx, graph in enumerate(graphs):
        if graph is None:
            encoded_graphs.append(None)
            continue

        encoded = encode_for_graphormer(graph)

        if encoded['num_nodes'] > 100:
            large_graph_count += 1

        encoded_graphs.append(encoded)

    if large_graph_count > 0:
        print(f"  WARNING: {large_graph_count} graphs have >100 nodes")

    print(f"✓ Encoded {len(encoded_graphs)} graphs for Graphormer")

    return encoded_graphs
