import torch
from transformers import AutoModel, AutoConfig

# ---------------------------------------------------------
# 1. Choose a GNN model from Hugging Face
#    Examples:
#      - "facebook/graphormer-base"
#      - "deepmind/graphnet"
#      - Any model with graph support
# ---------------------------------------------------------

model_name = "facebook/graphormer-base"

# Load config + model
config = AutoConfig.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)

# Put model on GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)


# ---------------------------------------------------------
# 2. Prepare input graph
#    This will vary depending on the model.
#    Most Hugging Face GNN models use:
#      - node_features: [num_nodes, feature_dim]
#      - edge_index: [2, num_edges]
#      - edge_features (optional)
# ---------------------------------------------------------

# Example dummy graph
node_features = torch.randn(10, config.hidden_size)  # 10 nodes
edge_index = torch.tensor(
    [
        [0, 1, 2, 3],  # source nodes
        [1, 2, 3, 4],  # target nodes
    ]
)

# Some models want batch dims
node_features = node_features.unsqueeze(0)
edge_index = edge_index.unsqueeze(0)

node_features = node_features.to(device)
edge_index = edge_index.to(device)

# ---------------------------------------------------------
# 3. Forward pass
# ---------------------------------------------------------
with torch.no_grad():
    outputs = model(
        node_features=node_features,
        edge_index=edge_index,
        # Some models want extra kwargs like:
        #   edge_attr=edge_features
        #   attention_mask=...
    )

print(outputs)
print("Output shape:", outputs.last_hidden_state.shape)
