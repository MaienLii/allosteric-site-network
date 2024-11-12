import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from enspara.info_theory.mutual_info import weighted_mi

def deviation_mi(devi, weights, threshold_low, threshold_mid1, threshold_mid2, threshold_high):
   
    states = np.zeros_like(devi, dtype=int)

    # Classify states into five categories
    states[devi > threshold_high] = 4  
    states[(devi > threshold_mid2) & (devi <= threshold_high)] = 3 
    states[(devi > threshold_mid1) & (devi <= threshold_mid2)] = 2
    states[(devi > threshold_low) & (devi <= threshold_mid1)] = 1  
    states[devi <= threshold_low] = 0

    devi_mi = weighted_mi(states, weights)

    return devi_mi

def calculate_pocket_mi_network(devi_mi, pocket_residues):
    n_pockets = len(pocket_residues)
    pocket_mi_matrix = np.zeros((n_pockets, n_pockets))

    # Calculate mutual information between pockets
    for i in range(n_pockets):
        for j in range(i, n_pockets):
            residues_i = pocket_residues[i]
            residues_j = pocket_residues[j]
            mi_values = devi_mi[np.ix_(residues_i, residues_j)]
            avg_mi = np.mean(mi_values)  # Calculate average mutual information
            pocket_mi_matrix[i, j] = avg_mi
            pocket_mi_matrix[j, i] = avg_mi  # Ensure matrix symmetry

    return pocket_mi_matrix

def rank_pockets_by_eigenvector_centrality(pocket_mi_matrix):
    # Create a network graph
    G = nx.Graph()
    n_pockets = pocket_mi_matrix.shape[0]

    # Add nodes starting from 1
    G.add_nodes_from(range(1, n_pockets + 1))

    # Add edges based on mutual information matrix
    for i in range(n_pockets):
        for j in range(i + 1, n_pockets):
            if pocket_mi_matrix[i, j] > 0:  # Add edges with non-zero mutual information
                G.add_edge(i + 1, j + 1, weight=pocket_mi_matrix[i, j])  # Nodes are indexed starting from 1

    # Calculate eigenvector centrality
    centrality = nx.eigenvector_centrality_numpy(G, weight='weight')

    # Sort pockets by centrality
    ranked_pockets = sorted(centrality.items(), key=lambda x: x[1], reverse=True)

    return ranked_pockets, G

# Updated function to select and plot six nodes
def plot_top_pockets_network(cluster_graph, centrality, k_value=0.1, custom_colors=None, font_size=10):
    pos = nx.spring_layout(cluster_graph, k=k_value)  # Use spring layout

    # Adjust node size based on centrality
    node_size = [v * 30000 for v in centrality.values()]

    # Check if color list is None, otherwise use default colors
    if custom_colors is None:
        custom_colors = ['#40E0D0', '#EE82EE', '#FF0000', '#FFFF00', '#FFA500', '#0000FF', '#0000FF']  # Color adjustments

    # Select the top 7 nodes
    top_seven_nodes = sorted(centrality, key=centrality.get, reverse=True)[:7]

    # Assign colors to nodes
    node_color = [custom_colors[top_seven_nodes.index(node)] if node in top_seven_nodes else 'gray'
                  for node in cluster_graph.nodes()]

    # Create a subgraph containing the top 7 nodes
    subgraph = cluster_graph.subgraph(top_seven_nodes).copy()

    # Extract and normalize edge weights for line widths
    all_weights = np.array([subgraph[u][v]['weight'] for u, v in subgraph.edges()])
    min_weight, max_weight = all_weights.min(), all_weights.max()
    edge_widths = 2 + 8 * (all_weights - min_weight) / (max_weight - min_weight)

    # Plot the graph
    fig, ax = plt.subplots(figsize=(12, 12))

    # Remove external borders
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.margins(0.1)

    # Draw edges
    nx.draw_networkx_edges(subgraph, pos, width=edge_widths, edge_color='gray')

    # Draw nodes
    nx.draw_networkx_nodes(
        subgraph, pos,
        node_size=[node_size[top_seven_nodes.index(node)] for node in top_seven_nodes],
        node_color=[custom_colors[top_seven_nodes.index(node)] for node in top_seven_nodes]
    )

    # Draw labels: ensure each node label matches its order in top_seven_nodes, label the sixth node as "AVP"
    labels = {node: (f"{i+1}" if i != 5 else "AVP") for i, node in enumerate(top_seven_nodes)}
    nx.draw_networkx_labels(subgraph, pos, labels=labels, font_size=font_size, font_color="black")

    # Save image
    plt.savefig('networkshi.tif', format='tiff')
    plt.show()

def plot_centrality_ranking(ranked_pockets, output_filename="pocket_centrality_ranking.tif"):
    # Extract pocket indices and centrality scores
    pockets, centrality_scores = zip(*ranked_pockets)

    # Create bar plot
    hex_color = '#FCE6C9'
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.bar(range(len(pockets)), centrality_scores, color=hex_color)
    
    # Add labels and title, adjust labelpad to space out labels from plot
    ax.set_xlabel('Pockets', fontsize=28, labelpad=30)  # Adjust x-axis label distance with labelpad
    ax.set_ylabel('Eigenvector Centrality', fontsize=28, labelpad=30)  # Adjust y-axis label distance with labelpad
    ax.set_title('Pocket Ranking by Eigenvector Centrality', fontsize=28, pad=30)

    ax.set_xticks(range(len(pockets)))
    ax.set_xticklabels(pockets, rotation=90, fontsize=25)
    ax.tick_params(axis='y', labelsize=25)

    # Save image
    plt.tight_layout()
    plt.savefig(output_filename, format='tiff', dpi=600)
    plt.show()


# Main program
devi = np.loadtxt('./repre2.txt')
damping = 0.85
weights = np.loadtxt('./weight2.txt')
threshold_low = 0.356
threshold_mid1 = 1.53
threshold_mid2 = 1.96
threshold_high = 10.84
pocket_residues = []

# Read pocket residues file
with open('pocket.txt', 'r') as file:
    for line in file:
        residues = line.strip().split()
        residues = [int(residue) for residue in residues]
        pocket_residues.append(residues)

devi_mi = deviation_mi(devi, weights, threshold_low, threshold_mid1, threshold_mid2, threshold_high)

pocket_mi_matrix = calculate_pocket_mi_network(devi_mi, pocket_residues)

ranked_pockets, G = rank_pockets_by_eigenvector_centrality(pocket_mi_matrix)

# Output ranked pockets, starting from 1
print("Ranked pockets by eigenvector centrality:")
for rank, (pocket, score) in enumerate(ranked_pockets):
    print(f"Rank {rank + 1}: Pocket {pocket}, Centrality Score: {score:.4f}")

# Plot the network graph of the top 6 pockets with custom parameters
custom_colors = ['#40E0D0', '#FFFF00', '#FFA500', '#EE82EE', '#FA7F6F', '#82B0D2', '#52BE80']
plot_top_pockets_network(G, dict(ranked_pockets), k_value=0.1, custom_colors=custom_colors, font_size=35)

plot_centrality_ranking(ranked_pockets, output_filename="pocket_centrality_ranking.tif")
