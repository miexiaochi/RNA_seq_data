import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import community as community_louvain
from collections import defaultdict
import numpy as np

# Read data
file_path = "node.weight.txt"  
data = pd.read_csv(file_path, sep='\t')

# Create graph
G = nx.from_pandas_edgelist(data, 'fromNode', 'toNode', ['weight'])

# Calculate node degrees
degree_dict = dict(G.degree())
degrees = pd.Series(degree_dict)

# Set degree threshold
#max=213
degree_threshold = 200  

# Remove nodes below the threshold
nodes_to_keep = degrees[degrees > degree_threshold].index
G_filtered = G.subgraph(nodes_to_keep).copy()

# Perform community detection using the Louvain algorithm
partition = community_louvain.best_partition(G_filtered)

# Extract community information
community_nodes = defaultdict(list)
for node, comm_id in partition.items():
    community_nodes[comm_id].append(node)

# Plotting
plt.figure(figsize=(15, 15))

# Set circular layout
circle_radius = 1  # Radius of each community circle
offset = 3  # Spacing between circles
community_positions = {}

# Assign positions to each community and plot nodes
for comm_id, nodes in community_nodes.items():
    # Calculate the center position of the circle
    angle = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)  # Generate angles
    center_x = offset * comm_id * np.cos(np.pi / len(community_nodes))  # Set X coordinate of the circle center
    center_y = offset * comm_id * np.sin(np.pi / len(community_nodes))  # Set Y coordinate of the circle center
    
    # Calculate the position of each node
    for i, node in enumerate(nodes):
        x = center_x + circle_radius * np.cos(angle[i])
        y = center_y + circle_radius * np.sin(angle[i])
        community_positions[node] = (x, y)  # Save node position

# Plot nodes
nx.draw_networkx_nodes(G_filtered, pos=community_positions, node_size=800, 
                         node_color=[('red' if partition[node] % 2 == 0 else 'green') for node in G_filtered.nodes()])

# Plot community boundary circles
for comm_id, nodes in community_nodes.items():
    center_x = offset * comm_id * np.cos(np.pi / len(community_nodes))
    center_y = offset * comm_id * np.sin(np.pi / len(community_nodes))
    circle = plt.Circle((center_x, center_y), circle_radius, 
                        color='red' if comm_id % 2 == 0 else 'green', fill=False, linewidth=2)
    plt.gca().add_patch(circle)

# Plot edges
for edge in G_filtered.edges():
    node1, node2 = edge
    # Determine edge color
    if partition[node1] == partition[node2]:  # Edges within the same community
        if partition[node1] % 2 == 0:  # Red community
            color = 'red'
        else:  # Green community
            color = 'green'
    else:  # Edges between different communities
        color = 'gray'  # Default color
        if partition[node1] % 2 == 0:  # If node1 belongs to the red community
            color = 'blue'
        elif partition[node1] % 2 == 1:  # If node1 belongs to the green community
            color = 'orange'
    
    # Only plot edges that connect nodes within communities
    if node1 in community_positions and node2 in community_positions:
        nx.draw_networkx_edges(G_filtered, pos=community_positions,
                                 edgelist=[edge], alpha=0.5, edge_color=color)

# Add gene name labels, larger font and brighter color
for node, pos in community_positions.items():
    plt.text(pos[0], pos[1], node, fontsize=20, ha='center', va='center', color='black')  # Larger font and brighter color

plt.axis('off')  # Turn off axis
plt.show()