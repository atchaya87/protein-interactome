import pandas as pd
import networkx as nx
import community as community_louvain
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml  
from pyvis.network import Network
import os

# 0. load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# 1. load the file
# 2. create a mapping from id to name
# 3. create a mapping from name to id for traversal 
info_df = pd.read_csv(config['info_path'], sep="\t")
id_to_name = dict(zip(info_df["#string_protein_id"], info_df["preferred_name"]))
name_to_id = {v: k for k, v in id_to_name.items()}

# interaction data
df = pd.read_csv(config['links_path'], sep=" ")

# filter for strong (using config value)
df = df[df["combined_score"] > config['min_confidence']]
print(f"Strong interactions > {config['min_confidence']}:", len(df))

# create graph
G = nx.from_pandas_edgelist(df, "protein1", "protein2", ["combined_score"])                        
print("Number of proteins:", G.number_of_nodes())
print("Number of interactions:", G.number_of_edges())

# create a subgraph
important_proteins = config['important_proteins']
important_ids = [name_to_id[name] for name in important_proteins if name in name_to_id]

neighborhood = set()
for prot in important_ids: 
    if prot in G:
        neighborhood.update(G.neighbors(prot))
        neighborhood.add(prot)

# convert to a full graph copy to avoid nx.draw issues
H = G.subgraph(neighborhood).copy()

# invert combined_score for weighted betweenness (stronger interactions = shorter paths)
for u, v, d in H.edges(data=True):
    d['weight'] = 1 / d['combined_score']

# Bridge Score/betweenness - how many paths go through the node
scores = nx.betweenness_centrality(H, weight='weight')

# Sort them from highest to lowest, sort by score and not id for the top 10
top_switches = sorted(scores.items(), key=lambda x: x[1], reverse=True)[:10]
top_switch_ids = [p[0] for p in top_switches]

# find the most connected proteins
centrality = nx.degree_centrality(G)
top_proteins = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:5]

# subgraph results
print("\nTop 10 Regulatory 'Switches' in Subgraph H:")
for p_id, score in top_switches:
    print(f"{id_to_name.get(p_id, p_id)}: {score:.4f}")

# community detection
partition = community_louvain.best_partition(H, weight='combined_score', random_state=config['random_seed'])

# pick a color for each community safely
node_colors = [partition[n] for n in H.nodes()]

# subgraph with communities highlighted: visualize communities via different colors and node sizes
plt.figure(figsize=(12, 10), facecolor='white')
pos = nx.spring_layout(H, k=config['layout_k'], iterations=config['layout_iterations'], seed=config['random_seed'])
node_sizes = [scores[n] * 10000 + 100 for n in H.nodes()]

n_labels_to_show = {n: id_to_name.get(n, n) for n in H.nodes() if n in top_switch_ids or n in important_ids}

#nodes
nx.draw_networkx_nodes(H, pos, node_size=node_sizes, node_color=node_colors, alpha=0.9, cmap='viridis')

#edges
nx.draw_networkx_edges(H, pos, alpha=0.1, edge_color='gray')

#labels for the top switches and proteins of interest
nx.draw_networkx_labels(H, pos, labels=n_labels_to_show, font_size=11, font_weight='bold')

plt.title("Important Switch Proteins in the HTT-STAT3 Network", fontsize=16)
plt.axis('off')
plt.show()

# pyvis
print("Generating interactive map...")
net = Network(height='750px', width='100%', bgcolor='#222222', font_color='white')

cmap = plt.get_cmap('viridis', max(partition.values()) + 1)

for n in H.nodes():
    comm_id = partition[n]
    color_hex = cm.colors.rgb2hex(cmap(comm_id))
    full_name = id_to_name.get(n, n)
    switch_score = scores.get(n, 0)
    
    # link to the protein's research page for further interaction
    # uniprot.org/uniprotkb/PROTEIN_NAME/entry
    uniprot_url = f"https://www.uniprot.org/uniprotkb?query={full_name}"
    
    # HTML
    hover_html = f"""
    <div style="font-family: sans-serif; padding: 10px;">
        <h3 style="margin: 0; color: {color_hex};">{full_name}</h3>
        <hr>
        <b>Community:</b> {comm_id}<br>
        <b>Bridge Score:</b> {switch_score:.4f}<br>
        <br>
        <a href="{uniprot_url}" target="_blank" style="color: #4db8ff;">View Protein Data ↗</a>
    </div>
    """
    
    display_label = full_name if (n in important_ids or n in top_switch_ids) else ""
    
    net.add_node(n, label=display_label, title=hover_html, color=color_hex, size=(switch_score * 500) + 15)
net.show_buttons(filter_=['physics', 'nodes'])
for u, v, d in H.edges(data=True):
    net.add_edge(u, v, value=d['combined_score']/1000, color="gray", alpha=0.3)

net.force_atlas_2based()
net.show("protein_interactome.html", notebook=False)

import os
print(f"\nOpen this file: {os.path.abspath('protein_interactome.html')}")
