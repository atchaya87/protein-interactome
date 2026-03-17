## Overview
The goal of this project is to identify regulatory switches within the interactions between the JAK/STAT pathway and Huntington (HTT) pathway. The project takes in STRING database human proteomics data, differentiating between locally functioning proteins and proteins acting as bridges between the signaling pathways of interest. 

## Background

### Proteins as Graphs
Proteins do not act in isolation--they interact with one another in signaling pathways and between various pathways. They form an interactome that can be modeled via a graph, where proteins are nodes and their interactions are edges. Therefore, graph theory can be used to connect various molecular mechanisms and find proteins-of-interest to research if a pathway is observed to fail.

### Bridge Proteins
Proteins exist in specialized families or communities, performing similar, highly-specific functions. However, a single protein may be crucial to initiating and participating in multiple pathways, acting as so-called bridges. Within graphs, bridge proteins can be identified using *Betweenness Centrality*, measuring how frequently a node appears in the shortest path between any two nodes. 

### JAK/STAT and HTT
- The JAK/STAT signaling pathway is key in regulating gene transcription for pro-inflammatory signaling, cell division, and cell death. Pathway initiation depends on a cytokine ligand binding to a a pair of cytokine receptors, causing receptor dimerization. The dimer is subsequently phosphorylated by JAK proteins, providing binding scaffolding for STAT proteins which are auto-phosphorylated by JAKs. The STATs can then dimerize and enter the nucleus, binding DNA and leading to specific gene transcription.
- Huntingtin, or HTT, is a multi-domain protein necessary for the development of the nervous system and cellular transport. Furthermore, it is vital for synaptic health by promoting the production of BNDF, a neurotrophic factor influential in sustained neuron growth and survival. Huntington's disease is strongly linked to the misfolding of HTT as caused by polyglutamine expansion in its associated gene.

In molecular biology, inflammatory signaling (JAK/STAT) and neurodegeneration (HTT malfunction) are often connected. By linking the bridge proteins between two important pathways, we can identify nodes where a signaling error triggers a collapse across the body, leading to disease such as Huntington's. 

## Implementation

### Data Sourcing and Pre-Processing
The project utilizes the STRING Human Proteomics Database (v12.0). Raw interaction datasets are parsed using a modular configuration system (config.yaml), allowing for targeted analysis of specific protein clusters. A subgraph containing proteins of interest associated with JAK/STAT and HTT is generated and used, only considering interactions above a strength threshold. 

### Manipulation of Confidence Scores
STRING Database confidence scores, S, represent the probability of a protein-protein interaction and range from 0-1000. In a graph, higher-probability interactions should cause nodes to be closer together. So, the pipeline converts scores into a probability-based distance, D: D = 1/S. 

Then, functional clusters are easily isolated by the force-directed physics engine as high-confidence interactor proteins are pulled together and visually separated.

### Betweenness Centrality 
The pipeline then calculates Betweenness Centrality to identify bridge proteins. For any node w, the Betweenness Centrality is given by the fraction of shortest paths between any two nodes that pass through w. The project focuses on the nodes with the highest betweenness, marking the associated proteins as supposed regulatory switches facilitating JAK/STAT and HTT interaction. 

## Results and Visualization 
### Interactome Dashboard
The final output is a D3.js-based dashboard allowing for real-time exploration of the calculated interactome. The pipeline also supports a simple 2D projection of the protein interactions through matplotlib, ideal for identifying general neighborhoods of proteins. The 3D visualization is a force-directed layout where the physical distance between nodes is reflects the probability-based distance, D.

The radius of each node is proportional to its Betweenness Centrality score, with larger nodes indicating a higher score and increased role as a bridge protein. The physics engine naturally segregates the network into functional clusters. In this specific analysis, the pro-inflammatory JAK/STAT signaling module and the HTT transport module appear as distinct communities, joined by bridge nodes. The top 10 proteins with the highest Betweenness Centrality are color-coded and scaled.

Every node is interactive. Clicking a protein node opens its associated UniProt database entry, providing access to its subcellular location, molecular function, and sequence data across species.

## Top bridge proteins identified:
1. STAT3
2. HTT
3. MAPK1
4. MAPK3
5. TP53
6. AKT1
7. EGFR
8. IL6
9. JAK2
10. STAT1

## Usage

1. Install the necessary dependencies via pip:
```
 pip install pandas networkx pyvis pyyaml
```

2. Adjust targets and confidence thresholds in config.yaml; refer to the following example:
```
targets: ["JAK2", "STAT3", "HTT"]
confidence_threshold: 700
top_n_highlights: 10
```

3. Execution
Run the pipeline to generate the interactome.

```
python main.py
```
4. Open the generated .html file in any browser:

- Bridge Nodes: The top 10 bridge proteins are scaled and color-coded by betweenness centrality.

- Modularity: Drag nodes to explore the physics-based cluster separation.

- UniProt: Click any node to open its corresponding biological profile.
  
- Customization: Adjust the font, color palette, and node opacity directly in the dashboard.

## Next Steps
While the current interactome models protein-protein interaction in terms of fixed, given confidence scores, it must be validated against live data. Furthermore, the pipeline could implement a targeted knock out feature on bridge proteins, measuring the protein's impact on the global interactome and highlight which pathways fail as a result. The intensity of the knockout could be assigned as scores to the proteins, demonstrating which bridge proteins are most important between which pathways. The modular system could further be configured to other neurodegenerative diseases such as Alzheimer's by creating a different subgraph.




