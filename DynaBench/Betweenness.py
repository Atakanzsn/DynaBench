# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 11:17:19 2024

@author: inant
"""

import networkx as nx
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import itertools
import csv

# Function to read and process each matrix file
def process_matrix_file(file_path):
    G = nx.Graph()
    with open(file_path, 'r') as file:
        for line in file:
            segid1, node1, segid2, node2, weight = line.split()
            node1 = (segid1, int(node1))
            node2 = (segid2, int(node2))
            weight = float(weight)
            G.add_edge(node1, node2, weight=1/weight)
    return G

# Read all matrix_final_i.out files
def calculate_betweenness(input_directory, table_path):
    file_pattern = os.path.join(input_directory, 'matrix_final_*.out')  ## burasi ana dosyadan farkli
    files = glob.glob(file_pattern)

# Aggregate betweenness centrality values across all files
    all_betweenness_centrality = []

    for file_path in files:
        G = process_matrix_file(file_path)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight') ## burayi kontrol et
        all_betweenness_centrality.extend(betweenness_centrality.values())

# Determine the 95th percentile value for betweenness centrality
    quantile_0_95 = np.quantile(all_betweenness_centrality, 0.95)

# Collect nodes with betweenness centrality above the 95th percentile
    top_betweenness_nodes = set()

    for file_path in files:
        G = process_matrix_file(file_path)
        betweenness_centrality = nx.betweenness_centrality(G, weight='weight') ## burayi kontrol et
        for node, centrality in betweenness_centrality.items():
            if centrality >= quantile_0_95:
                top_betweenness_nodes.add((node, centrality))

# Write nodes with 95th percentile betweenness centrality to a file
    top_betweenness_file_path= os.path.join(table_path, 'top_betweenness_nodes.csv')
    with open(top_betweenness_file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Node", "Betweenness Centrality"])
        for node, centrality in sorted(top_betweenness_nodes, key=lambda x: x[1], reverse=True):
            writer.writerow([f"{node[0]}_{node[1]}", centrality])

# Calculate the frequency of each node being in the top 95th percentile
    node_frequencies = {}
    for node, centrality in top_betweenness_nodes:
        if node in node_frequencies:
            node_frequencies[node] += 1
        else:
            node_frequencies[node] = 1

# Write nodes with 95th percentile betweenness centrality and their frequencies to a file
    frequency_file_path= os.path.join(table_path, 'top_betweenness_node_frequencies.csv')
    with open(frequency_file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Node", "Frequency"])
        for node, frequency in node_frequencies.items():
            writer.writerow([f"{node[0]}_{node[1]}",frequency])
# Print confirmation message
    print(f"Nodes and their frequencies written to {frequency_file_path}")

# Write names of the residues (nodes) to a separate CSV file
    residues_file_path = os.path.join(table_path, 'residue_names.csv')
    with open(residues_file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Residue"])
        for node in node_frequencies.keys():
            writer.writerow([f"{node[0]} {node[1]}"])

# Plot the frequencies
#    nodes = list(node_frequencies.keys())
#    frequencies = list(node_frequencies.values())


