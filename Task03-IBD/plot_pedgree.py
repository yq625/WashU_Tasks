# Save this script as plot_pedigree.py

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def read_and_plot_genome(file_path):
    # Read the .genome file into a pandas DataFrame
    genome_df = pd.read_csv(file_path, delim_whitespace=True)

    # Initialize the graph
    G = nx.Graph()

    # Add edges to the graph based on relationships
    for _, row in genome_df.iterrows():
        # Add both nodes and the edge between them
        G.add_edge(row['IID1'], row['IID2'], weight=row['PI_HAT'])

    # Use a layout that spreads nodes out more, such as kamada_kawai_layout
    pos = nx.kamada_kawai_layout(G)

    # Create the plot
    plt.figure(figsize=(15, 15))  # Increase the figure size to spread out nodes

    # Draw the nodes and edges
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color='lightblue')
    nx.draw_networkx_labels(G, pos, font_size=10)

    # Define edge widths based on PI_HAT values, with a cap on the maximum width
    # Edge widths are scaled down for thinner lines
    edge_widths = [d['weight'] * 2 for (u, v, d) in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, width=edge_widths)

    # Annotate the edges with PI_HAT values, adjusting for clarity
    for (u, v, d) in G.edges(data=True):
        edge_x, edge_y = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
        plt.text(edge_x, edge_y, f"{d['weight']:.3f}", fontsize=8, ha='center', va='center', color='darkred')

    # Show the plot
    plt.title('Pedigree Chart Based on IBD Analysis')
    plt.axis('off')

    # Save the plot to a file
    plt.savefig('pedigree_chart_full.png', format='png', bbox_inches='tight', dpi=300)

    plt.show()

if __name__ == '__main__':
    genome_file_path = 'task3_ibd.genome'  # Path to your .genome file
    read_and_plot_genome(genome_file_path)

