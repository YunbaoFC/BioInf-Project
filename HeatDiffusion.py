#!/usr/bin/env python

"""
Heat Diffusion Algorithm
A diffusion-based method for disease gene prioritization using heat kernel
"""

import sys
import csv
import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


def print_usage():
    print(' ')
    print('        usage: python3 HeatDiffusion.py network_file seed_file n t(optional) outfile_name(optional)')
    print('        -----------------------------------------------------------------')
    print('        network_file : The edgelist must be provided as any delimiter-separated')
    print('                       table. Make sure the delimiter does not exist in gene IDs')
    print('                       The first two columns will be interpreted as gene1 <==> gene2')
    print('        seed_file    : table containing the seed genes (first column will be used)')
    print('        n            : desired number of top-ranked genes to output')
    print('        t            : diffusion time parameter (float), default is 0.5')
    print('                       small t (0.1-0.3) = stays close to seeds')
    print('                       large t (0.8-2.0) = diffuses far in network')
    print('        outfile_name : results will be saved under this file name')
    print('                       default: "HeatDiffusion_top_n_genes_time_X.txt"')
    print(' ')


def check_input_style(input_list):
    try:
        network_file = input_list[1]
        seeds_file = input_list[2]
        top_n = int(input_list[3])
    except:
        print_usage()
        sys.exit(0)
        
    t = 0.5
    outfile = f'HeatDiffusion_top_{top_n}_genes_time_{t}.txt'
    
    if len(input_list) == 5:
        try:
            t = float(input_list[4])
            outfile = f'HeatDiffusion_top_{top_n}_genes_time_{t}.txt'
        except:
            outfile = input_list[4]
    
    if len(input_list) == 6:
        try:
            t = float(input_list[4])
            outfile = input_list[5]
        except:
            print_usage()
            sys.exit(0)
            
    return network_file, seeds_file, top_n, t, outfile


def read_input(network_file, seed_file):
    """
    Reads network and seed genes from files
    """
    # Detect delimiter
    sniffer = csv.Sniffer()
    line_delimiter = None
    for line in open(network_file, 'r'):
        if line[0] == '#':
            continue
        dialect = sniffer.sniff(line)
        line_delimiter = dialect.delimiter
        break
    
    if line_delimiter is None:
        print('network_file format not correct')
        sys.exit(0)
    
    # Read network
    G = nx.Graph()
    for line in open(network_file, 'r'):
        if line[0] == '#':
            continue
        line_data = line.strip().split(line_delimiter)
        node1 = line_data[0]
        node2 = line_data[1]
        weight = float(line_data[2]) if len(line_data) > 2 else 1.0
        G.add_edge(node1, node2, weight=weight)
    
    # Read seed genes
    seed_genes = set()
    for line in open(seed_file, 'r'):
        if line[0] == '#':
            continue
        line_data = line.strip().split('\t')
        seed_gene = line_data[0]
        seed_genes.add(seed_gene)
    
    return G, seed_genes


def HeatDiffusion(G_original, seed_genes, top_n, t, outfile=None):
    """
    Heat Diffusion algorithm using heat kernel
    
    Parameters:
    -----------
    G_original : networkx.Graph
        The network
    seed_genes : set
        Set of seed genes
    top_n : int
        Number of top genes to output
    t : float
        Diffusion time parameter
    outfile : str
        Output filename
        
    Returns:
    --------
    results : list of tuples
        (rank, gene, score)
    """
    
    # Filter seeds that are in the network
    all_genes = set(G_original.nodes())
    seed_genes = set(seed_genes)
    valid_seeds = seed_genes & all_genes
    
    if len(valid_seeds) != len(seed_genes):
        print(f"HeatDiffusion(): ignoring {len(seed_genes - all_genes)} of {len(seed_genes)} seed genes not in network")
    
    if len(valid_seeds) == 0:
        print("ERROR: No valid seed genes in the network!")
        sys.exit(1)
    
    print(f"Network: {G_original.number_of_nodes()} nodes, {G_original.number_of_edges()} edges")
    print(f"Valid seeds: {len(valid_seeds)}")
    print(f"Running Heat Diffusion with t={t}...")
    
    # Prepare node list
    all_genes_list = list(all_genes)
    gene_to_idx = {gene: idx for idx, gene in enumerate(all_genes_list)}
    
    # Get the normalized Laplacian matrix
    L = nx.normalized_laplacian_matrix(G_original, nodelist=all_genes_list, weight='weight')
    
    # Initialize heat vector (seeds)
    h0 = np.zeros(len(all_genes_list))
    for seed in valid_seeds:
        h0[gene_to_idx[seed]] = 1.0
    h0 /= len(valid_seeds)  # Normalize
    
    # Compute heat diffusion: h(t) = exp(-t * L) * h0
    # Using expm_multiply for efficient computation
    print("Computing heat kernel...")
    ht = expm_multiply(-t * L, h0)
    
    # Normalize scores to sum to 1
    ht = np.abs(ht)
    if np.sum(ht) > 0:
        ht /= np.sum(ht)
    
    print("Heat diffusion completed")
    
    # Get results: exclude seeds from output
    results = []
    for idx, gene in enumerate(all_genes_list):
        if gene not in valid_seeds:  # Exclude seeds
            results.append((gene, ht[idx]))
    
    # Sort by score (descending)
    results.sort(key=lambda x: x[1], reverse=True)
    
    # Take top N
    top_results = results[:top_n]
    
    # Save results
    with open(outfile, 'w') as fout:
        fout.write('\t'.join(['#rank', 'gene', 'heat_score']) + '\n')
        for rank, (gene, score) in enumerate(top_results, 1):
            fout.write(f"{rank}\t{gene}\t{score:.6e}\n")
    
    print(f"\nResults saved to '{outfile}'")
    
    return [(rank, gene, score) for rank, (gene, score) in enumerate(top_results, 1)]


if __name__ == '__main__':
    
    # Check input
    input_list = sys.argv
    network_file, seeds_file, top_n, t, outfile = check_input_style(input_list)
    
    # Read network and seeds
    G_original, seed_genes = read_input(network_file, seeds_file)
    
    # Run Heat Diffusion
    results = HeatDiffusion(G_original, seed_genes, top_n, t, outfile=outfile)
    
    print(f"\nTop 10 predicted genes:")
    for rank, gene, score in results[:10]:
        print(f"  {rank}. {gene}: {score:.6e}")