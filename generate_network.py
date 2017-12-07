# author: Roman Schulte-Sasse
# email: sasse@molgen.mpg.de
# generate_network.py


import networkx as nx
import numpy as np
import pandas as pd
import argparse, random
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Generate Random Network with Graph Motifs')
    parser.add_argument('--node_num', help='Number of Nodes in Network',
                        dest='node_num',
                        type=int
                        )
    parser.add_argument('--min_edge_num',
                        help='Number of Edges in Network',
                        dest='min_edge_num',
                        type=int
                        )
    parser.add_argument('--subnetworks',
                        help='Path to File with subnetworks',
                        dest='subnet_path',
                        type=str
                        )
    parser.add_argument('--out',
                        help='Path to output file',
                        dest='out_path',
                        type=str
                        )
    args = parser.parse_args()
    return args.node_num, args.min_edge_num, args.subnet_path, args.out_path

def get_random_neighbor(G, position):
    all_nbs = nx.neighbors(G, position)
    return random.choice(all_nbs)


def _implant_rec(to_assign, G, G_m, pos, mapping):
    if len(to_assign) is 0:
        return mapping
    pos_gm = nodes_sub[0]
    mapping[pos] = pos_gm
    nbs_gm = [i for i in nx.all_neighbors(G_m, pos_gm) if i in to_assign]
    nbs_g = nx.all_neighbors(G, pos)
    # TODO: Check if nbs_gm is longer than nbs_g and add nodes to g until same size
    next_nb = 0
    for node_gm in nbs_gm:
        _implant_rec(to_assign[1:], G, G_m, nbs_g[next_nb], mapping)
        next_nb += 1


def _get_neighbors(G, pos, closed_list):
    return [i for i in nx.all_neighbors(G, pos) if not i in closed_list]

def _add_edges_if_required(G, pos_g, neighbors_G, neighbors_sub, next_node_count):
    changed = False
    if len(neighbors_sub) > len(neighbors_G):
        changed = True
        diff = len(neighbors_sub) - len(neighbors_G)
        for i in range(diff):
            G.add_node(next_node_count)
            G.add_edge(pos_g, next_node_count)
            print ("new edge: {} -> {}".format(pos_g, next_node_count))
            next_node_count += 1
    return G, next_node_count, changed


def calculate_node_mapping(G, G_m, position):
    """Calculates a mapping of nodes between network and subnetwork.

    This method calculates a mapping between nodes in the real graph
    and a subnetwork, given a position at which to implant the subnetwork.
    It adds nodes and edges to the original graph if required for implantation.
    """
    next_node_name = nx.number_of_nodes(G)
    pos_gm = list(G_m.nodes())[0]
    pos_g = position
    print (pos_gm, pos_g)
    mapping = {pos_gm: pos_g}
    to_explore = [pos_gm]
    mapped_g = [pos_g]
    mapped_gm = [pos_gm]
    while not len(to_explore) is 0:
        pos_gm = to_explore[0]
        pos_g = mapping[pos_gm]
        nbs_g = _get_neighbors(G, pos_g, mapped_g)
        nbs_gm = _get_neighbors(G_m, pos_gm, mapped_gm)
        # add node if more neighbors in motif than in G
        G, next_node_name, changed = _add_edges_if_required(G,
                                                            pos_g,
                                                            nbs_g,
                                                            nbs_gm,
                                                            next_node_name
                                                            )
        if changed:
            nbs_g = _get_neighbors(G, pos_g, mapped_g)
            nbs_gm = _get_neighbors(G_m, pos_gm, mapped_gm)
        to_explore.extend(nbs_gm) # add neighbors to open list

        # map neighbors of current node
        for i in range(len(nbs_gm)):
            mapping[nbs_gm[i]] = nbs_g[i]
            mapped_g.append(nbs_g[i])
            mapped_gm.append(nbs_gm[i])
        to_explore.pop(0) # remove from open list
    return G, mapping

def implant_single_motif(G, G_m, position):
    """Implant a motif (G_m) at position in G.

    This method incorporates a graph motif G_m into the graph G at position.
    It is done recursively.
    """
    # construct a mapping from G_m to G
    G, mapping = calculate_node_mapping(G, G_m, position)

    # now, add edges to graph according to subnetwork
    for from_gm, to_gm in G_m.edges():
        if not G.has_edge(mapping[from_gm], mapping[to_gm]):
            G.add_edge(mapping[from_gm], mapping[to_gm])
    return G


def implant_motifs(G, subnetworks):
    """Implant subnetworks in a given graph.

    This method implants subnetworks into a given network.
    This is done by iteratively mapping nodes from the subnetwork
    to nodes from the original graph and adding nodes if required.
    """
    implant_positions = []
    for sub_idx in range(len(subnetworks)):
        position = random.randrange(0, nx.number_of_nodes(G))
        print ("sub_idx: {}\tposition: {}".format(sub_idx, position))
        G = implant_single_motif(G, subnetworks[sub_idx], position)
        implant_positions.append(position)
    return G, implant_positions

def generate_network(node_num, min_edge_num, subnetworks):
    """Generate a biological network with implanted subnetworks.

    This method generates a network that has biological properties
    like following a power law distribution for the node degrees and
    containing hubs and bottlenecks as well as a whole lot of low-degree
    nodes.
    Furthermore, subnetworks will be implanted to the graph randomly, disturbing
    the network structure as little as possible while still containing the
    desired subnetwork topology.
    """
    # Create random network.
    for i in range(1, node_num):
        G = nx.powerlaw_cluster_graph(node_num, i, 0.8, seed=42)
        if nx.number_of_edges(G) > min_edge_num:
            break

    # plot node degree distribution
    degrees = pd.Series(G.degree())
    fig = plt.figure(figsize=(14, 8))
    plt.hist(degrees, bins=np.arange(1, 100, 1))
    fig.savefig('random_network_distribution.png')

    # output some statistics
    print ("Nodes: {}".format(len(G.nodes())))
    print ("Edges: {}".format(len(G.edges())))

    # implant motifs and return
    return implant_motifs(G, subnet_path)

if __name__ == "__main__":
    #node_num, min_edge_num, subnet_path, out_path = parse_args()
    #generate_network(node_num=node_num, min_edge_num=min_edge_num)
    subnetworks = []
    # star
    A1 = np.array([[0,1,1,1],
                   [1,0,0,0],
                   [1,0,0,0],
                   [1,0,0,0]])
    subnetworks.append(nx.from_numpy_matrix(A1))
    # complex
    A2 = np.array([[0,1,0,0,0,0,0,0],
                   [1,0,1,1,0,0,0,0],
                   [0,1,0,1,1,0,0,0],
                   [0,1,1,0,0,0,0,1],
                   [0,0,1,0,0,1,1,0],
                   [0,0,0,0,1,0,0,0],
                   [0,0,0,0,1,0,0,0],
                   [0,0,0,1,0,0,0,0]])
    #subnetworks.append(nx.from_numpy_matrix(A2))
    # clique
    A3 = np.array([[0,1,1,1],
                   [1,0,1,1],
                   [1,1,0,1],
                   [1,1,1,0]])
    subnetworks.append(nx.from_numpy_matrix(A3))

    # network
    A = np.array([[0,1,0,0,0,0,0,0],
                  [1,0,1,1,0,0,0,0],
                  [0,1,0,1,1,0,0,0],
                  [0,1,1,0,0,0,0,1],
                  [0,0,1,0,0,1,1,0],
                  [0,0,0,0,1,0,0,0],
                  [0,0,0,0,1,0,0,0],
                  [0,0,0,1,0,0,0,0]])
    G, positions = implant_motifs(nx.from_numpy_matrix(A), subnetworks)
    print (G.nodes())
    print (G.edges())
    print (positions)
    #nx.draw(G)
