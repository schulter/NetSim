# author: Roman Schulte-Sasse
# email: sasse@molgen.mpg.de
# generate_network.py


import networkx as nx
import numpy as np
import pandas as pd
import argparse, random, os, operator
import matplotlib.pyplot as plt


class NetworkGenerator:
    """Easily simulate biological networks with implanted subnetworks.

    This class provides an easy-to-use interface to simulate large networks
    and to implant subnetworks into them.
    """
    def __init__(self, graph=None, num_nodes=None, min_num_edges=None,
                 insert_strategy='random'):
        # graph already specified
        if not graph is None:
            self.graph = graph
            self.number_of_nodes = nx.number_of_nodes(graph)
            self.minimum_number_of_edges = nx.number_of_edges(graph)
        # graph is to be generated
        elif not num_nodes is None and not min_num_edges is None:
            self.number_of_nodes = num_nodes
            self.minimum_number_of_edges = min_num_edges
            self.graph = None
        self.supported_formats = ['edgelist', 'gml']
        self.mappings = None
        self.insert_strategy_name = insert_strategy
        if insert_strategy is 'random':
            self.calculate_insert_positions = self.random_insert_strategy
        elif insert_strategy is 'pagerank':
            self.calculate_insert_positions = self.pagerank_insert_strategy
        else:
            raise Exception("Unknown insert strategy {}".format(insert_strategy))

    def random_insert_strategy(self, G):
        return random.sample(range(0, self.number_of_nodes), self.number_of_nodes)

    def pagerank_insert_strategy(self, G):
        scores = nx.pagerank(G)
        sorted_scores = sorted(scores.items(), key=operator.itemgetter(1))[::-1]
        return [k for k, v in sorted_scores]

    def _get_neighbors(self, G, pos, closed_list):
        return [i for i in nx.all_neighbors(G, pos) if not i in closed_list]

    def _insert_strategy(self, G, num_of_inserts, min_distance=2):
        # get insert positions (as many as nodes in network)
        # TODO: This method does not have any failsafe if we run out of positions
        scores = self.calculate_insert_positions(G)

        # add highest ranks as insert positions
        insert_positions = []
        found = 0
        attempts = 0
        while num_of_inserts > found:
            insert_pos = scores[found + attempts]
            allowed = True
            for pos in insert_positions:
                if nx.shortest_path_length(G, insert_pos, pos) < min_distance:
                    allowed = False
                attempts += 1
            if allowed:
                insert_positions.append(insert_pos)
                found += 1
                attempts = 0
        assert (len(insert_positions) == num_of_inserts)
        return insert_positions

    def _add_edges_if_required(self, G, pos_g, neighbors_G, neighbors_sub,
                               next_node_count):
        changed = False
        if len(neighbors_sub) > len(neighbors_G):
            changed = True
            diff = len(neighbors_sub) - len(neighbors_G)
            for i in range(diff):
                G.add_node(next_node_count)
                G.add_edge(pos_g, next_node_count)
                #print ("new edge: {} -> {}".format(pos_g, next_node_count))
                next_node_count += 1
        return G, next_node_count, changed

    def calculate_node_mapping(self, G, G_m, position):
        """Calculates a mapping of nodes between network and subnetwork.

        This method calculates a mapping between nodes in the real graph
        and a subnetwork, given a position at which to implant the subnetwork.
        It adds nodes and edges to the original graph if required for implantation.
        """
        next_node_name = nx.number_of_nodes(G)
        pos_gm = list(G_m.nodes())[0]
        pos_g = position
        mapping = {pos_gm: pos_g}
        to_explore = [pos_gm]
        mapped_g = [pos_g]
        mapped_gm = [pos_gm]
        while not len(to_explore) is 0:
            pos_gm = to_explore[0]
            pos_g = mapping[pos_gm]
            nbs_g = self._get_neighbors(G, pos_g, mapped_g)
            nbs_gm = self._get_neighbors(G_m, pos_gm, mapped_gm)
            # add node if more neighbors in motif than in G
            G, next_node_name, changed = self._add_edges_if_required(G,
                                                                     pos_g,
                                                                     nbs_g,
                                                                     nbs_gm,
                                                                     next_node_name
                                                                     )
            if changed:
                nbs_g = self._get_neighbors(G, pos_g, mapped_g)
                nbs_gm = self._get_neighbors(G_m, pos_gm, mapped_gm)
            to_explore.extend(nbs_gm) # add neighbors to open list

            # map neighbors of current node
            for i in range(len(nbs_gm)):
                mapping[nbs_gm[i]] = nbs_g[i]
                mapped_g.append(nbs_g[i])
                mapped_gm.append(nbs_gm[i])
            to_explore.pop(0) # remove from open list
        return G, mapping

    def _implant_single_motif(self, G, G_m, position):
        """Implant a motif (G_m) at position in G.

        This method incorporates a graph motif G_m into the graph G at position.
        It is done recursively.
        """
        # construct a mapping from G_m to G
        G, mapping = self.calculate_node_mapping(G, G_m, position)
        # now, add edges to graph according to subnetwork
        for from_gm, to_gm in G_m.edges():
            if not G.has_edge(mapping[from_gm], mapping[to_gm]):
                G.add_edge(mapping[from_gm], mapping[to_gm])
        return G, mapping

    def implant_motifs(self, G, subnetworks):
        """Implant subnetworks in a given graph.

        This method implants subnetworks into a given network.
        This is done by iteratively mapping nodes from the subnetwork
        to nodes from the original graph and adding nodes if required.
        """
        implant_positions = self._insert_strategy(G, len(subnetworks))
        mappings = []
        for sub_idx in range(len(subnetworks)):
            position = implant_positions[sub_idx]
            print ("Implanting network {} at position: {}".format(sub_idx, position))
            G, mapping = self._implant_single_motif(G,
                                                    subnetworks[sub_idx],
                                                    position)
            mappings.append(mapping)
        return G, mappings

    def plot_distributions(self, out_dir='.'):
        """Plot the distribution of node degrees and shortest paths.
        """
        # plot node degree distribution
        degrees = np.array([v for k, v in list(self.graph.degree())])
        fig = plt.figure(figsize=(14, 8))
        plt.hist(degrees, bins=np.arange(0, degrees.max(), 1))
        plt.xlabel('Node Degree')
        plt.ylabel('Frequency')
        plt.title('Node Degree Distribution')
        fig.savefig(os.path.join(out_dir, 'network_node_distribution.png'))

        # plot distribution of shortest paths
        paths = pd.DataFrame(dict(nx.shortest_path_length(self.graph)))
        fig = plt.figure(figsize=(14, 8))
        plt.hist(paths.values.flatten(), bins=np.arange(0, 8, 1))
        plt.xlabel('Shortest Path Length')
        plt.ylabel('Frequency')
        plt.title('Distribution of Distance Between Nodes')
        fig.savefig(os.path.join(out_dir, 'network_paths_distribution.png'))

    def generate_network(self, subnetworks):
        """Generate a biological network with implanted subnetworks.

        This method generates a network that has biological properties
        like following a power law distribution for the node degrees and
        containing hubs and bottlenecks as well as a whole lot of low-degree
        nodes.
        Furthermore, subnetworks will be implanted to the graph randomly, disturbing
        the network structure as little as possible while still containing the
        desired subnetwork topology.

        Parameters:
        -----------
        subnetworks:            The subnetworks to insert into the network.

        Returns:
        A simulated network with biological properties and implanted subnetworks
        as well as as list of the positions at which insertion happened.
        """
        # Create random network.
        print ("Creating Network...")
        for i in range(1, self.number_of_nodes):
            G = nx.barabasi_albert_graph(n=self.number_of_nodes, m=i)
            if nx.number_of_edges(G) > self.minimum_number_of_edges:
                break

        # implant motifs and return
        print ("Network built! Implanting {} subnetworks using {} strategy...".format(len(subnetworks),
                                                                                      self.insert_strategy_name))
        G, mappings = self.implant_motifs(G, subnetworks)

        # output some statistics and return
        print ("Done! Created network with {} nodes and {} edges".format(nx.number_of_nodes(G),
                                                                         nx.number_of_edges(G))
               )
        self.graph = G
        self.mappings = mappings
        return G, mappings

    def read_subnetworks(self, subnet_dir, f_format='edgelist'):
        """Read the subnetworks to implant from disk.

        Supported formats are edgelist and gml for the moment.
        The method traverses the given directory and searches for files with
        the endings corresponding to the file format. Once, such a file is
        found, it is read to a networkx graph object.
        Parameters:
        ----------
        subnet_dir:             The directory in which the subnetworks are
                                located. Each file corresponds to one subnetwork
        f_format:               The file format of the subnetworks.

        Returns:
        A list with the subnetworks.
        """
        # check if file format is supported
        if not f_format in self.supported_formats:
            print ("{} is not supported. Please provide your subnetworks as one of {}".format(f_format, self.supported_formats))

        # read subnetworks and return them
        subnetworks = []
        for files in os.listdir(subnet_dir):
            f_path = os.path.join(subnet_dir, files)
            if files.endswith(f_format):
                if f_format is 'edgelist':
                    print ("Found valid edgelist subnetwork in {}".format(f_path))
                    subnetworks.append(nx.read_edgelist(f_path))
                elif f_format is 'gml':
                    print ("Found valid gml subnetwork in {}".format(f_path))
                    subnetworks.append(nx.read_gml(f_path))
        return subnetworks


    def draw_network(self, out_dir):
        if not self.graph is None:
            fig = plt.figure(figsize=(14, 8))
            nx.draw(self.graph, with_labels=True)
            fig.savefig(os.path.join(out_dir, 'network.png'))


    def save_network(self, out_dir, f_format='edgelist'):
        """Write the generated network to file.
        """
        # check that file format is supported
        if not f_format in self.supported_formats:
            print ("{} is not a supported file format. Try one of {}".format(f_format, self.supported_formats))
            return

        if f_format is 'edgelist':
            nx.write_edgelist(self.graph, os.path.join(out_dir, 'network.edgelist'))
        elif f_format is 'gml':
            nx.write_gml(self.graph, os.path.join(out_dir, 'network.gml'))

        if not self.mappings is None:
            with open(os.path.join(out_dir, 'implant_positions.txt'), 'w') as f:
                count = 1
                f.write('# subnetwork insert positions. Each row corresponds to a subnetwork and each column to a position in it.\n')
                for m in self.mappings:
                    f.write('Subnetwork {}: '.format(count))
                    for fro, to in m.items():
                        f.write('{}\t'.format(to))
                    f.write('\n')
                    count += 1


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
    parser.add_argument('--subnetpath',
                        help='Path to dir with subnetworks',
                        dest='subnetpath',
                        default=None,
                        type=str
                        )
    parser.add_argument('--insert_strategy',
                        help='How to insert the subnetworks',
                        dest='insert_strategy',
                        default='random',
                        type=str
                        )
    parser.add_argument('--outdir',
                        help='Path to output directory (writes network and plots)',
                        dest='outdir',
                        type=str
                        )
    args = parser.parse_args()
    return args.node_num, args.min_edge_num, args.subnetpath, args.insert_strategy, args.outdir


if __name__ == "__main__":
    node_num, min_edge_num, subnet_dir, insert_strat, out_dir = parse_args()
    simulator = NetworkGenerator(num_nodes=node_num,
                                 min_num_edges=min_edge_num,
                                 insert_strategy=insert_strat)
    if not subnet_dir is None:
        subnetworks = simulator.read_subnetworks(subnet_dir)
    else:
        A = np.array([[0,1,1,1,1],
                      [1,0,1,1,1],
                      [1,1,0,1,1],
                      [1,1,1,0,1],
                      [1,1,1,1,0]]) # clique
        subnetworks = [nx.from_numpy_matrix(A)] * 30
    G, insert_pos = simulator.generate_network(subnetworks)
    if not out_dir is None:
        simulator.plot_distributions(out_dir)
        simulator.save_network(out_dir)
    else:
        simulator.plot_distributions()
        simulator.save_network('.')
    if node_num < 300: # draw when not so many nodes
        simulator.draw_network('.')
