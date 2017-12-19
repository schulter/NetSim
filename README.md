# NetSim
This package will simulate realistic biological networks and implant subnetworks into them.
NetSim consists of two parts that let you create realistic biological networks with implanted subnetworks like cancer modules or similar. It uses the Barabasi-Albert graph to create a network and then inserts subnetworks into it.
NetSim offers two insertion strategies (strategies to find nodes at which the subnetworks are inserted):
1. Random
2. PageRank (the subnetwork will be inserted at positions that have the highest PageRank indices in the network)

## Usage
To simulate your network, you can run the script or access NetSim programmatically.
### Run Script
On the command line, you can simply type:

```python generate_network.py --node_num <#nodes> --min_edge_num <#edges> --subnetpath <path> --outdir <path>```

*Note: The number of edges can not be determined exactly by the Barabasi-Albert-Graph. So the graph is run with increasing numbers of neighbors for any new node until at least `min_edge_num` edges are in the graph.*

* `--node_num`: The number of nodes in the simulated network
* `--min_edge_num`: The minimum number of edges in the simulated network
* `--subnetpath`: The directory in which the subnetworks are contained. From that folder, each file with the corresponding file ending will be treated as valid subnetwork.
* `--outdir`: The directory to which output is written. Output is the network is edgelist format as well as a textfile containing the insertion positions for the subnetworks and plots of node degree distribution and shortest path distribution in the network.

### Programmatic Access
You can also directly use the code as follows:

```python
from generate_network import NetworkGenerator
import numpy as np
import networkx as nx
simulator = NetworkGenerator(num_nodes=10, min_num_edges=10)
# some simple clique subnetwork
A = np.array([[0,1,1,1],
              [1,0,1,1],
              [1,1,0,1],
              [1,1,1,0]]) # clique
subnetworks = [nx.from_numpy_matrix(A)]
# generate network
G, insert_pos = simulator.generate_network(subnetworks)
# plot
simulator.plot_distributions('.')
simulator.save_network('.')
simulator.draw_network('.')
```
## Requirements
This package is implemented in python using the following packages:
* `networkx`
* `pandas`
* `numpy`
