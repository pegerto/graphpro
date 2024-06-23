import networkx as nx
import numpy as np
import torch

from torch_geometric.data import Data
from dataclasses import dataclass
from .model import Target, ProteinMetadata

class Graph():
    """ Graph provides a representation of a graph and required helpers.
    """

    def __init__(self, name: str,
                 distances: np.array,
                 positions: np.array,
                 res_map: dict[int, dict],
                 metadata: ProteinMetadata = None):

        self.name = name
        self.distances = distances
        self.positions = positions
        self._n_attr = {i: {"resid": res_attr['resid']}
                        for i, res_attr in enumerate(res_map)}
        self._resid_to_node = {
            res_attr['resid']: i for i,
            res_attr in enumerate(res_map)}
        self.metadata = metadata

    def __eq__(self, other):
        """Compare two graphs for equality"""
        # TODO: may need to compare more than the distances
        if not other:
            return False
        return (self.distances == other.distances).any()

    def node_attr(self, node_id: int):
        return self._n_attr.get(node_id)

    def node_attr_add(self, node_id: int, attribute_name: str, attribute: any):
        """Adds a specific attribute to a noode in the graph"""
        attrs = self.node_attr(node_id)
        if attrs:
            attrs[attribute_name] = attribute

    def get_node_by_resid(self, resid: int) -> int:
        """Returns the node number using the residue id, None if the residue id is not known"""
        return self._resid_to_node.get(resid)

    def communities(self) -> list[tuple[float, list[set]]]:
        """ Perform Girvan Newman communinity detection returning the list of communities.
            The algorithm is perform all the way until no more edges are left to be removed.

        """
        from networkx.algorithms import community

        c_iter = community.girvan_newman(self.to_networkx())
        return [(community.modularity(self.to_networkx(), com), com)
                for com in c_iter]

    def to_networkx(self) -> nx.Graph:
        """ Returns a networkx G undirected graph with populated attributes """
        G = nx.from_numpy_array(self.distances)
        nx.set_node_attributes(G, self._n_attr)
        return G

    def to_data(self, node_encoders = [], target: Target = None) -> Data:
        """ Return a PyG object from this existing graph"""
        directed = torch.tensor([[edge[0],edge[1]] for edge in self.to_networkx().edges], dtype=torch.long)
        inversed = torch.tensor([[edge[1],edge[0]] for edge in self.to_networkx().edges], dtype=torch.long)
        cco = torch.cat((directed, inversed), 0).t().contiguous()
    
        x = None
        y = None

        # Concat a list of node features into a X tensor
        for encoder in node_encoders:
            ecoded_attr = encoder.encode(self)
            if isinstance(x,torch.Tensor):
                x = torch.concat((x, ecoded_attr), 1) # concat to 1 dim
            else:
                x = ecoded_attr
        if target:
            y = target.encode(self)

        return Data(x=x, edge_index=cco, y=y)

    def nodes(self) -> list[int]:
        """ Return node list """
        return self._resid_to_node.values()
    
    def plot(self,
             figsize: tuple[int, int] = (8, 10),
             communities: list[set[int]] = [],
             show = True
             ) -> None:
        """ Plot the graph represention in 3D using residue positions.
        """
        import matplotlib.pyplot as plt

        node_xyz = np.array([self.positions[v] for v in sorted(self.to_networkx())])
        edge_xyz = np.array([(self.positions[u], self.positions[v])
                            for u, v in self.to_networkx().edges()])

        node_colors = None
        if len(communities) > 0:
            community_node = sorted(
                [(n, i) for i, c in enumerate(communities) for n in c])
            node_colors = [n[1] for n in community_node]

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(*node_xyz.T, s=100, ec="w", c=node_colors)

        # Plot the edges
        for vizedge in edge_xyz:
            ax.plot(*vizedge.T, color="tab:gray")

        def _format_axes(ax):
            ax.grid(False)
            for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
                dim.set_ticks([])

        _format_axes(ax)
        fig.tight_layout()
        if show:
            plt.show()

    def __repr__(self) -> str:
        return self.name
