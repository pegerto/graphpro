import networkx as nx 
import numpy as np

from .util.residues import one_letter_res

class Graph():
    """ Graph provides a representation of a graph and required helpers.
    """
    def __init__(self, distances: np.array , positions: np.array, res_map: dict[int, dict]):
        self.distances = distances
        self.positions = positions
        self._n_attr = {i: {"resid": res_attr['resid'], "resname": one_letter_res(res_attr['resname'])} for i, res_attr in enumerate(res_map)}
        self._resid_to_node = {res_attr['resid']: i for i, res_attr in enumerate(res_map)} 

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

        c_iter = community.girvan_newman(self.graph())
        return  [(community.modularity(self.graph(), com), com) for com in c_iter]

    def graph(self) -> nx.Graph:
        """ Returns a networkx G undirected graph with populated attributes """
        G = nx.from_numpy_matrix(self.distances)
        nx.set_node_attributes(G,self._n_attr)
        return G
    
    def plot(self, 
        figsize: tuple[int,int] = (8,10),
        communities: list[set[int]] = []
        ) -> None: 
        """ Plot the graph represention in 3D using real 3D possitions.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        node_xyz = np.array([self.positions[v] for v in sorted(self.graph())])
        edge_xyz = np.array([(self.positions[u], self.positions[v]) for u, v in self.graph().edges()])

        node_colors = None
        if len(communities) > 0: 
            community_node = [(n,i) for i,c in enumerate(communities) for n in c]
            community_node.sort()
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
        plt.show()