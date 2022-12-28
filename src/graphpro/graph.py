import networkx as nx 
import numpy as np

class Graph():
    """ Graph provides a representation of a graph and required helpers.
    """
    def __init__(self, graph: nx.Graph , positions: np.array, res_map: dict[int, int]):
        self.graph = graph
        self.positions = positions
        self.res_map = res_map
    
    
    def communities(self) -> list[tuple[float, list[set]]]:
        """ Perform Girvan Newman communinity detection returning the list of communities.
            The algorithm is perform all the way until no more edges are left to be removed.

        """
        from networkx.algorithms import community

        c_iter = community.girvan_newman(self.graph)
        return  [(community.modularity(self.graph, com), com) for com in c_iter]

    
    def plot(self, 
        figsize: tuple[int,int] = (8,10),
        communities: list[set[int]] = []
        ) -> None: 
        """ Plot the graph represention in 3D using real 3D possitions.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        node_xyz = np.array([self.positions[v] for v in sorted(self.graph)])
        edge_xyz = np.array([(self.positions[u], self.positions[v]) for u, v in self.graph.edges()])

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