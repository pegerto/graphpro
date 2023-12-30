from .graph import Graph
from .collection import GraphCollection


class RepresentationMethod():
    def res_map(self, ag, chain=None):
        pass

    def generate(self, ag, name: str):
        pass


class ContactMap(RepresentationMethod):
    def __init__(self, cutoff, chain=None):
        self.cutoff = cutoff
        self.chain = chain

    def generate(self, ag, name: str):
        from scipy.spatial import distance
        
        ca_position = ag.c_alphas_positions(self.chain)
        dist = distance.squareform(distance.pdist(ca_position))
        dist[dist > self.cutoff] = 0
        return Graph(name, dist, ca_position, ag.c_alphas_residues(self.chain))


class ProGraphGenerator:
    def __init__(self, ag, trajectory = None, name=''):
        self.ag = ag
        self.trajectory = trajectory
        self.name = name

    def _generate(self, ag, rep, node_annotations=[]):
        G = rep.generate(self.ag, self.name)
        if len(G.nodes()) > 0:
            for node_annotation in node_annotations:
                node_annotation.generate(G, self.ag)
        return G
    
    def generate(self, rep, node_annotations=[]) -> Graph:
        return self._generate(self.ag, rep, node_annotations)
    
    def generate_trajectory(self, rep, node_annotations=[]):
        return GraphCollection([self._generate(ag, rep, node_annotations) for ag in self.trajectory])
