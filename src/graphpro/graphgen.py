import networkx as nx 
from scipy.spatial import distance
from .graph import Graph

class RepresentationMethod():
    def res_map(self, ag):
        pass
        
    def generate(self, ag, name: str):
        pass


class ContactMap(RepresentationMethod):
    def __init__(self, cutoff, chain=None):
        self.cutoff = cutoff
        self.chain = chain
    
    def generate(self, ag, name: str):
        ca_position = ag.c_alphas_positions(self.chain)
        dist = distance.squareform(distance.pdist(ca_position))
        dist[dist > self.cutoff] = 0
        return Graph(name, dist, ca_position, ag.c_alphas_residues(self.chain))


class ProGraphGenerator:
    def __init__(self, ag, name=''):
        self.ag = ag
        self.name = name

    def generate(self, rep):
        G = rep.generate(self.ag, self.name)
        return G