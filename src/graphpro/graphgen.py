import networkx as nx 
from scipy.spatial import distance
from .graph import Graph

class RepresentationMethod():
    def res_map(self, ag):
        pass
        
    def generate(ag):
        pass


class ContactMap(RepresentationMethod):
    def __init__(self, cutoff, chain=None):
        self.cutoff = cutoff
        self.chain = chain
    
    def generate(self, ag):
        ca_position = ag.c_alphas_positions(self.chain)
        dist = distance.squareform(distance.pdist(ca_position))
        dist[dist > self.cutoff] = 0
        return Graph(dist, ca_position, ag.c_alphas_residues())


class ProGraphGenerator:
    def __init__(self, ag):
        self.ag = ag

    def generate(self, rep):
        G = rep.generate(self.ag)
        return G