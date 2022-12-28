import networkx as nx 
import numpy as np

from scipy.spatial import distance
from .model import AtomGroup
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
    
    def res_map(self, ag):
        return {i: e for i, e in enumerate(ag.c_alphas_residues())}
    
    def generate(self, ag):
        ca_position = ag.c_alphas_positions(self.chain)
        dist = distance.squareform(distance.pdist(ca_position))
        dist[dist > self.cutoff] = 0
        return Graph(nx.from_numpy_matrix(dist), ca_position, self.res_map(ag))


class ProGraphGenerator:
    def __init__(self, ag):
        self.ag = ag

    def generate(self, rep):
        G = rep.generate(self.ag)
        return G