import numpy as np
from .graph import Graph
from .collection import GraphCollection



class RepresentationMethod():
    """ This interface defines a generation strategy and can be extended to implement new strategies for transforming 
        a collection of atoms into a graph representation.
    """
    def res_map(self, ag, chain=None):
        pass

    def generate(self, ag, name: str):
        pass


class ContactMap(RepresentationMethod):
    """ Illustrates the spatial proximity between amino acids in a protein structure. 
    """
    def __init__(self, cutoff, chain=None):
        self.cutoff = cutoff
        self.chain = chain

    def generate(self, ag, name: str):
        from scipy.spatial import distance
        
        ca_position = ag.c_alphas_positions(self.chain)
        dist = distance.squareform(distance.pdist(ca_position))
        dist[dist > self.cutoff] = 0
        return Graph(name, dist, ca_position, ag.c_alphas_residues(self.chain))

class KNN(RepresentationMethod):
    """ Generate the structure form a defined number of neighbours
    """
    def __init__(self, k, chain=None):
        self.k = k
        self.chain = chain

    def generate(self, ag, name: str):
        from scipy.spatial import  KDTree

        ca_position = ag.c_alphas_positions(self.chain)
        kdtree = KDTree(ca_position)
        residue_num = len(ca_position)
        adjacency = np.zeros((residue_num, residue_num))
        for i, pos in enumerate(ca_position):
            _, neig = kdtree.query(pos, k= self.k + 1)
            for j in neig:
                adjacency[i,j] = 1

        return Graph(name, adjacency, ca_position, ag.c_alphas_residues(self.chain))



class GraphProGenerator:
    """ Graph Pro Generator
        
        Generate both a graph or a graph colection from a structure of a trajectory.
    """
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
