import numpy as np

from graphpro.graph import Graph
from graphpro.collection import GraphCollection

DISTANCES = np.array([[4,2],[2,4]])
POSITIONS_3D = np.array([[1,2,3], [-1,-2,-3]])
RESIDUES = [{"resid": 1, "resname": "GLY"}, {"resid": 2, "resname": "ALA"},]
SIMPLE_G = Graph("test", DISTANCES, POSITIONS_3D, RESIDUES)

def test_graph_collection():
    col = GraphCollection([SIMPLE_G])
    assert(len(col) == 1)