import MDAnalysis as mda
import numpy as np
import networkx as nx
import os

from graphpro import md_analisys
from graphpro.graph import Graph
from graphpro.graphgen import ContactMap

DISTANCES = np.array([[4,2],[2,4]])
POSITIONS_3D = np.array([[1,2,3], [-1,-2,-3]])
RESIDUES = [{"resid": 1, "resname": "GLY"}, {"resid": 2, "resname": "ALA"},]
SIMPLE_G = Graph(DISTANCES, POSITIONS_3D, RESIDUES)
HETNAM = mda.Universe(os.path.dirname(os.path.realpath(__file__)) + '/../testdata/hetnam.pdb')
HETNAM_G = md_analisys(HETNAM).generate(ContactMap(cutoff=6))

def test_graph():
    assert (len(SIMPLE_G.graph().nodes) == 2)
    assert (SIMPLE_G.node_attr(0)['resname'] == "G")

def test_nx_graph_contains_attributes():
    G = SIMPLE_G.graph()
    assert (nx.get_node_attributes(G,"resname")[0] == "G")

def test_graph_allow_add_node_attributes():
    HETNAM_G.node_attr_add(1,"test", 1)
    assert (HETNAM_G.node_attr(1)['test'] == 1)
    
def test_graph_allow_retrieve_nodes_by_resid():
    assert(HETNAM_G.get_node_by_resid(13920) == 5749)