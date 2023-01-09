import MDAnalysis as mda
import os

from graphpro import md_analisys
from graphpro.graphgen import ContactMap

HETNAM = mda.Universe(os.path.dirname(os.path.realpath(__file__)) + '/../testdata/hetnam.pdb')
G = md_analisys(HETNAM).generate(ContactMap(cutoff=6))

def test_graph_allow_add_node_attributes():
    G.node_attr_add(1,"test", 1)
    assert (G.node_attr(1)['test'] == 1)
    
def test_graph_allow_retrieve_nodes_by_resid():
    assert(G.get_node_by_resid(13920) == 5749)