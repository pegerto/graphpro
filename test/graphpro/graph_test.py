import os

import MDAnalysis as mda
import numpy as np
import networkx as nx
import torch 

from torch_geometric.data import Data

from graphpro import md_analisys
from graphpro.graph import Graph
from graphpro.graphgen import ContactMap
from graphpro.annotations import ResidueType

DISTANCES = np.array([[4, 2], [2, 4]])
POSITIONS_3D = np.array([[1, 2, 3], [-1, -2, -3]])
RESIDUES = [{"resid": 1, "resname": "GLY"}, {"resid": 2, "resname": "ALA"},]
SIMPLE_G = Graph("test", DISTANCES, POSITIONS_3D, RESIDUES)
HETNAM = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../testdata/hetnam.pdb')
HETNAM_G = md_analisys(HETNAM).generate(ContactMap(cutoff=6), [ResidueType()])


def test_graph():
    assert (len(SIMPLE_G.nodes()) == 2)


def test_graph_allow_add_node_attributes():
    HETNAM_G.node_attr_add(1, "test", 1)
    assert (HETNAM_G.node_attr(1)['test'] == 1)


def test_graph_allow_retrieve_nodes_by_resid():
    assert(HETNAM_G.get_node_by_resid(13920) == 5749)


def test_graph_str_informed():
    assert(str(SIMPLE_G) == "test")

def test_graph_plot():
    SIMPLE_G.plot(show=False)

def test_to_data_index():
    edge_index = torch.tensor([[0, 0, 1, 0, 1, 1],
                  [0, 1, 1, 0, 0, 1]], dtype=torch.int32)
    
    assert(torch.allclose(SIMPLE_G.to_data().edge_index,  edge_index))
    