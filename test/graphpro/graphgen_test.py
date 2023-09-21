import os

import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.collection import GraphCollection

from MDAnalysis.tests.datafiles import PDB, XTC

u1 = mda.Universe(PDB, XTC)


def test_graph_generation_from_mdanalysis():
    G = md_analisys(u1).generate(ContactMap(cutoff=6))
    assert(len(G.nodes()) == 214)


def test_graph_generation_from_mdanalysis_custom_residue():
    hetnam = mda.Universe(
        os.path.dirname(
            os.path.realpath(__file__)) +
        '/../testdata/hetnam.pdb')
    G = md_analisys(hetnam).generate(ContactMap(cutoff=6))
    assert(len(G.nodes()) == 5752)


def test_graph_generation_collection():
    graph_col = md_analisys(u1).generate_trajectory(ContactMap(cutoff=6))
    
    assert type(graph_col) == GraphCollection
    assert len(graph_col) == 10 # XTC has 10 frames