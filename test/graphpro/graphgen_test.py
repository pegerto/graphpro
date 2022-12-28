import pytest

import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap

from MDAnalysis.tests.datafiles import PDB

u1 = mda.Universe(PDB)

def test_graph_generation_from_mdanalysis():
    G = md_analisys(u1).generate(ContactMap(cutoff=1))
    assert(len(G.graph.nodes) == 214)