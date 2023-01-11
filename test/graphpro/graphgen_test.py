import pytest
import os 

import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap

from MDAnalysis.tests.datafiles import PDB

u1 = mda.Universe(PDB)

def test_graph_generation_from_mdanalysis():
    G = md_analisys(u1).generate(ContactMap(cutoff=6))
    assert(len(G.graph().nodes) == 214)

def test_graph_generation_from_mdanalysis_custom_residue():
    hetnam = mda.Universe(os.path.dirname(os.path.realpath(__file__)) + '/../testdata/hetnam.pdb')
    G = md_analisys(hetnam).generate(ContactMap(cutoff=6))
    assert(len(G.graph().nodes) == 5752)
    assert(G.node_attr(5751)['resname'] == 'X')