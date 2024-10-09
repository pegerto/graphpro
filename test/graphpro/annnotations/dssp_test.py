import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import DSSP

from MDAnalysis.tests.datafiles import PDB_helix

u1 = mda.Universe(PDB_helix)

def test_dssp():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [DSSP()])
    assert len(G.nodes()) == 13
    assert G.node_attr(0)['dssp'] == '-'
    assert G.node_attr(1)['dssp'] == 'H'

def test_encode_enconde():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [DSSP()])
    data = G.to_data(node_encoders=[DSSP()])
    assert data.x.size() == (13,3)