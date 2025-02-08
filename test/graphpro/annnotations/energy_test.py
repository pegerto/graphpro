import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import BT_Potential

from MDAnalysis.tests.datafiles import PDB_small

u1 = mda.Universe(PDB_small)

def test_bt_potential_calculation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [BT_Potential()])
    assert len(G.nodes()) == 214
    assert G.node_attr(0)['bt_potential'] == -702.0824887338479