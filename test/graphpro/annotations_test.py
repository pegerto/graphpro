import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import ResidueType

from MDAnalysis.tests.datafiles import PDB

u1 = mda.Universe(PDB)

def test_resname_annotation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ResidueType()])
    
    assert len(G.graph().nodes) == 214
    assert G.node_attr(0)['resname'] == 'M'