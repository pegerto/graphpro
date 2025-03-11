import torch
import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import BTPotential, BTEigenCentrality

from MDAnalysis.tests.datafiles import PDB_small

u1 = mda.Universe(PDB_small)

def test_bt_potential_calculation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [BTPotential()])
    assert len(G.nodes()) == 214
    assert round(G.node_attr(0)['bt_potential'],2) == -702.08

def test_gnm_encoding():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [BTPotential()])
    data = G.to_data(node_encoders=[BTPotential()])

    assert data.x.size() == (214, 1)
    assert data.x.dtype == torch.float
    
def test_bt_potential_calculation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [BTEigenCentrality()])
    assert len(G.nodes()) == 214
    assert round(G.node_attr(90)['bt_eigen_centrality'],5) == 0.00022
