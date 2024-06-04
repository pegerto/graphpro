import torch
import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import ResidueType, NodeTargetBinaryAttribute, SASAResArea, Polarity

from MDAnalysis.tests.datafiles import PDB_small

u1 = mda.Universe(PDB_small)

def test_node_target_binary():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ResidueType()])
    G.node_attr_add(1, "test_attr", "x")

    target = NodeTargetBinaryAttribute("test_attr")
    y = target.encode(G)
    assert y.size() == (214,2)
    assert torch.all(y[0].eq(torch.tensor([1., 0.])))
    assert torch.all(y[1].eq(torch.tensor([0., 1.])))


def test_resname_annotation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ResidueType()])

    assert len(G.nodes()) == 214
    assert G.node_attr(0)['resname'] == 'M'

def test_resname_encoded():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ResidueType()])
    data = G.to_data(node_encoders=[ResidueType()])

    assert data.x.size() == (214, 22)
    assert data.x.dtype == torch.float

def test_sasa_annotation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [SASAResArea(chain_id='')])
    assert G.node_attr(0)['sasa_area'] != 0

def test_sasa_encoding():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [SASAResArea(chain_id='')])
    data = G.to_data(node_encoders=[ SASAResArea(chain_id='')])
    
    assert data.x.size() == (214,1)

def test_polarity_generation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [Polarity(), ResidueType()])
    assert G.node_attr(0)['resname'] == 'M'
    assert G.node_attr(0)['polarity'] == 'a'
    assert G.node_attr(30)['resname'] == 'T'
    assert G.node_attr(30)['polarity'] == 'p'

def test_polarity_encoding():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [Polarity()])
    data = G.to_data(node_encoders=[Polarity()])
    
    assert data.x.size() == (214, 5)
    
def test_encode_multiple_attributes():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ResidueType(), SASAResArea(chain_id='')])
    data = G.to_data(node_encoders=[ResidueType(), SASAResArea(chain_id='')])

    assert data.x.size() == (214,23)