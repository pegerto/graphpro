import torch
import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import ResidueType, NodeTargetBinaryAttribute

from MDAnalysis.tests.datafiles import PDB

u1 = mda.Universe(PDB)

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
    data = G.to_data(node_encoders=  [ResidueType()])

    assert data.x.size() == (214, 22)
    assert data.x.dtype == torch.float