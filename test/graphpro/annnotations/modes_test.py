import MDAnalysis as mda
import torch

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import GNMSlowModes, ANMSlowModes

from MDAnalysis.tests.datafiles import PDB_small

u1 = mda.Universe(PDB_small)

def test_gnm_slow_modes_annotation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [GNMSlowModes(modes=3)])
    assert len(G.nodes()) == 214
    assert G.node_attr(0)['gnm_slow_mode_0'] == -0.038
    assert G.node_attr(0)['gnm_slow_mode_1'] == -0.047
    assert G.node_attr(0)['gnm_slow_mode_2'] == 0.074

def test_gnm_slow_modes_is_present_for_all_nodes():
     G = md_analisys(u1).generate(ContactMap(cutoff=6), [GNMSlowModes(modes=3)])
     for n in G.nodes():
          assert 'gnm_slow_mode_0' in G.node_attr(n)
          assert 'gnm_slow_mode_1' in G.node_attr(n)

def test_gnm_encoding():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [GNMSlowModes(modes=3)])
    data = G.to_data(node_encoders=[GNMSlowModes(modes=3)])

    assert data.x.size() == (214, 3)
    assert data.x.dtype == torch.float


def test_anm_slow_modes_annotation():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [ANMSlowModes(modes=1)])
    assert len(G.nodes()) == 214
    assert G.node_attr(0)['anm_slow_mode_0'] == 0.014