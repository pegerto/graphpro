import os
import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import ConservationScore


u3aw0 = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../../testdata/4aw0.pdb')


def test_conservation():
    G = md_analisys(u3aw0, '4AW0').generate(ContactMap(cutoff=6, chain='A'), [ConservationScore()])
    assert len(G.nodes()) == 283
    print(G.node_attr(0))
    print(G.node_attr(1))
    print(G.node_attr(2))
    assert G.node_attr(0)['cons_shannon'] == 0.8285890860167484
    assert G.node_attr(1)['cons_shannon'] == 0.8487017361472784
    assert G.node_attr(2)['cons_shannon'] == 1.330480269847359

def test_conservation_encoding():
    G = md_analisys(u3aw0, '4AW0').generate(ContactMap(cutoff=6, chain='A'), [ConservationScore()])
    data = G.to_data(node_encoders=[ConservationScore()])
    assert data.x.size() == (283, 1)