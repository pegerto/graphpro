import os
import MDAnalysis as mda

from graphpro import md_analisys
from graphpro.graphgen import ContactMap
from graphpro.annotations import DSSP

from MDAnalysis.tests.datafiles import PDB_helix

u1 = mda.Universe(PDB_helix)
u2 = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../../testdata/4aw0.pdb')

u3 = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../../testdata/1w96.pdb')

u4 = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../../testdata/1qp0.pdb')

u5 = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../../testdata/4a07.pdb')

def test_dssp():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [DSSP()])
    assert len(G.nodes()) == 13
    assert G.node_attr(0)['dssp'] == '-'
    assert G.node_attr(1)['dssp'] == 'H'
    assert G.node_attr(12)['dssp'] == '-'

def test_dssp_multiple_occupancies():
    G = md_analisys(u2).generate(ContactMap(cutoff=6), [DSSP()])
    assert len(G.nodes()) == 283

def test_dssp_incompleted_resiude():
    G = md_analisys(u3).generate(ContactMap(cutoff=6), [DSSP()])
    assert len(G.nodes()) == 554

def test_dssp_enconde():
    G = md_analisys(u1).generate(ContactMap(cutoff=6), [DSSP()])
    data = G.to_data(node_encoders=[DSSP()])
    assert data.x.size() == (13, 4)

def test_dssp_encode_incompleted():
    G = md_analisys(u3).generate(ContactMap(cutoff=6), [DSSP()])
    data = G.to_data(node_encoders=[DSSP()])
    assert data.x.size() == (554, 4)

def test_dssp_encode_incompleted_residue_atoms():
    G = md_analisys(u4).generate(ContactMap(cutoff=6), [DSSP()])
    data = G.to_data(node_encoders=[DSSP()])
    assert data.x.size() == (338, 4)
    
    
def test_dssp_encode_partial_atom_residue_atoms():
    G = md_analisys(u5).generate(ContactMap(cutoff=6), [DSSP()])
    data = G.to_data(node_encoders=[DSSP()])
    assert data.x.size() == (282, 4)