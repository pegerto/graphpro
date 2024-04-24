import os
import MDAnalysis as mda

from graphpro.model import MDAnalisysAtomGroup

FIVEHTC = mda.Universe(
    os.path.dirname(
        os.path.realpath(__file__)) +
    '/../testdata/5htc.pdb')


def test_model_from_mdanalysis_partial_occupancy():
    """ This protein has partial ocupancy on chain C """
    ag = MDAnalisysAtomGroup(FIVEHTC)
    assert(len(ag.c_alphas_residues(chain="C")) ==
           len(ag.c_alphas_positions(chain="C")))


def test_transform_to_prody_contains_same_atoms():
    ag = MDAnalisysAtomGroup(FIVEHTC)
    pdAG = ag.to_prody()
    assert(len(pdAG) == 3346)

def test_transform_to_prody_cotains_atom_names():
    ag = MDAnalisysAtomGroup(FIVEHTC)
    pdAG = ag.to_prody()
    assert(pdAG.getNames()[0] == 'N')

def test_transform_to_prody_cotains_residue_names():
    ag = MDAnalisysAtomGroup(FIVEHTC)
    pdAG = ag.to_prody()
    assert(pdAG.getResnames()[0] == 'GLY')

def test_transform_to_prody_cotains_residue_number():
    ag = MDAnalisysAtomGroup(FIVEHTC)
    pdAG = ag.to_prody()
    assert(pdAG.getResnums()[0] == 471)

def test_transform_to_prody_contains_chain_id():
    ag = MDAnalisysAtomGroup(FIVEHTC)
    pdAG = ag.to_prody()
    assert(pdAG.getChids()[0] == 'A')