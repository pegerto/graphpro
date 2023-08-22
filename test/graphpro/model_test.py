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
