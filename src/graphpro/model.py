import torch
import prody as pdy
import MDAnalysis as mda

from functools import reduce
from dataclasses import dataclass

class Target():
    def __init__(self):
        pass

    def encode(self, G) -> torch.FloatTensor:
        pass

class NodeTarget(Target):
    def __init__(self, attr_name: str):
        self.attr_name = attr_name
    
    def encode(self, G) -> torch.FloatTensor:
        pass

@dataclass
class ProteinMetadata:
    """ Additional information associated to the protein
    """
    uniprot_id: str
    chain: str

class Trajectory():
    def __init__(self):
        pass

class AtomGroup():
    def n_atoms(self) -> int:
        pass

    def c_alphas_positions(self, chain: None):
        pass

    def c_alphas_residues(self, chain: None):
        pass

    def positions(self):
        pass

    def names(self):
        """ Return atoms names
        """
        pass

    def resnames(self):
        """ Returns resnames
        """
        pass

    def resnums(self):
        """ Returns resnumbers
        """
        pass

    def chain_ids(self):
        """ Returns the chainId for each individual atom
        """
        pass

    def to_prody(self) -> pdy.AtomGroup:
        ag = pdy.AtomGroup('name')
        ag.setCoords(self.positions())
        ag.setNames(self.names())
        ag.setResnames(self.resnames())
        ag.setResnums(self.resnums())
        ag.setChids(self.chain_ids())
        return ag


class MDAnalisysAtomGroup(AtomGroup):
    def __init__(self, atom_group):
        if isinstance(atom_group, mda.Universe):
            self.atoms = atom_group.atoms
        else:
            self.atoms = atom_group

    def n_atoms(self):
        return self.n_atoms


    def positions(self):
        return self.atoms.positions


    def names(self):
        return self.atoms.names

    def resnames(self):
        return self.atoms.resnames

    def resnums(self):
        return self.atoms.resnums
    
    def chain_ids(self):
        return [a.chainID for a in self.atoms]

    def _c_alphas(self, chain=None):
        chain_sel = ''
        if chain is not None:
            chain_sel = f'chainid {chain} and '
        atoms = self.atoms.select_atoms(chain_sel + 'name CA')

        # select only atoms with large occupancy
        atoms_filtered = [ag[ag.occupancies.tolist().index(
            max(ag.occupancies))] for ag in atoms.groupby('resnums').values()]
        if len(atoms_filtered) < 2:
            return atoms
        else:
            return reduce(lambda a, b: a + b, atoms_filtered)

    def c_alphas_positions(self, chain=None):
        return self._c_alphas(chain).positions

    def c_alphas_residues(self, chain=None):
        return [{"resid": a.resid, "resname": a.resname}
                for a in self._c_alphas(chain)]

    def __repr__(self):
        return f'AtomGroup {self.n_atoms()} atoms'


class MDAnalisysTrajectory(Trajectory):
    def __init__(self, trajectory):
        self.trajectory = trajectory


    def __iter__(self):
        for ag in self.trajectory:
            yield MDAnalisysAtomGroup(ag)