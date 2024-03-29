import torch

from functools import reduce

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


class MDAnalisysAtomGroup(AtomGroup):
    def __init__(self, atom_goup):
        self.atoms = atom_goup

    def n_atoms(self):
        return self.n_atoms

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