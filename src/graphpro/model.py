from functools import reduce

class AtomGroup():
    def n_atoms(self) -> int:
        pass

    def c_alphas_positions(self, chain: None):
        pass

    def c_alphas_residues(self, chain: None):
        pass
    

class MDAnalisysAtomGroup():
    def __init__(self, u):
        self.u = u

    def n_atoms(self):
        return self.u.atoms.n_atoms

    def _c_alphas(self, chain=None):
        chain_sel = ''
        if chain != None:
            chain_sel = f'chainid {chain} and '
        atoms = self.u.select_atoms(chain_sel + 'name CA')
        
        # select only residues with large occupancy
        atoms_filtered = [ag[ag.occupancies.tolist().index(max(ag.occupancies))] for ag in atoms.groupby('resnums').values()] 
        if len(atoms_filtered) == 0:
            return atoms
        elif len(atoms_filtered) == 1:
            return atoms
        else:
            return reduce(lambda a,b: a+b, atoms_filtered)

    def c_alphas_positions(self, chain = None):
         return self._c_alphas(chain).positions
    
    def c_alphas_residues(self, chain = None):
         return [{"resid": res.resid, "resname": res.resname } for res in self._c_alphas(chain).residues]
         
    def __repr__(self):
        return f'AtomGroup {self.n_atoms()} atoms'