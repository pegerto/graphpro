import torch 
import numpy as np

from graphpro.model import AtomGroup

DSSP_ATOM_NUM  = {'N':0, 'CA': 1, 'C': 2, 'O': 3}

# simplification over DSSP implemenntation so not including G, I, S, T
# U is a custom made, where the calculation fail or there is not information
DSSP_CLASS = ['-', 'H', 'E', 'U']

def compute_dssp(atom_group: AtomGroup):
    import pydssp
    resids = []
    coords = []
    atoms = []
    visited = []
    resid_old = None

    # Asumes atom order 
    for atom in atom_group.atoms:
        iatom = DSSP_ATOM_NUM.get(atom.name, None)
        if resid_old != atom.resid:
            visited = []
            if len(atoms) > 0:
                resids.append(resid_old)
                coords.append(np.array(atoms))
                atoms = []
            resid_old = atom.resid

        if iatom is not None:
            if iatom not in visited:
                visited.append(iatom)
                atoms.append(atom.position)
    if atoms and len(atoms) >= 4: # not adding partial residue atoms
        coords.append(np.array(atoms))
        resids.append(resid_old)
    
    npcoords = np.array(coords)
    secondary = pydssp.assign(torch.tensor(npcoords))
    return (resids, secondary)