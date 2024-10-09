import torch 
import numpy as np

from graphpro.model import AtomGroup

DSSP_ATOM_NUM  = {'N':0, 'CA': 1, 'C': 2, 'O': 3}

# simplification over DSSP implemenntation so not including G, I, S, T
DSSP_CLASS = ['-', 'H', 'E']

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
            resid_old = atom.resid
            resids.append(atom.resid)
            if atoms:
                coords.append(np.array(atoms))
                atoms = []     
        
        if iatom is not None:
            if iatom not in visited:
                visited.append(iatom)
                atoms.append(atom.position)
    if atoms and len(atoms) >= 4: # not adding partial residue atoms
        coords.append(np.array(atoms))
    
    npcoords = np.array(coords)
    secondary = pydssp.assign(torch.tensor(npcoords))
    return (resids, secondary)