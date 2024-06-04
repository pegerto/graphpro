"""
    Define per residue polarity
    p: polar
    a: apolar
    nc: negative charged
    pc: positive charged
    unk: unknown
"""
POLARITY = {
    'CYS': 'p', 
    'HIS': 'p', 
    'ASN': 'p', 
    'GLN': 'p', 
    'SER': 'p', 
    'THR': 'p', 
    'TYR': 'p', 
    'TRP': 'p',
    'ALA': 'p', 
    'PHE': 'a', 
    'GLY': 'a', 
    'ILE': 'a', 
    'VAL': 'a', 
    'MET': 'a', 
    'PRO': 'a', 
    'LEU': 'a',
    'GLU': 'nc', 
    'ASP': 'nc', 
    'LYS': 'nc', 
    'ARG': 'pc',
    'HSD': 'pc',
    '': 'unk' }

def residue_polarity(resname):
    return POLARITY.get(resname, 'unk')

POLARITY_CLASSES = list(set(POLARITY.values()))
