_RES_NAME = {
    'GLY': 'G',
    'ALA': 'A',
    'LEU': 'L',
    'MET': 'M',
    'PHE': 'F',
    'TRP': 'W',
    'LYS': 'K',
    'GLN': 'Q',
    'GLU': 'E',
    'SER': 'S',
    'PRO': 'P',
    'VAL': 'V',
    'ILE': 'I',
    'CYS': 'C',
    'TYR': 'Y',
    'HSD': 'H',
    'HIS': 'H',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'THR': 'T',   
}

def one_letter_res(res: str) -> str:
    return _RES_NAME.get(res,'X')