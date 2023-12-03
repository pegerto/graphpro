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
    '': 'X',
}


def one_letter_res(res: str) -> str:
    # TODO: Add warning here, include a guesser in the annotation
    return _RES_NAME.get(res, 'X')

def res_letters() -> list[str]:
    return list(_RES_NAME.values())