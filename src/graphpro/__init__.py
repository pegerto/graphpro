from .model import MDAnalisysAtomGroup
from .graphgen import ProGraphGenerator

import MDAnalysis as mda

def md_analisys(u: mda.Universe, name: str = ""):
    return ProGraphGenerator(MDAnalisysAtomGroup(u), name)
