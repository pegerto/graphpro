from .model import MDAnalisysAtomGroup, MDAnalisysTrajectory
from .graphgen import ProGraphGenerator

import MDAnalysis as mda


def md_analisys(u: mda.Universe, name: str = ""):
    return ProGraphGenerator(MDAnalisysAtomGroup(u.atoms), MDAnalisysTrajectory(u.trajectory), name)
