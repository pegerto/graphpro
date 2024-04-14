from .model import MDAnalisysAtomGroup, MDAnalisysTrajectory
from .graphgen import ProGraphGenerator

import MDAnalysis as mda

def md_analisys(u: mda.Universe, name: str = ""):
    """ This helper connects MDAnalysis with Graphpro, returning a 
        generator that builds a graph representation from a static 3D graph 
        structure or a trajectory.

    Parameters
        ----------
        u : MDAnalysis univererse
        name: protein name for reference
            Default: empty string

    """
    return ProGraphGenerator(MDAnalisysAtomGroup(u.atoms), MDAnalisysTrajectory(u.trajectory), name)
