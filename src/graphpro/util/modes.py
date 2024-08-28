from graphpro.model import AtomGroup
import prody as pdy

def compute_gnm_slow_modes(ag: AtomGroup, n_modes: int = 20):
    pdy_ag = ag.to_prody()
    calphas = pdy_ag.select('calpha')
    gnm = pdy.GNM()

    gnm.buildKirchhoff(calphas)
    gnm.calcModes(n_modes)

    return ([a.getResnum() for a in calphas], gnm.getEigvecs())


def compute_anm_slow_modes(ag: AtomGroup, n_modes: int = 20):
    pdy_ag = ag.to_prody()
    calphas = pdy_ag.select('calpha')
    anm = pdy.ANM()
    anm.buildHessian(calphas)
    anm.calcModes(n_modes)

    return ([a.getResnum() for a in calphas], anm.getEigvecs())