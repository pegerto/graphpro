from graphpro.model import AtomGroup

def compute_sasa(ag: AtomGroup):
    import freesasa
    
    # Avoid excessive logging during classifier guessing 
    freesasa.setVerbosity(freesasa.silent)

    structure = freesasa.Structure()  
    for a in ag.atoms:
        x,y,z = a.position
        structure.addAtom(a.type.rjust(2), a.resname, a.resnum.item(), a.chainID, x, y, z)
    return freesasa.calc(structure)