import MDAnalysis as mda

from MDAnalysis.tests.datafiles import PDB, XTC

from graphpro import md_analisys
from graphpro.graphgen import KNN

u1 = mda.Universe(PDB, XTC)

def test_graph_generation_knn():
    G = md_analisys(u1).generate(KNN(k=3))
    assert(len(G.nodes()) == 214)
