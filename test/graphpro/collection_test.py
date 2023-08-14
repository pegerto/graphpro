# pylint: disable=C

import numpy as np
import tempfile

from graphpro.graph import Graph
from graphpro.collection import GraphCollection

DISTANCES = np.array([[4,2],[2,4]])
POSITIONS_3D = np.array([[1,2,3], [-1,-2,-3]])
RESIDUES = [{"resid": 1, "resname": "GLY"}, {"resid": 2, "resname": "ALA"},]
SIMPLE_G = Graph("test", DISTANCES, POSITIONS_3D, RESIDUES)

def test_graph_collection():
    col = GraphCollection([SIMPLE_G])
    assert len(col) == 1


def test_load_and_save():
    col = GraphCollection([SIMPLE_G])
    f = tempfile.NamedTemporaryFile()
    print(f.name)
    col.save(f.name)
    col_res = GraphCollection.load(f.name)

    assert col == col_res


def test_compare_collections():
    col = GraphCollection([SIMPLE_G])
    assert col == col
