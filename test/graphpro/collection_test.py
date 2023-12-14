# pylint: disable=C

import numpy as np
import tempfile
import torch

from graphpro.graph import Graph
from graphpro.collection import GraphCollection

DISTANCES = np.array([[4, 2], [2, 4]])
POSITIONS_3D = np.array([[1, 2, 3], [-1, -2, -3]])
RESIDUES = [{"resid": 1, "resname": "GLY"}, {"resid": 2, "resname": "ALA"},]
SIMPLE_G = Graph("test", DISTANCES, POSITIONS_3D, RESIDUES)


def test_graph_collection():
    col = GraphCollection([SIMPLE_G])
    assert len(col) == 1


def test_load_and_save():
    col = GraphCollection([SIMPLE_G])
    f = tempfile.NamedTemporaryFile()
    col.save(f.name)
    col_res = GraphCollection.load(f.name)

    assert col == col_res


def test_compare_collections():
    col = GraphCollection([SIMPLE_G])
    assert col == col


def test_graphs_are_iterable():
    col = GraphCollection([SIMPLE_G, SIMPLE_G])
    for graph in col:
        assert graph is not None


def test_graph_collection_split():
    collection = GraphCollection([SIMPLE_G] * 100)
    train, test = collection.split(seed=42)    
    assert len(train) == 80
    assert len(test) == 20

def test_graphs_dataset():
    col = GraphCollection([SIMPLE_G, SIMPLE_G])
    ds = col.to_dataset('.')
    assert torch.all(ds[0].edge_index.eq(SIMPLE_G.to_data().edge_index))
    