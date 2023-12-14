""" This collection module allow a set of utilities to manage a collection of graphs.
"""
import pickle
import random

from typing import Callable, Optional
from torch_geometric.data import InMemoryDataset

from .graph import Graph
from .model import Target

class GraphCollection():
    """This ultility provides a way to organise and distribute multiple graphpro graphs

        Researches will normally work with a collection of protein graphs that
        requires organising, this utility class enable handle the colection.
    """

    def __init__(self, graphs: list[Graph], metadata: dict[str, str] = {}):
        self._graphs = graphs
        self._metadata = metadata

    def __len__(self):
        """Number of graphs in the collection"""
        return len(self._graphs)

    def __eq__(self, other):
        return self._graphs == other._graphs

    def __next__(self):
        for graph in self._graphs:
            yield graph

    def __iter__(self):
        return self.__next__()

    def __repr__(self):
        return f'GraphCollection size {len(self)}'

    def save(self, filename: str):
        """Serialise the collection into a file name, allowing a pipeline to be
        compose

        Args
            filename: location of the file to be stored.
        """
        with open(filename, "wb") as outfile:
            pickle.dump(self, outfile)

    def to_dataset(self, root: str,  node_encoders = [], target: Target = None) -> InMemoryDataset:
        """ Return the collection as InMemoryDataset
        """
        return GraphProDataset(root, self, node_encoders, target)


    def split(self,
              test_size: float = 0.8,
              seed: int = None):
        """ Split the graph collection into trainning and validation sets.
        """
        random.seed(seed)
        #avoid use random.choice rather decide the splits base on size
        test = []
        val = []
        for graph in self._graphs:
            if random.random() < test_size:
                test.append(graph)
            else:
                val.append(graph)
   
        return GraphCollection(test), GraphCollection(val)
   

    @staticmethod
    def load(filename: str):
        """ Loads a collection from a stored file, restoring the collection
            of graph to process
        """
        with open(filename, 'rb') as input:
            col = pickle.load(input)
            return col

class GraphProDataset(InMemoryDataset):
    def __init__(
        self,
        root: str,
        collection: GraphCollection,
        node_encoders = [],
        target: Target = None,
        transform: Optional[Callable] = None,
        pre_transform: Optional[Callable] = None,
        pre_filter: Optional[Callable] = None,
    ):
        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = self.collate([g.to_data(node_encoders, target) for g in collection])