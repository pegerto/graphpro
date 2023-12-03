""" This collection module allow a set of utilities to manage a collection of graphs.
"""
import pickle

from torch_geometric.data import InMemoryDataset
from typing import Callable, Optional

from .graph import Graph

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

    def to_dataset(self, root: str,  node_encoders = []) -> InMemoryDataset:
        """ Return the collection as InMemoryDataset
        """
        return GraphProDataset(root, self, node_encoders)

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
        transform: Optional[Callable] = None,
        pre_transform: Optional[Callable] = None,
        pre_filter: Optional[Callable] = None,
    ):
        super().__init__(root, transform, pre_transform, pre_filter)
        self.data, self.slices = self.collate([g.to_data(node_encoders) for g in collection])
