import torch
import torch.nn.functional as F

from graphpro.graph import Graph
from graphpro.model import AtomGroup
from graphpro.util.residues import one_letter_res, res_letters


class NodeAnnotation():
    def __init__(self):
        pass

    def generate(self, G: Graph, atom_group: AtomGroup):
        pass

    def encode(self) -> torch.tensor:
        """ Encode the node property to a tensor

            :returns a tensor the first dimension as number of nodes in the graph, 
            and dtype float (to be able to cat to other node properties.)
        """
        pass


class ResidueType(NodeAnnotation):
    """ Generates a one letter annotation of the residue type and add this to the node
        with the specific `attr_name`
    """

    def __init__(self, attr_name: str = 'resname'):
        """ Inits the residue generation
            Args:
                attr_name: name of the attribute
        """
        self.attr_name = attr_name
        self.res_letters = res_letters()

    def generate(self, G: Graph, atom_group: AtomGroup):
        for res in atom_group.c_alphas_residues():
            node_id = G.get_node_by_resid(res['resid'])
            G.node_attr_add(node_id, self.attr_name, one_letter_res(res['resname']))
    
    def encode(self, G: Graph) -> torch.tensor:
        res_names = [G.node_attr(n)['resname'] for n in G.nodes()]
        res_ids = [self.res_letters.index(name) for name in res_names]
        return F.one_hot(torch.tensor(res_ids, dtype=torch.int64), num_classes=len(self.res_letters)).to(torch.float)
