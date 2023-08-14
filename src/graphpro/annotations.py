from graphpro.graph import Graph
from graphpro.model import AtomGroup
from graphpro.util.residues import one_letter_res


class NodeAnnotation():
    def __init__(self):
        pass

    def generate(self, G: Graph, atom_group: AtomGroup): 
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


    def generate(self, G: Graph, atom_group: AtomGroup):
        for res  in atom_group.c_alphas_residues():
            node_id = G.get_node_by_resid(res['resid'])
            G.node_attr_add(node_id, self.attr_name, one_letter_res(res['resname']))
