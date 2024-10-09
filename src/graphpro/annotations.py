import math

import torch
import torch.nn.functional as F
import numpy as np

from graphpro.graph import Graph
from graphpro.model import AtomGroup, NodeTarget
from graphpro.util.residues import one_letter_res, res_letters
from graphpro.util.sasa import compute_sasa
from graphpro.util.modes import compute_gnm_slow_modes, compute_anm_slow_modes
from graphpro.util.dssp import compute_dssp, DSSP_CLASS
from graphpro.util.polarity import POLARITY_CLASSES, residue_polarity

class NodeTargetBinaryAttribute(NodeTarget):
    """ Binary target, creates a binary one_hot encoding of the property
        been present on the node or not.
    """
    def encode(self, G: Graph) -> torch.Tensor:
        present = [self.attr_name in G.node_attr(n) for n in G.nodes()]
        return F.one_hot(torch.tensor(present, dtype=torch.int64), num_classes=2).to(torch.float)


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
        res_names = [G.node_attr(n)[self.attr_name] for n in G.nodes()]
        res_ids = [self.res_letters.index(name) for name in res_names]
        return F.one_hot(torch.tensor(res_ids, dtype=torch.int64), num_classes=len(self.res_letters)).to(torch.float)

class SASAResArea(NodeAnnotation):
    """ Generates a node an attribute with SASA area for individual 
        node residues.
    """
    def __init__(self, attr_name: str = 'sasa_area', chain_id: str = ''):
        """ Inits the residue generation
            Args:
                attr_name: name of the attribute
        """
        self.attr_name  = attr_name
        self.chain_id = chain_id
    
    def generate(self, G: Graph, atom_group: AtomGroup):
        sasa = compute_sasa(atom_group)
        chain_res = sasa.residueAreas()[self.chain_id]

        for k in chain_res:
            node_id = G.get_node_by_resid(int(k))
            # Protect agains nan in sasa computation
            total_area = 0 if math.isnan(chain_res[k].total) else chain_res[k].total
            G.node_attr_add(node_id, self.attr_name, total_area)
    
    def encode(self, G: Graph) -> torch.tensor:
        total_area = [G.node_attr(n)[self.attr_name] if self.attr_name in G.node_attr(n) else 0 for n in G.nodes()]
        return F.normalize(torch.tensor([total_area], dtype=torch.float).T, dim=(0,1))
    

class GNMSlowModes(NodeAnnotation):
    """ Calculate GNM and select specific slow modes per residue.s
    """
    def __init__(self, attr_name: str = 'gnm_slow_mode', modes: int = 3):
        """ Inits the residue generation
            Args:
                attr_name: name of the attribute prefix
                modes: number of modes.
        """
        self.attr_name  = attr_name
        self.modes = modes
    
    def generate(self, G: Graph, atom_group: AtomGroup):
        resids, modes_eigenvec = compute_gnm_slow_modes(atom_group, self.modes)
        for i, resid  in enumerate(resids):
             node_id = G.get_node_by_resid(int(resid))
             for m in range(0, self.modes):
                 G.node_attr_add(node_id, f"{self.attr_name}_{m}", modes_eigenvec[i,m].round(3))
    
    def encode(self, G: Graph) -> torch.tensor:
        values = [
            [G.node_attr(n)[f"{self.attr_name}_{m}"] for m in range(0, self.modes)] for n in G.nodes()
        ]
        return torch.tensor(values ,dtype=torch.float)

class ANMSlowModes(NodeAnnotation):
    """ Calculate ANM and select specific slow modes per residue
    """
    def __init__(self, attr_name: str = 'anm_slow_mode', modes: int = 3):
        """ Inits the residue generation
            Args:
                attr_name: name of the attribute prefix
                modes: number of modes.
        """
        self.attr_name  = attr_name
        self.modes = modes
    
    def generate(self, G: Graph, atom_group: AtomGroup):
        resids, modes_eigenvec = compute_anm_slow_modes(atom_group, self.modes)
        for i, resid  in enumerate(resids):
             node_id = G.get_node_by_resid(int(resid))
             for m in range(0, self.modes):
                 G.node_attr_add(node_id, f"{self.attr_name}_{m}", modes_eigenvec[i,m].round(3))
    
    def encode(self, G: Graph) -> torch.tensor:
        values = [
            [G.node_attr(n)[f"{self.attr_name}_{m}"] for m in range(0, self.modes)] for n in G.nodes()
        ]
        return torch.tensor(values ,dtype=torch.float)


class Polarity(NodeAnnotation):
    """ Anotates and encode per residue polarity
    """
    def __init__(self, attr_name: str = 'polarity'):
        """ Inits the residue generation
            Args:
                attr_name: name of the attribute prefix
        """
        self.attr_name  = attr_name

    def generate(self, G: Graph, atom_group: AtomGroup):
        for res in atom_group.c_alphas_residues():
            node_id = G.get_node_by_resid(res['resid'])
            G.node_attr_add(node_id, self.attr_name, residue_polarity(res['resname']))

    def encode(self, G: Graph) -> torch.tensor:
        polarity = [G.node_attr(n)[self.attr_name] for n in G.nodes()]
        polarity_class = [POLARITY_CLASSES.index(p) for p in polarity]
        return F.one_hot(torch.tensor(polarity_class, dtype=torch.int64), num_classes=len(POLARITY_CLASSES)).to(torch.float)
    

class NodeAttributeEncoder(NodeAnnotation):
    """ Encodes a node attribute that has been added manually to 
        the graph.
    """
    def __init__(self, attr_name: str):
        """ Attribute name
        """
        self.attr_name  = attr_name
    
    def generate(self, G: Graph, atom_group: AtomGroup):
        pass

    def encode(self, G: Graph) -> torch.tensor:
        values = [G.node_attr(n)[self.attr_name] if self.attr_name in G.node_attr(n) else 0 for n in G.nodes()]
        return F.normalize(torch.tensor([values], dtype=torch.float).T, dim=(0,1))
   
class DSSP(NodeAnnotation):
    """Encodes the specific secondary structure description per residues using DSSP"""
    def __init__(self, attr_name: str = 'dssp'):
        """ Attribute name
        """
        self.attr_name  = attr_name

    def generate(self, G: Graph, atom_group: AtomGroup):
        (resids, secondary) = compute_dssp(atom_group)
        for resid, sec in zip(resids, secondary):
            node_id = G.get_node_by_resid(resid)
            G.node_attr_add(node_id,self.attr_name, sec)

    def encode(self, G: Graph) -> torch.tensor:
        secondary = [G.node_attr(n)[self.attr_name] if self.attr_name in G.node_attr(n) else 'U' for n in G.nodes()]
        secondary_class = [DSSP_CLASS.index(p) for p in secondary]
        return F.one_hot(torch.tensor(secondary_class, dtype=torch.int64), num_classes=len(DSSP_CLASS)).to(torch.float)