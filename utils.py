import torch
import torch.nn.functional as F
import torch_geometric
from torch_geometric.data import Data
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops
from rdkit import RDLogger
from typing import Any
from mapping_conf import ATOM, EDGE
from tqdm import tqdm

def to_smiles(data: 'torch_geometric.data.Data',
              kekulize: bool = True, data_name: str = 'MUTAG', add_H = False) -> Any:
    """Converts a :class:`torch_geometric.data.Data` instance to a SMILES
    string.

    Args:
        data (torch_geometric.data.Data): The molecular graph.
        kekulize (bool, optional): If set to :obj:`True`, converts aromatic
            bonds to single/double bonds. (default: :obj:`False`)
        data_name: The name of dataset
    """

    mol = Chem.RWMol()

    for i in range(data.num_nodes):
        # Some dataset does not have 
        if data_name in ["COX2", "BZR", "NCI1"]:
            atom = rdchem.Atom(torch.argmax(data.x[i]).item()+1)
        else:
            atom = rdchem.Atom(ATOM[data_name][torch.argmax(data.x[i]).item()])
        mol.AddAtom(atom)
    edges = [tuple(i) for i in data.edge_index.t().tolist()]
    visited = set()
    deleted = []
    
    for i in range(len(edges)):
        src, dst = edges[i]
        if tuple(sorted(edges[i])) in visited:
            continue
        if "edge_attr" in data.keys():
            bond_type = EDGE[data_name][torch.argmax(data.edge_attr[i]).item()]
            if bond_type == None:
                deleted.append(tuple(edges[i]))
            else:
                mol.AddBond(src, dst, bond_type)
        else:
            mol.AddBond(src, dst)

        visited.add(tuple(sorted(edges[i])))

    mol = mol.GetMol()

    mol = sanitize_mol(mol, add_H)
    if mol is None:
        return None

    # Chem.AssignStereochemistry(mol)

    return sanitize_smiles(get_smiles(mol), add_H)