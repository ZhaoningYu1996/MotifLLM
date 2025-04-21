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


bond_decoder_m = {1: Chem.rdchem.BondType.SINGLE, 2: Chem.rdchem.BondType.DOUBLE, 3: Chem.rdchem.BondType.TRIPLE}

def check_valency(mol):
    """
    Checks that no atoms in the mol have exceeded their possible valency

    Return:
        True if no valency issues, False otherwise
    """
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return True, None
    except ValueError as e:
        e = str(e)
        p = e.find('#')
        e_sub = e[p:]
        atomid_valence = list(map(int, re.findall(r'\d+', e_sub)))
        return False, atomid_valence
    
def correct_mol(mol):
    no_correct = False
    flag, _ = check_valency(mol)
    if flag:
        no_correct = True

    while True:
        flag, atomid_valence = check_valency(mol)
        if flag:
            break
        else:
            # Error message is one of the form: 
            # 'Explicit valence for atom # 0 O, 3, is greater than permitted
            # 'Explicit valence for atom # 15 Rn greater than permitted'
            # assert len(atomid_valence) == 2
            idx = atomid_valence[0]
            queue = []

            for b in mol.GetAtomWithIdx(idx).GetBonds():
                queue.append(
                    (b.GetIdx(), int(b.GetBondType()), b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                )
            queue.sort(key=lambda tup: tup[1], reverse=True)

            if len(queue) > 0:
                start = queue[0][2]
                end = queue[0][3]
                t = queue[0][1] - 1
                mol.RemoveBond(start, end)
                if t >= 1:
                    mol.AddBond(start, end, bond_decoder_m[t])

    return mol, no_correct


def sanitize_smiles(smiles, kekulize=True):
    mol = Chem.MolFromSmiles(smiles)
    mol = sanitize_mol(mol)
    if mol is None:
        return None
    if kekulize:
        Chem.Kekulize(mol, clearAromaticFlags=True)

    return Chem.MolToSmiles(mol)

def sanitize_mol(mol):
    try:
        smiles = Chem.MolToSmiles(mol)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
    except:
        return None
    return mol