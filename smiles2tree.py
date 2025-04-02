from rdkit import Chem
import json
from enum import Enum

# bond_type_mapping = {
#     "SINGLE": 1,
#     "DOUBLE": 2,
#     "TRIPLE": 3,
# }

def smiles_to_tree(atom_format, bond_format, atom_type_enum, smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    # Chem.Kekulize(mol, clearAromaticFlags=True)
    atom_idx_to_atom = {}
    visited_atoms = set()
    queue = []
    added_bonds = set([])

    # Initialize with the first atom
    for atom in mol.GetAtoms():
        if atom.GetIdx() == 0:  # Start with the first atom
            atom_molecule = atom_format(atom_id=atom.GetIdx(), atom_type=atom_type_enum(atom.GetSymbol()), bonds=[])
            atom_idx_to_atom[atom.GetIdx()] = atom_molecule
            queue.append((None, atom.GetIdx(), None))  # (parent_idx, current_idx, bond_type)
            break

    # Process the queue using BFS
    while queue:
        parent_idx, current_idx, bond_type = queue.pop(0)
        current_atom = atom_idx_to_atom[current_idx]
        if parent_idx is not None:
            if current_idx in visited_atoms:
                current_atom = atom_format(atom_id=current_idx, atom_type=atom_type_enum(mol.GetAtomWithIdx(current_idx).GetSymbol()), bonds=[])
            parent_atom = atom_idx_to_atom[parent_idx]
            # print(bond_type)
            bond = bond_format(atom=current_atom, bond_type=bond_type)
            parent_atom.bonds.append(bond)
            visited_atoms.add(current_idx)
            visited_atoms.add(parent_idx)
            
        # Explore the neighbors
        for bond in mol.GetAtomWithIdx(current_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(current_idx)
            if (current_idx, neighbor_idx) not in added_bonds:
                neighbor_bond_type = (str(bond.GetBondType()).replace("BondType.", ""))
                # neighbor_bond_type = bond_type_mapping[str(bond.GetBondType()).replace("BondType.", "")]
                if neighbor_bond_type == "DATIVE":
                    neighbor_bond_type = "SINGLE"
                if neighbor_idx not in atom_idx_to_atom:
                    neighbor_atom = atom_format(atom_id=neighbor_idx, atom_type=atom_type_enum(mol.GetAtomWithIdx(neighbor_idx).GetSymbol()), bonds=[])
                    atom_idx_to_atom[neighbor_idx] = neighbor_atom
                queue.append((current_idx, neighbor_idx, neighbor_bond_type))
                added_bonds.add((current_idx, neighbor_idx))
        
    return atom_idx_to_atom[0]

# When serializing, convert enum to its value
def custom_json_serializer(obj):
    if isinstance(obj, Enum):
        return obj.value
    raise TypeError("Type not serializable")

# Example usage
if __name__ == "__main__":
    smiles = "C1CC1CO"  # Cyclopropane with a hydroxymethyl group
    molecule = smiles_to_tree(smiles)
    molecule_json = json.dumps(molecule.dict(), indent=0, default=custom_json_serializer)
    print("========================================")
    print(molecule_json)