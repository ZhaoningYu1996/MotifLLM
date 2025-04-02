import torch
import torch.nn as nn
from torch.nn import Linear
import torch.nn.functional as F
from tree import Tree, TreeDataset
from torch_geometric.data import Data, Batch
import numpy as np
from collections import defaultdict, deque
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdmolops
from sklearn.model_selection import train_test_split
import os

def get_tree(self, data, test=False):
    # Get the tree from the dataset
    tree = Tree(data, self.data_name, self.add_H)
    tree.transform()
    # print(stop)
    node_motif_map = defaultdict(set)
    node_indices = []
    for i, motif in enumerate(tree.fragments):
        node_indices.append(tree.atom_list[i])
        if motif not in self.motif_id:
            if test:
                return None
            self.motif_id[motif] = len(self.motif_id)
            self.id_motif[self.motif_id[motif]] = motif
        for node in tree.atom_list[i]:
            node_motif_map[node].add(i)
    # Create x for the tree, each node is a motif, the node feature is the motif id
    x = torch.tensor([[self.motif_id[motif]] for motif in tree.fragments], dtype=torch.long)
    edge_index = []

    # Need to think about singlton nodes.
    for bond in tree.bond_list:
        set1 = node_motif_map[bond[0]]
        set2 = node_motif_map[bond[1]]

        for node1 in set1:
            for node2 in set2:
                edge_index.append([node1, node2])
                edge_index.append([node2, node1])
    if len(edge_index) == 0:
        edge_index = torch.tensor([[], []], dtype=torch.long)
    else:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    tree_data = Data(x=x, edge_index=edge_index, node_ori_map=node_indices, data=data)
    
    return tree_data