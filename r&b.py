
import argparse
from collections import defaultdict
import csv
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import time

def load_smiles_dataset(file_path):
    """
    Loads SMILES strings from a CSV file.
    
    Args:
        file_path (str): Path to the CSV file containing SMILES strings.
        
    Returns:
        list: List of SMILES strings.
    """
    data = pd.read_csv(file_path)
    smiles_list = data['smiles'].tolist()  
    return smiles_list

def extract_motifs(mol):
    """
    Extracts ring and non-ring bond motifs from a molecule.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object.
        
    Returns:
        list: List of extracted motifs (rings and non-ring bonds) as SMILES strings.
    """
    motifs = []
    
    # Extract ring motifs
    rings = Chem.GetSymmSSSR(mol)
    for ring in rings:
        ring_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(ring))
        motifs.append(ring_smiles)
    
    # Extract non-ring bond motifs
    for bond in mol.GetBonds():
        if not bond.IsInRing():
            bond_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            motifs.append(bond_smiles)
    
    return motifs

def generate_motif_vocabulary(smiles_list, freq_threshold=10):
    """
    Generates the motif vocabulary from a list of SMILES strings by extracting ring and non-ring bond motifs.
    
    Args:
        smiles_list (list): List of SMILES strings.
        freq_threshold (int): Frequency threshold for motif selection (default: 10).
        
    Returns:
        dict: Dictionary representing the motif vocabulary, where keys are motifs and values are frequencies.
    """
    motif_counts = defaultdict(int)
    
    for smiles in tqdm(smiles_list, desc="Generating Motif Vocabulary"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        motifs = extract_motifs(mol)
        for motif in motifs:
            motif_counts[motif] += 1
    
    motif_vocab = {motif: count for motif, count in motif_counts.items() if count >= freq_threshold}
    
    return motif_vocab

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', default='dataset/tox21/raw/tox21.csv', help='Path to the SMILES data CSV file')
    # parser.add_argument('--output', required=True, help='Path to save the motif vocabulary JSON file')
    parser.add_argument('--freq_threshold', type=int, default=10, help='Frequency threshold for motif selection (default: 10)')
    args = parser.parse_args()
    
    smiles_list = load_smiles_dataset(args.data)

    start_time = time.time()
    motif_vocab = generate_motif_vocabulary(smiles_list, args.freq_threshold)

    end_time = time.time()
    total_time = end_time - start_time
    avg_time_per_mol = total_time / len(smiles_list)
    
    print(f'Generated {len(motif_vocab)} motifs in {total_time:.2f} seconds')
    print(f'Average time per molecule: {avg_time_per_mol:.4f} seconds')

    print(motif_vocab)

    