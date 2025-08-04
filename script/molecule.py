from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity

import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import pandas as pd
import numpy as np

# Substructure SMARTS for aromatic ring
AROMATIC_RING_SMARTS = "c1ccccc1"

GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def read_molecules_from_smi(filepath):
    mols = []
    with open(filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else smi
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol.SetProp("name", name)
                mols.append(mol)
    return mols


def passes_lipinski(mol):
    """Return True if molecule passes Lipinski's rule of 5."""
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    return (
        mw <= 500 and
        logp <= 5 and
        hbd <= 5 and
        hba <= 10
    )


def has_substructure(mol, smarts):
    """Return True if molecule contains the given substructure."""
    substructure = Chem.MolFromSmarts(smarts)
    return mol.HasSubstructMatch(substructure)


def get_fingerprint(mol):
    return GENERATOR.GetFingerprint(mol)


def compute_similarity_matrix(fps):
    size = len(fps)
    matrix = np.zeros((size, size))
    for i in range(size + 1):
        matrix[i, i] = 1.0
        for j in range(i + 1, size + 1):
            sim = TanimotoSimilarity(fps[i], fps[j])
            matrix[i, j] = sim
            matrix[j, i] = sim
    return matrix


def cluster_fingerprints(fps, labels):
    sim_matrix = compute_similarity_matrix(fps)
    distance_matrix = 1 - sim_matrix
    condensed = squareform(distance_matrix)
    linkage_matrix = linkage(condensed, method='average')
    
    # Plot dendrogram
    plt.figure(figsize=(8, 5))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=90)
    plt.title("Hierarchical Clustering of Molecules")
    plt.tight_layout()
    plt.savefig("molecular_dendrogram.png")
    plt.close()

    # Assign clusters (arbitrary threshold)
    clusters = fcluster(linkage_matrix, t=0.4, criterion='distance')
    return clusters


def main():
    molecules = read_molecules_from_smi("files/mols.smi")
    print(f"Total molecules: {len(molecules)}")

    lipinski_passed = [mol for mol in molecules if mol and passes_lipinski(mol)]
    print(f"Molecules passing Lipinski: {len(lipinski_passed)}")

    final_selection = [
        mol for mol in lipinski_passed
        if has_substructure(mol, AROMATIC_RING_SMARTS)
    ]
    print(f"Molecules with aromatic ring: {len(final_selection)}")

    if not final_selection:
        print("No molecules passed the filters.")
        return

    fingerprints = [get_fingerprint(mol) for mol in final_selection]
    smiles = [Chem.MolToSmiles(mol) for mol in final_selection]
    names = [mol.GetProp("id") for mol in final_selection]
    clusters = cluster_fingerprints(fingerprints, labels=names)

    # Save results
    df = pd.DataFrame({
        "SMILES": smiles,
        "Name": names,
        "Cluster": clusters
    })
    df.to_csv("clustered_molecules.csv", index=False)
    print("Clustering complete. Results saved to 'clustered_molecules.csv' and 'dendrogram.png'.")


if __name__ == "__main__":
    main()
