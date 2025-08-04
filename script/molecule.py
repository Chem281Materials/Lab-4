from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski

# Sample list of SMILES strings
SMILES_LIST = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
    "CC1=CC(=O)NC(C)=C1",            # Paracetamol
    "CCN(CC)CCCC(C)NC1=NC=NC2=C1C=CN2",  # Caffeine
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"   # Theobromine
]

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

def main():
    molecules = [Chem.MolFromSmiles(smiles) for smiles in SMILES_LIST]
    print(f"Total molecules: {len(molecules)}")

    lipinski_passed = [mol for mol in molecules if mol and passes_lipinski(mol)]
    print(f"Molecules passing Lipinski: {len(lipinski_passed)}")


if __name__ == "__main__":
    main()