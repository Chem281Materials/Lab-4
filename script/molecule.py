from rdkit import Chem

# Sample list of SMILES strings
SMILES_LIST = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
    "CC1=CC(=O)NC(C)=C1",            # Paracetamol
    "CCN(CC)CCCC(C)NC1=NC=NC2=C1C=CN2",  # Caffeine
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"   # Theobromine
]

def main():
    molecules = [Chem.MolFromSmiles(smiles) for smiles in SMILES_LIST]
    print(f"Total molecules: {len(molecules)}")


if __name__ == "__main__":
    main()