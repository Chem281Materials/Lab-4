from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import rdFingerprintGenerator

# Sample list of SMILES strings
SMILES_LIST = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",      # Aspirin
    "CC1=CC(=O)NC(C)=C1",            # Paracetamol
    "CCN(CC)CCCC(C)NC1=NC=NC2=C1C=CN2",  # Caffeine
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"   # Theobromine
]

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


if __name__ == "__main__":
    main()
