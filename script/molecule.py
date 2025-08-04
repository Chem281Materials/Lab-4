from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski

# Sample list of SMILES strings
SMILES_LIST = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",                 # Aspirin
    "CC1=CC(=O)NC(C)=C1",                       # Paracetamol
    "CCN(CC)CCCC(C)NC1=NC=NC2=C1C=CN2",         # Caffeine
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"              # Theobromine
    "CN1CCC(CC1)C2=CN=CC=C2",                   # Nicotine
    "CC(C)NCC(O)COC1=CC=CC=C1",                 # Propranolol
    "CN(C)C(=O)C1=CC=CC=C1Cl",                  # Lidocaine
    "CC(C)(C)OC(=O)N1CCCC1C(=O)NC2=CC=CC=C2",   # Atorvastatin fragment
    "CCOC(=O)C1=CC=CC=C1Cl",                    # Ethyl 4-chlorobenzoate
    "CNC1=CC=CC=C1Cl",                          # Chlorpheniramine
    "CC(C)COC(=O)C1=CC=CC=C1",                  # Ibuprofen
    "CC(C)NCC(O)COC1=CC=CC=C1Cl",               # Metoprolol
    "CC(C)N1CCC(CC1)NC2=NC=NC3=C2C=CN3",        # Theophylline derivative
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",             # Repeated (Theobromine)
    "CCOC(=O)C1=CC=CC=C1OC",                    # Methyl salicylate
    "COC1=CC=CC=C1OC",                          # Anisole
    "CC1=CC=C(C=C1)C(C)C(=O)O",                 # Naproxen
    "COC(=O)C1=CC=CC=C1Cl",                     # Methyl 4-chlorobenzoate
]

# Substructure SMARTS for aromatic ring
AROMATIC_RING_SMARTS = "c1ccccc1"

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
        hba <= 0
    )

def has_substructure(mol, smarts):
    """Return True if molecule contains the given substructure."""
    substructure = Chem.MolFromSmarts(smarts)
    return mol.HasSubstructMatch(substructure)

def main():
    molecules = [Chem.MolFromSmiles(smiles) for smiles in SMILES_LIST]
    print(f"Total molecules: {len(molecules)}")

    lipinski_passed = [mol for mol in molecules if mol and passes_lipinski(mol)]
    print(f"Molecules passing Lipinski: {len(lipinski_passed)}")

    final_selection = [
        mol for mol in lipinski_passed
        if has_substructure(mol, AROMATIC_RING_SMARTS)
    ]
    print(f"Molecules with aromatic ring: {len(final_selection)}")

if __name__ == "__main__":
    main()
