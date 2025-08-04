from rdkit import Chem
from script.molecule import passes_lipinski, has_substructure

aspirin = Chem.MolFromSmiles("CC(=O)OC1=CC=CC=C1C(=O)O")
methane = Chem.MolFromSmiles("C")
aromatic_smarts = "c1ccccc1"

def test_lipinski():
    assert passes_lipinski(aspirin)
    assert passes_lipinski(methane)  # Small molecule, should pass

def test_substructure():
    assert has_substructure(aspirin, aromatic_smarts)
    assert not has_substructure(methane, aromatic_smarts)
