from rdkit import Chem
from rdkit.Chem import rdFMCS


def get_different_atoms(sdf_ligand1, sdf_ligand2):
    # Get Mol object from SDF files
    r1 = Chem.SDMolSupplier(sdf_ligand1)
    for m1 in r1:
        mol1 = m1
    r2 = Chem.SDMolSupplier(sdf_ligand2)
    for m2 in r2:
        mol2 = m2
    # Find the minimal common substructure, returning indices of different atoms per ligand
    mcs = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    match1 = mol1.GetSubstructMatch(mcs_mol)
    target_atm1 = []
    for atom in mol1.GetAtoms():
        if atom.GetIdx() not in match1:
            target_atm1.append(atom.GetIdx())
    match2 = mol2.GetSubstructMatch(mcs_mol)
    target_atm2 = []
    for atom in mol2.GetAtoms():
        if atom.GetIdx() not in match2:
            target_atm2.append(atom.GetIdx())
    # if return_frag: THIS IS TO RETURN THE SMILES SYMBOLS
    #    return Chem.MolFragmentToSmiles(mol1, target_atm1), Chem.MolFragmentToSmiles(mol2, target_atm2) , mcs.smartsString

    return target_atm1, target_atm2, mcs.smartsString
