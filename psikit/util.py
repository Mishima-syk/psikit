from rdkit import Chem

def mol2xyz(mol, multiplicity=1):
    charge = Chem.GetFormalCharge(mol)
    xyz_string = "\n{} {}\n".format(charge, multiplicity)
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
    return xyz_string