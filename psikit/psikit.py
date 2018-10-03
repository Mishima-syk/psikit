# -*- coding: utf-8 -*-
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class Psikit(object):
    def __init__(self, threads=2, memory=4):
        import psi4
        self.psi4 = psi4
        self.psi4.core.set_output_file("psikit_out.dat", True)
        self.psi4.set_memory("{} GB".format(memory))
        self.psi4.set_num_threads(threads)
        self.wfn = None

    def rdkit_optimize(self, mol, addHs=True):
        if addHs:
            mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
        AllChem.UFFOptimizeMolecule(mol)
        return mol

    def geometry(self, mol):
        if type(mol) == rdkit.Chem.rdchem.Mol:
            xyz = mol2xyz(mol)
        self.psi4.geometry(xyz)        
            
    def energy(self, basis_sets= "scf/6-31g**", return_wfn=True):
        scf_energy, wfn = self.psi4.energy(basis_sets, return_wfn=return_wfn)
        self.wfn = wfn
        return scf_energy

    def optimize(self, basis_sets= "scf/6-31g**", return_wfn=True):
        scf_energy, wfn = self.psi4.optimize(basis_sets, return_wfn=return_wfn)
        self.wfn = wfn
        return scf_energy

    @property
    def dipolemoment(self, basis_sets="scf/6-31g**", return_wfn=True):
        #  The three components of the SCF dipole [Debye]
        x = self.psi4.get_variable('SCF DIPOLE X')
        y = self.psi4.get_variable('SCF DIPOLE Y')
        z = self.psi4.get_variable('SCF DIPOLE Z')
        total = np.sqrt(x * x + y * y + z * z)
        return (x, y, z, total)

    @property
    def HOMO(self):
        return self.wfn.epsilon_a_subset('AO', 'ALL').np[self.wfn.nalpha()-1]

    @property
    def LUMO(self):
        return self.wfn.epsilon_a_subset('AO', 'ALL').np[self.wfn.nalpha()]



def mol2xyz(mol):
    xyz_string = "\n"
    for _, atom in enumerate(mol.GetAtoms()):
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
    xyz_string += "units angstrom\n"
    return xyz_string
