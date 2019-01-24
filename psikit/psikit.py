# -*- coding: utf-8 -*-
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class Psikit(object):
    def __init__(self, threads=4, memory=4):
        import psi4
        self.psi4 = psi4
        self.psi4.core.set_output_file("psikit_out.dat", True)
        self.psi4.set_memory("{} GB".format(memory))
        #self.psi4.set_options({"save_jk": True})  # for JK calculation
        self.psi4.set_num_threads(threads)
        self.wfn = None
        self.mol = None
        self.psi4.core.initialize()

    def read_from_smiles(self, smiles_str, opt=True):
        self.mol = Chem.MolFromSmiles(smiles_str)
        if opt:
            self.rdkit_optimize()   

    def rdkit_optimize(self, addHs=True):
        if addHs:
            self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
        AllChem.UFFOptimizeMolecule(self.mol)

    def geometry(self):
        xyz = self.mol2xyz()
        self.psi4.geometry(xyz)

    def energy(self, basis_sets= "scf/6-31g**", return_wfn=True):
        self.geometry()
        scf_energy, wfn = self.psi4.energy(basis_sets, return_wfn=return_wfn)
        self.wfn = wfn
        return scf_energy

    def optimize(self, basis_sets= "scf/6-31g**", return_wfn=True):
        self.geometry()
        scf_energy, wfn = self.psi4.optimize(basis_sets, return_wfn=return_wfn)
        self.wfn = wfn
        self.mol = self.xyz2mol()
        return scf_energy

    def mol2xyz(self):
        xyz_string = "\n"
        for _, atom in enumerate(self.mol.GetAtoms()):
            pos = self.mol.GetConformer().GetAtomPosition(atom.GetIdx())
            xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
        xyz_string += "units angstrom\n"
        return xyz_string

    def xyz2mol(self, confId=0):
        natom = self.wfn.molecule().natom()
        mol_array = self.wfn.molecule().geometry().to_array()
        nmol = Chem.Mol(self.mol)
        conf = nmol.GetConformer(confId)
        for i in range(natom):
            conf.SetAtomPosition(i, tuple(mol_array[i]))
        return nmol

    
    def clone_mol(self):
        return Chem.Mol(self.mol)

    def getMOview(self, gridspace=0.15):
        if self.wfn == None:
            print('please run optimze() at first!')
            return None
        else:
            a = self.wfn.nalpha()  # HOMO
            b = a + 1  # LUMO
            self.psi4.set_options({"cubeprop_tasks":['orbitals'],
                                   "cubeprop_orbitals":[a, b, -a, -b],
                                   "cubic_grid_spacing":[gridspace, gridspace, gridspace]})
            Chem.MolToMolFile(self.mol, 'target.mol')
            self.psi4.cubeprop(self.wfn)
            print('Done!')
            
    @property
    def resp_charge(self):
        if self.wfn.molecule() == None:
            print('please run optimze() at first!')
            return None
        try:
            import resp
        except:
            print('please install resp at first')
            print('conda install -c psi4 resp')
            return None
        # https://www.cgl.ucsf.edu/chimerax/docs/user/radii.html
        options = {'N_VDW_LAYERS'       : 4,
                   'VDW_SCALE_FACTOR'   : 1.4,
                   'VDW_INCREMENT'      : 0.2,
                   'VDW_POINT_DENSITY'  : 1.0,
                   'resp_a'             : 0.0005,
                   'RESP_B'             : 0.1,
                   'RADIUS'             : {'Br':1.98, 'I':2.09, }
                   }
        charges = resp.resp([self.wfn.molecule()], [options])
        atoms = self.mol.GetAtoms()
        for idx, atom in enumerate(atoms):
            atom.SetProp("EP_C", str(charges[0][0][idx]))
            atom.SetProp("RESP_C", str(charges[0][1][idx]))
        return { 'Electrostatic Potential Charges':charges[0][0], 
                 'Restrained Electrostatic Potential Charges':charges[0][1]}

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

    @property
    def coulomb_matrix(self):
        return self.wfn.jk().J[0].to_array()

    @property
    def exchange_matrix(self):
        return self.wfn.jk().K[0].to_array()


