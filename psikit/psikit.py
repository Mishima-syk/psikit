# -*- coding: utf-8 -*-
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import uuid

class Psikit(object):
    def __init__(self, threads=4, memory=4, debug=False):
        import psi4
        self.psi4 = psi4
        self.psi4.set_memory("{} GB".format(memory))
        #self.psi4.set_options({"save_jk": True})  # for JK calculation
        self.psi4.set_num_threads(threads)
        self.wfn = None
        self.mol = None
        self.debug = debug
        if self.debug:
            self.psi4.core.set_output_file("psikit_out.dat", True)
        else:
            self.psi4.core.be_quiet()

    def read_from_smiles(self, smiles_str, opt=True):
        self.mol = Chem.MolFromSmiles(smiles_str)
        if opt:
            self.rdkit_optimize()   

    def read_from_molfile(self, molfile, opt=True, removeHs=False):
        self.mol = Chem.MolFromMolFile(molfile, removeHs=removeHs)
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
        self.psi4.core.clean()
        self.wfn = wfn
        return scf_energy

    def optimize(self, basis_sets= "scf/6-31g**", return_wfn=True, name=None, maxiter=50):
        if not name:
            name = uuid.uuid4().hex
        self.psi4.core.IO.set_default_namespace(name)
        self.geometry()
        self.psi4.set_options({'GEOM_MAXITER':maxiter})
        try:
            scf_energy, wfn = self.psi4.optimize(basis_sets, return_wfn=return_wfn)
            self.wfn = wfn
        except self.psi4.OptimizationConvergenceError as cError:
            print('Convergence error caught: {0}'.format(cError))
            self.wfn = cError.wfn
            scf_energy = self.wfn.energy()
            self.psi4.core.clean()
        self.mol = self.xyz2mol()

        if not self.debug:
            self.psi4.core.opt_clean() # Seg fault will occured when the function is called before optimize.
        return scf_energy

    def set_options(self, **kwargs):
        """
        http://www.psicode.org/psi4manual/1.2/psiapi.html
        IV. Analysis of Intermolecular Interactions
        and 
        http://forum.psicode.org/t/how-can-i-change-max-iteration-in-energy-method/1238/2
        """
        self.psi4.set_options(kwargs)

    def mol2xyz(self):
        xyz_string = "\n{} 1\n".format(Chem.GetFormalCharge(self.mol))
        for atom in self.mol.GetAtoms():
            pos = self.mol.GetConformer().GetAtomPosition(atom.GetIdx())
            xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
        # the "no_com" stops Psi4 from moving your molecule to its center of mass, 
        # "no_reorient" stops it from spinning to align with axis of inertia
        xyz_string += "no_reorient\n"
        xyz_string += "no_com\n"
        xyz_string += "units angstrom\n"
        return xyz_string

    def xyz2mol(self, confId=0):
        natom = self.wfn.molecule().natom()
        mol_array_bohr = self.wfn.molecule().geometry().to_array()
        mol_array = mol_array_bohr * 0.52917721092
        nmol = Chem.Mol(self.mol)
        conf = nmol.GetConformer(confId)
        for i in range(natom):
            conf.SetAtomPosition(i, tuple(mol_array[i]))
        return nmol

    
    def clone_mol(self):
        return Chem.Mol(self.mol)

    def getMOview(self, gridspace=0.15):
        if self.wfn == None:
            print('please run optimze()/energy() at first!')
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
    
    def view_on_pymol(self):
        '''
        To use the function, user need to install pymol and run the pymol for server mode
        The command is pymol -R
        '''
        import sys
        import xmlrpc.client as xmlrpc
        filepath = os.getcwd()
        nalpha = self.wfn.nalpha()
        srv = xmlrpc.ServerProxy('http://localhost:9123')
        srv.do('delete *')
        srv.do('load '+os.path.join(filepath, 'target.mol'))
        srv.do('load '+os.path.join(filepath, f'Psi_a_{nalpha}_{nalpha}-A.cube'))
        srv.do('load '+os.path.join(filepath, f'Psi_b_{nalpha}_{nalpha}-A.cube'))
        srv.do('load '+os.path.join(filepath, f'Psi_a_{nalpha+1}_{nalpha+1}-A.cube'))
        srv.do('load '+os.path.join(filepath, f'Psi_b_{nalpha+1}_{nalpha+1}-A.cube'))
        srv.do(f'isomesh HOMO_A, Psi_a_{nalpha}_{nalpha}-A, -0.02')
        srv.do(f'isomesh HOMO_B, Psi_b_{nalpha}_{nalpha}-A, 0.02')
        srv.do(f'isomesh LUMO_A, Psi_a_{nalpha+1}_{nalpha+1}-A, -0.02')
        srv.do(f'isomesh LUMO_B, Psi_b_{nalpha+1}_{nalpha+1}-A, 0.02')
        srv.do('color blue, HOMO_A')
        srv.do('color red, HOMO_B')
        srv.do('color blue, LUMO_A')
        srv.do('color red, LUMO_B')
        outputpath = os.path.join(filepath, 'mo.pse')
        srv.do(f'save {outputpath}')
        print('finished !')

    def save_frontier(self, gridspace=0.15):
        if self.wfn == None:
            print('please run optimze()/energy() at first!')
        else:
            a = self.wfn.nalpha()  # HOMO
            b = a + 1  # LUMO
            self.psi4.set_options({"cubeprop_tasks":['orbitals'],
                                   "cubeprop_orbitals":[a, b, -a, -b],
                                   "cubic_grid_spacing":[gridspace, gridspace, gridspace]})
            Chem.MolToMolFile(self.mol, 'target.mol')
            self.psi4.cubeprop(self.wfn)

            homo_a = "Psi_a_{0}_{0}-A".format(a)
            homo_b = "Psi_b_{0}_{0}-A".format(a)
            lumo_a = "Psi_a_{0}_{0}-A".format(b)
            lumo_b = "Psi_b_{0}_{0}-A".format(b)
            with open("frontier.py", "w") as f:
                f.write('from pymol import *\n')
                f.write('cmd.load("{0}.cube")\n'.format(homo_a))
                f.write('cmd.load("{0}.cube")\n'.format(homo_b))
                f.write('cmd.load("{0}.cube")\n'.format(lumo_a))
                f.write('cmd.load("{0}.cube")\n'.format(lumo_b))
                f.write('cmd.load("target.mol")\n')               
                f.write('cmd.isomesh("HOMO_A", "{0}", -0.02)\n'.format(homo_a))
                f.write('cmd.isomesh("HOMO_B", "{0}", 0.02)\n'.format(homo_b))
                f.write('cmd.isomesh("LUMO_A", "{0}", 0.02)\n'.format(lumo_a))
                f.write('cmd.isomesh("LUMO_B", "{0}", -0.02)\n'.format(lumo_b))
                f.write('cmd.color("blue", "HOMO_A")\n')       
                f.write('cmd.color("red", "HOMO_B")\n')       
                f.write('cmd.color("blue", "LUMO_A")\n')
                f.write('cmd.color("red", "LUMO_B")\n')
                f.write('cmd.disable("LUMO_A")\n')              
                f.write('cmd.disable("LUMO_B")\n')         

    def save_fchk(self, filename="output.fchk"):
        fchk_writer = self.psi4.core.FCHKWriter(self.wfn)
        fchk_writer.write(filename)

    def save_cube(self):
        self.psi4.cubeprop(self.wfn)

    def calc_resp_charges(self):
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
            atom.SetDoubleProp("EP", charges[0][0][idx])
            atom.SetDoubleProp("RESP", charges[0][1][idx])
        return charges[0][1]

    def calc_mulliken_charges(self):
        '''
        Compute Mulliken Charges
        And return the results as numpy array.
        '''
        if self.wfn.molecule() == None:
            print('please run optimze() at first!')
            return None
        self.psi4.oeprop(self.wfn, 'MULLIKEN_CHARGES')
        mulliken_acp = self.wfn.atomic_point_charges()
        atoms = self.mol.GetAtoms()
        for idx, atom in enumerate(atoms):
            atom.SetDoubleProp("MULLIKEN", mulliken_acp.np[idx])
        return mulliken_acp.np

    def calc_lowdin_charges(self):
        '''
        Compute Lowdin Charges
        And return the results as numpy array.
        '''
        if self.wfn.molecule() == None:
            print('please run optimze() at first!')
            return None
        self.psi4.oeprop(self.wfn, 'LOWDIN_CHARGES')
        lowdin_acp = self.wfn.atomic_point_charges()
        atoms = self.mol.GetAtoms()
        for idx, atom in enumerate(atoms):
            atom.SetDoubleProp("LOWDIN", lowdin_acp.np[idx])
        return lowdin_acp.np


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


