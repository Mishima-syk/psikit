# -*- coding: utf-8 -*-
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import glob
import os
import uuid
import warnings
from tempfile import mkdtemp
from shutil import rmtree
from debtcollector import moves
from .util import mol2xyz
from .pymol_helper import run_pymol_server, save_pyscript
warnings.simplefilter('ignore')


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
        self.tempdir = mkdtemp()
        if self.debug:
            self.psi4.core.set_output_file("psikit_out.dat", True)
        else:
            self.psi4.core.be_quiet()

    def clean(self):
        rmtree(self.tempdir)
    
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

    def geometry(self, multiplicity=1):
        xyz = self.mol2xyz(multiplicity=multiplicity)
        self.psi4.geometry(xyz)

    def energy(self, basis_sets= "scf/6-31g**", return_wfn=True, multiplicity=1):
        self.geometry(multiplicity=multiplicity)
        scf_energy, wfn = self.psi4.energy(basis_sets, return_wfn=return_wfn)
        self.psi4.core.clean()
        self.wfn = wfn
        self.mol = self.xyz2mol()
        return scf_energy

    def optimize(self, basis_sets= "scf/6-31g**", return_wfn=True, name=None, multiplicity=1, maxiter=50):
        if not name:
            name = uuid.uuid4().hex
        self.psi4.core.IO.set_default_namespace(name)
        self.geometry(multiplicity=multiplicity)
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

    def mol2xyz(self, multiplicity=1):
        return mol2xyz(self.mol)

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

    def create_cube_files(self, gridspace=0.3):
        if self.wfn == None:
            print('wfn not found. run optimze()/energy()')
            return None
        else:
            a = self.wfn.nalpha()  # HOMO
            b = a + 1  # LUMO
            self.psi4.set_options({"cubeprop_tasks": ['ESP', 'FRONTIER_ORBITALS', 'Density', 'DUAL_DESCRIPTOR'],
                                   "cubic_grid_spacing": [gridspace, gridspace, gridspace],
                                   "cubeprop_filepath": self.tempdir
                                   })
            Chem.MolToMolFile(self.mol, os.path.join(self.tempdir, 'target.mol'))
            self.psi4.cubeprop(self.wfn)
    
    getMOview = moves.moved_function(create_cube_files, 'getMOview', __name__)

    def view_on_pymol(self, target='FRONTIER', maprange=0.05, gridspace=0.3):
        self.create_cube_files(gridspace=gridspace)
        run_pymol_server(self.tempdir, target=target, maprange=maprange)

    def save_frontier(self, gridspace=0.3, isotype="isosurface"):
        self.create_cube_files(gridspace=gridspace)
        save_pyscript(self.tempdir, isotype=isotype)  

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
        options = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
                   'VDW_POINT_DENSITY'  : 1.0,
                   'RESP_A'             : 0.0005,
                   'RESP_B'             : 0.1,
                   'RESTRAINT'          : True,
                   'RADIUS'             : {'Br':1.98, 'I':2.09,}
                   }
        charges = resp.resp([self.wfn.molecule()], options)
        #breakpoint()

        options['resp_a'] = 0.001
        resp.set_stage2_constraint(self.wfn.molecule(), charges[1], options)
        options['grid']=['%i_%s_grid.dat'%(1, self.wfn.molecule().name())]
        options['esp']=['%i_%s_grid_esp.dat'%(1, self.wfn.molecule().name())]

        charges2 = resp.resp([self.wfn.molecule()], options)

        atoms = self.mol.GetAtoms()
        for idx, atom in enumerate(atoms):
            atom.SetDoubleProp("EP", charges2[0][idx])
            atom.SetDoubleProp("RESP", charges2[1][idx])
        return charges2[1]


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
