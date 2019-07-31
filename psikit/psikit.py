# -*- coding: utf-8 -*-
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import glob
import os
import uuid
import warnings
from debtcollector import moves
warnings.simplefilter('always')

def mol2xyz(mol):
    #xyz_string = "\n"
    xyz_string = "\n{} 1\n".format(Chem.GetFormalCharge(mol))
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
    # the "no_com" stops Psi4 from moving your molecule to its center of mass, 
    # "no_reorient" stops it from spinning to align with axis of inertia
    return xyz_string


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

    def create_cube_files(self, gridspace=0.3):
        if self.wfn == None:
            print('please run optimze()/energy() at first!')
            return None
        else:
            a = self.wfn.nalpha()  # HOMO
            b = a + 1  # LUMO
            self.psi4.set_options({"cubeprop_tasks":['ESP', 'FRONTIER_ORBITALS', 'Density', 'DUAL_DESCRIPTOR'],
                                   "cubic_grid_spacing":[gridspace, gridspace, gridspace]})
            Chem.MolToMolFile(self.mol, 'target.mol')
            self.psi4.cubeprop(self.wfn)
            print('Done!')
    
    getMOview = moves.moved_function(create_cube_files, 'getMOview', __name__)

    def view_on_pymol(self, target='FRONTIER', maprange=0.05):
        '''
        To use the function, user need to install pymol and run the pymol for server mode
        The command is pymol -R
        target ['ESP', 'FRONTIER', 'DUAL']
        '''
        targetlist = ['ESP', 'FRONTIER', 'DUAL']
        if target not in targetlist:
            raise Exception(f'Please set target from ESP, FRONTIER, DENSITY!!')
        import sys
        import xmlrpc.client as xmlrpc
        filepath = os.getcwd()
        nalpha = self.wfn.nalpha()
        srv = xmlrpc.ServerProxy('http://localhost:9123')
        srv.do('delete *')
        srv.do('load '+os.path.join(filepath, 'target.mol'))
        srv.do('as sticks, target')
        if target == 'FRONTIER':
            homof = glob.glob(f'Psi_a_{nalpha}*_HOMO.cube')[0]
            lumof = glob.glob(f'Psi_a_{nalpha+1}*_LUMO.cube')[0]
            srv.do('load '+os.path.join(filepath, homof) + ',HOMO')
            srv.do('load '+os.path.join(filepath, lumof)+ ',LUMO')
            srv.do(f'isosurface HOMO_A, HOMO, -0.02')
            srv.do(f'isosurface HOMO_B, HOMO, 0.02')
            srv.do(f'isosurface LUMO_A, LUMO, -0.02')
            srv.do(f'isosurface LUMO_B, LUMO, 0.02')
            srv.do('color blue, HOMO_A')
            srv.do('color red, HOMO_B')
            srv.do('color green, LUMO_A')
            srv.do('color yellow, LUMO_B')
            srv.do('set transparency, 0.2')
            srv.do('disable HOMO_A')
            srv.do('disable HOMO_B')
            abb = 'frontier_'
        elif target == 'ESP':
            srv.do('load '+ 'ESP.cube' + ', ESP')
            srv.do('show surface, target')
            srv.do(f'ramp_new cmap, ESP, [-{maprange}, {maprange}]')
            srv.do('color cmap, target')
            srv.do('set transparency, 0.2')
            abb = 'esp_'
        elif target == 'DUAL':
            dualf = glob.glob('DUAL*.cube')[0]
            srv.do('load '+ dualf + ', DUAL_DESC')
            srv.do('show surface, target')
            srv.do(f'ramp_new cmap, DUAL_DESC, [-{maprange}, {maprange}]')
            srv.do('color cmap, target')
            srv.do('set transparency, 0.2')
            abb = 'dual_'

        outputpath = os.path.join(filepath, abb + 'mo.pse')
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

class Sapt():
    def __init__(self, threads=4, memory=4, debug=False):
        import psi4
        from . import helper_SAPT

        self.psi4 = psi4
        self.helper_SAPT = helper_SAPT
        self.psi4.set_memory("{} GB".format(memory))
        #self.psi4.set_options({"save_jk": True})  # for JK calculation
        self.psi4.set_num_threads(threads)
        self.wfn = None
        self.monomer1 = None
        self.monomer2 = None
        self.dimer = None
        self.debug = debug
        if self.debug:
            self.psi4.core.set_output_file("psikit_out.dat", True)
        else:
            self.psi4.core.be_quiet()

    def monomer1_from_molfile(self, molfile, opt=True, removeHs=False):
        self.monomer1 = Chem.MolFromMolFile(molfile, removeHs=removeHs)

    def monomer2_from_molfile(self, molfile, opt=True, removeHs=False):
        self.monomer2 = Chem.MolFromMolFile(molfile, removeHs=removeHs)

    def make_dimer(self):
        xyz1 = mol2xyz(self.monomer1)
        xyz2 = mol2xyz(self.monomer2)
        self.dimer = "{}--\n{}".format(xyz1, xyz2)
        self.dimer += "no_reorient\n"
        self.dimer += "no_com\n"
        self.dimer += "units angstrom\n"
        
 
    def run_sapt(self, basis='aug-cc-pvdz', e_convergence=1e-8, d_convergence=1e-8, memory=4):
        self.psi4.set_options({'basis':basis,
                               'e_convergence':e_convergence,
                               'd_convergence':d_convergence})
        dimer = self.psi4.geometry(self.dimer)
        sapt = self.helper_SAPT.helper_SAPT(dimer, memory=memory)
        ### Start E100 Electostatics
        elst_timer = self.helper_SAPT.sapt_timer('electrostatics')
        Elst10 = 4 * np.einsum('abab', sapt.vt('abab'))
        elst_timer.stop()
        ### End E100 Electrostatics
        
        ### Start E100 Exchange
        exch_timer =self.helper_SAPT.sapt_timer('exchange')
        vt_abba = sapt.vt('abba')
        vt_abaa = sapt.vt('abaa')
        vt_abbb = sapt.vt('abbb')
        vt_abab = sapt.vt('abab')
        s_ab = sapt.s('ab')
        
        Exch100 = np.einsum('abba', vt_abba)
        
        tmp = 2 * vt_abaa - vt_abaa.swapaxes(2, 3)
        Exch100 += np.einsum('Ab,abaA', s_ab, tmp)
        
        tmp = 2 * vt_abbb - vt_abbb.swapaxes(2, 3)
        Exch100 += np.einsum('Ba,abBb', s_ab.T, tmp)
        
        Exch100 -= 2 * np.einsum('Ab,BA,abaB', s_ab, s_ab.T, vt_abab)
        Exch100 -= 2 * np.einsum('AB,Ba,abAb', s_ab, s_ab.T, vt_abab)
        Exch100 += np.einsum('Ab,Ba,abAB', s_ab, s_ab.T, vt_abab)
        
        Exch100 *= -2
        exch_timer.stop()
        ### End E100 (S^2) Exchange
        
        ### Start E200 Disp
        disp_timer = self.helper_SAPT.sapt_timer('dispersion')
        v_abrs = sapt.v('abrs')
        v_rsab = sapt.v('rsab')
        e_rsab = 1/(-sapt.eps('r', dim=4) - sapt.eps('s', dim=3) + sapt.eps('a', dim=2) + sapt.eps('b'))
        
        Disp200 = 4 * np.einsum('rsab,rsab,abrs->', e_rsab, v_rsab, v_abrs)
        ### End E200 Disp
        
        ### Start E200 Exchange-Dispersion
        
        # Build t_rsab
        t_rsab = np.einsum('rsab,rsab->rsab', v_rsab, e_rsab)
        
        # Build h_abrs
        vt_abar = sapt.vt('abar')
        vt_abra = sapt.vt('abra')
        vt_absb = sapt.vt('absb')
        vt_abbs = sapt.vt('abbs')
        
        tmp = 2 * vt_abar - vt_abra.swapaxes(2, 3)
        h_abrs = np.einsum('as,AbAr->abrs', sapt.s('as'), tmp)
        
        tmp = 2 * vt_abra - vt_abar.swapaxes(2, 3)
        h_abrs += np.einsum('As,abrA->abrs', sapt.s('as'), tmp)
        
        tmp = 2 * vt_absb - vt_abbs.swapaxes(2, 3)
        h_abrs += np.einsum('br,aBsB->abrs', sapt.s('br'), tmp)
        
        tmp = 2 * vt_abbs - vt_absb.swapaxes(2, 3)
        h_abrs += np.einsum('Br,abBs->abrs', sapt.s('br'), tmp)
        
        # Build q_abrs
        vt_abas = sapt.vt('abas')
        q_abrs =      np.einsum('br,AB,aBAs->abrs', sapt.s('br'), sapt.s('ab'), vt_abas)
        q_abrs -= 2 * np.einsum('Br,AB,abAs->abrs', sapt.s('br'), sapt.s('ab'), vt_abas)
        q_abrs -= 2 * np.einsum('br,aB,ABAs->abrs', sapt.s('br'), sapt.s('ab'), vt_abas)
        q_abrs += 4 * np.einsum('Br,aB,AbAs->abrs', sapt.s('br'), sapt.s('ab'), vt_abas)
        
        vt_abrb = sapt.vt('abrb')
        q_abrs -= 2 * np.einsum('as,bA,ABrB->abrs', sapt.s('as'), sapt.s('ba'), vt_abrb)
        q_abrs += 4 * np.einsum('As,bA,aBrB->abrs', sapt.s('as'), sapt.s('ba'), vt_abrb)
        q_abrs +=     np.einsum('as,BA,AbrB->abrs', sapt.s('as'), sapt.s('ba'), vt_abrb)
        q_abrs -= 2 * np.einsum('As,BA,abrB->abrs', sapt.s('as'), sapt.s('ba'), vt_abrb)
        
        vt_abab = sapt.vt('abab')
        q_abrs +=     np.einsum('Br,As,abAB->abrs', sapt.s('br'), sapt.s('as'), vt_abab)
        q_abrs -= 2 * np.einsum('br,As,aBAB->abrs', sapt.s('br'), sapt.s('as'), vt_abab)
        q_abrs -= 2 * np.einsum('Br,as,AbAB->abrs', sapt.s('br'), sapt.s('as'), vt_abab)
        
        vt_abrs = sapt.vt('abrs')
        q_abrs +=     np.einsum('bA,aB,ABrs->abrs', sapt.s('ba'), sapt.s('ab'), vt_abrs)
        q_abrs -= 2 * np.einsum('bA,AB,aBrs->abrs', sapt.s('ba'), sapt.s('ab'), vt_abrs)
        q_abrs -= 2 * np.einsum('BA,aB,Abrs->abrs', sapt.s('ba'), sapt.s('ab'), vt_abrs)
        
        # Sum it all together
        xd_absr = sapt.vt('absr')
        xd_absr += h_abrs.swapaxes(2, 3)
        xd_absr += q_abrs.swapaxes(2, 3)
        ExchDisp20 = -2 * np.einsum('absr,rsab->', xd_absr, t_rsab)
        
        disp_timer.stop()
        ### End E200 Exchange-Dispersion
        
        
        ### Start E200 Induction and Exchange-Induction
        
        # E200Induction and CPHF orbitals
        ind_timer = self.helper_SAPT.sapt_timer('induction')
        
        CPHF_ra, Ind20_ba = sapt.chf('B', ind=True)
        self.helper_SAPT.sapt_printer('Ind20,r (A<-B)', Ind20_ba)
        
        CPHF_sb, Ind20_ab = sapt.chf('A', ind=True)
        self.helper_SAPT.sapt_printer('Ind20,r (A->B)', Ind20_ab)
        
        Ind20r = Ind20_ba + Ind20_ab
        
        # Exchange-Induction
        
        # A <- B
        vt_abra = sapt.vt('abra')
        vt_abar = sapt.vt('abar')
        ExchInd20_ab  =     np.einsum('ra,abbr', CPHF_ra, sapt.vt('abbr'))
        ExchInd20_ab += 2 * np.einsum('rA,Ab,abar', CPHF_ra, sapt.s('ab'), vt_abar)
        ExchInd20_ab += 2 * np.einsum('ra,Ab,abrA', CPHF_ra, sapt.s('ab'), vt_abra)
        ExchInd20_ab -=     np.einsum('rA,Ab,abra', CPHF_ra, sapt.s('ab'), vt_abra)
        
        vt_abbb = sapt.vt('abbb')
        vt_abab = sapt.vt('abab')
        ExchInd20_ab -=     np.einsum('ra,Ab,abAr', CPHF_ra, sapt.s('ab'), vt_abar)
        ExchInd20_ab += 2 * np.einsum('ra,Br,abBb', CPHF_ra, sapt.s('br'), vt_abbb)
        ExchInd20_ab -=     np.einsum('ra,Br,abbB', CPHF_ra, sapt.s('br'), vt_abbb)
        ExchInd20_ab -= 2 * np.einsum('rA,Ab,Br,abaB', CPHF_ra, sapt.s('ab'), sapt.s('br'), vt_abab)
        
        vt_abrb = sapt.vt('abrb')
        ExchInd20_ab -= 2 * np.einsum('ra,Ab,BA,abrB', CPHF_ra, sapt.s('ab'), sapt.s('ba'), vt_abrb)
        ExchInd20_ab -= 2 * np.einsum('ra,AB,Br,abAb', CPHF_ra, sapt.s('ab'), sapt.s('br'), vt_abab)
        ExchInd20_ab -= 2 * np.einsum('rA,AB,Ba,abrb', CPHF_ra, sapt.s('ab'), sapt.s('ba'), vt_abrb)
        
        ExchInd20_ab +=     np.einsum('ra,Ab,Br,abAB', CPHF_ra, sapt.s('ab'), sapt.s('br'), vt_abab)
        ExchInd20_ab +=     np.einsum('rA,Ab,Ba,abrB', CPHF_ra, sapt.s('ab'), sapt.s('ba'), vt_abrb)
        
        ExchInd20_ab *= -2
        self.helper_SAPT.sapt_printer('Exch-Ind20,r (A<-B)', ExchInd20_ab)
        
        # B <- A
        vt_abbs = sapt.vt('abbs')
        vt_absb = sapt.vt('absb')
        ExchInd20_ba  =     np.einsum('sb,absa', CPHF_sb, sapt.vt('absa'))
        ExchInd20_ba += 2 * np.einsum('sB,Ba,absb', CPHF_sb, sapt.s('ba'), vt_absb)
        ExchInd20_ba += 2 * np.einsum('sb,Ba,abBs', CPHF_sb, sapt.s('ba'), vt_abbs)
        ExchInd20_ba -=     np.einsum('sB,Ba,abbs', CPHF_sb, sapt.s('ba'), vt_abbs)
        
        vt_abaa = sapt.vt('abaa')
        vt_abab = sapt.vt('abab')
        ExchInd20_ba -=     np.einsum('sb,Ba,absB', CPHF_sb, sapt.s('ba'), vt_absb)
        ExchInd20_ba += 2 * np.einsum('sb,As,abaA', CPHF_sb, sapt.s('as'), vt_abaa)
        ExchInd20_ba -=     np.einsum('sb,As,abAa', CPHF_sb, sapt.s('as'), vt_abaa)
        ExchInd20_ba -= 2 * np.einsum('sB,Ba,As,abAb', CPHF_sb, sapt.s('ba'), sapt.s('as'), vt_abab)
        
        vt_abas = sapt.vt('abas')
        ExchInd20_ba -= 2 * np.einsum('sb,Ba,AB,abAs', CPHF_sb, sapt.s('ba'), sapt.s('ab'), vt_abas)
        ExchInd20_ba -= 2 * np.einsum('sb,BA,As,abaB', CPHF_sb, sapt.s('ba'), sapt.s('as'), vt_abab)
        ExchInd20_ba -= 2 * np.einsum('sB,BA,Ab,abas', CPHF_sb, sapt.s('ba'), sapt.s('ab'), vt_abas)
        
        ExchInd20_ba +=     np.einsum('sb,Ba,As,abAB', CPHF_sb, sapt.s('ba'), sapt.s('as'), vt_abab)
        ExchInd20_ba +=     np.einsum('sB,Ba,Ab,abAs', CPHF_sb, sapt.s('ba'), sapt.s('ab'), vt_abas)
        
        ExchInd20_ba *= -2
        self.helper_SAPT.sapt_printer('Exch-Ind20,r (A->B)', ExchInd20_ba)
        ExchInd20r = ExchInd20_ba + ExchInd20_ab
        
        ind_timer.stop()
        ### End E200 Induction and Exchange-Induction
        
        print('\nSAPT0 Results')
        print('-' * 70)
        self.helper_SAPT.sapt_printer('Exch10 (S^2)', Exch100)
        self.helper_SAPT.sapt_printer('Elst10', Elst10)
        self.helper_SAPT.sapt_printer('Disp20', Disp200)
        self.helper_SAPT.sapt_printer('Exch-Disp20', ExchDisp20)
        self.helper_SAPT.sapt_printer('Ind20,r', Ind20r)
        self.helper_SAPT.sapt_printer('Exch-Ind20,r', ExchInd20r)
        
        print('-' * 70)
        sapt0 = Exch100 + Elst10 + Disp200 + ExchDisp20 + Ind20r + ExchInd20r
        self.helper_SAPT.sapt_printer('Total SAPT0', sapt0)        
        return sapt0, Exch100, Elst10, Disp200, ExchDisp20, Ind20r, ExchInd20r

    def run_fisapt(self, basis='jun-cc-pvdz', scf_type='df', d_convergence=1e-8, memory=4, fisapt_path='fsapt/', return_wfn=False):
        import shutil
        from distutils.dir_util import copy_tree
        from . import fsapt_helper
        self.psi4.set_options({'basis':basis,
                               'scf_type':scf_type,
                               'd_convergence':d_convergence,
                               'FISAPT_FSAPT_FILEPATH':fisapt_path,
                               'FISAPT_DO_PLOT': 'true'
                               })
        self.psi4.geometry(self.dimer)
        res = self.psi4.energy('fisapt0', return_wfn=return_wfn)
        copy_tree(self.psi4.core.get_datadir()+'/fsapt', fisapt_path)
        #sapt = self.helper_SAPT.helper_SAPT(dimer, memory=memory)
        feats1 = fsapt_helper.make_feat_data(self.monomer1, 1)
        feats2 = fsapt_helper.make_feat_data(self.monomer2, 1 + self.monomer1.GetNumAtoms())
        with open(os.path.join(fisapt_path, 'fA.dat'), 'w') as fA:
            for feat1 in feats1:
                fA.write(' '.join(feat1) + '\n')

        with open(os.path.join(fisapt_path,'fB.dat'), 'w') as fB:
            for feat2 in feats2:
                fB.write(' '.join(feat2) + '\n')

        print('done')
        return res
