from glob import glob
import os

def run_pymol_server(nalpha=None, target='FRONTIER', maprange=0.05):
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
    srv = xmlrpc.ServerProxy('http://localhost:9123')
    srv.do('delete *')
    srv.do('load '+os.path.join(filepath, 'target.mol'))
    srv.do('as sticks, target')
    if target == 'FRONTIER':
        homof = glob(f'Psi_a_{nalpha}*_HOMO.cube')[0]
        lumof = glob(f'Psi_a_{nalpha+1}*_LUMO.cube')[0]
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


deg save_pyscript(nalpha=None):
    a = self.wfn.nalpha()  # HOMO
    b = a + 1  # LUMO
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