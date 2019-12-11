from glob import glob
import os

def run_pymol_server(tmpdir, target='FRONTIER', maprange=0.05):
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
    srv = xmlrpc.ServerProxy('http://localhost:9123')
    srv.do('delete *')
    srv.do('load '+os.path.join(tmpdir, 'target.mol'))
    srv.do('as sticks, target')
    if target == 'FRONTIER':
        homof = glob(os.path.join(tmpdir, 'Psi*_HOMO.cube'))[0]
        lumof = glob(os.path.join(tmpdir, 'Psi*_LUMO.cube'))[0]
        srv.do('load ' + homof + ',HOMO')
        srv.do('load ' + lumof + ',LUMO')
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

    outputpath = abb + 'mo.pse'
    srv.do(f'save {outputpath}')
    print('finished !')


def save_pyscript(tmpdir, isotype="isosurface"):
    homof = glob(os.path.join(tmpdir, 'Psi*_HOMO.cube'))[0]
    lumof = glob(os.path.join(tmpdir, 'Psi*_LUMO.cube'))[0]
    with open("frontier.py", "w") as f:
        f.write('from pymol import *\n')
        f.write('cmd.load("{0}", "HOMO")\n'.format(homof))
        f.write('cmd.load("{0}", "LUMO")\n'.format(lumof))
        f.write('cmd.load("{0}")\n'.format(os.path.join(tmpdir, "target.mol")))
        if isotype == "isomesh":        
            f.write('cmd.isomesh("HOMO_A", "HOMO", -0.02)\n')
            f.write('cmd.isomesh("HOMO_B", "HOMO", 0.02)\n')
            f.write('cmd.isomesh("LUMO_A", "LUMO", 0.02)\n')
            f.write('cmd.isomesh("LUMO_B", "LUMO", -0.02)\n')
        else:
            f.write('cmd.isosurface("HOMO_A", "HOMO", -0.02)\n')
            f.write('cmd.isosurface("HOMO_B", "HOMO", 0.02)\n')
            f.write('cmd.isosurface("LUMO_A", "LUMO", 0.02)\n')
            f.write('cmd.isosurface("LUMO_B", "LUMO", -0.02)\n')            
        f.write('cmd.color("blue", "HOMO_A")\n')       
        f.write('cmd.color("red", "HOMO_B")\n')       
        f.write('cmd.color("blue", "LUMO_A")\n')
        f.write('cmd.color("red", "LUMO_B")\n')
        f.write('cmd.disable("LUMO_A")\n')              
        f.write('cmd.disable("LUMO_B")\n')       