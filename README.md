# Psikit: a thin wrapper library for Psi4 and RDKit

Inspired from this blog entry [Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

## Install RDKit and Psi4 from Conda

    conda install -c psi4 psi4
    conda install -c rdkit rdkit
    conda install -c psi4 resp # optional

## Install Psikit from Github

Psikit is under development.
We haven't uploaded psikit to PyPI yet, so plz install from github.

    pip install git+https://github.com/Mishima-syk/psikit


## Usage

### Single point calcuration

    from psikit import Psikit
    
    pk = Psikit()
    pk.read_from_smiles("c1ccccc1")
    print("SCF Energy: ", pk.energy())
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    x, y, z, total = pk.dipolemoment
    print("SCF Total Dipole Moment: {}".format(total))
    
    # SCF Energy:  -230.712279648862
    # HOMO:  -0.32848562009092513
    # LUMO:  0.1456515222506689
    # SCF Total Dipole Moment: 3.292464934070545e-05

### Structure optimization

    pk = Psikit()
    pk.read_from_smiles("c1ccccc1")
    print("Optimized SCF Energy: ", pk.optimize())
    # Optimizer: Optimization complete!
    # Optimized SCF Energy:  -230.71352354223438

### Calculate RESP Charge

    # REF http://ambermd.org/tutorials/advanced/tutorial1/files/resp_paper_jacs.pdf
    pk = Psikit()
    pk.read_from_smiles("CC(=O)O")
    pk.optimize()
    # Optimizer: Optimization complete!
    # -227.82180859253418
    pk.resp_charge
    # {'Electrostatic Potential Charges': array([-0.48689846,  0.89569528, -0.62426239, -0.67814945,  0.13380356,
    #    0.15316056,  0.1537424 ,  0.45290849]),
    # 'Restrained Electrostatic Potential Charges':
    #  array([-0.29156449,  0.81676972, -0.60648253, -0.66743813,  0.08243827, 0.10522699,  0.10565565,  0.45539453])}
    # pk.resp_charge method binds ESP and RESP charge to pk.mol object.

    for atom in pk.mol.GetAtoms(): 
        print(atom.GetSymbol(), "ESP:{}\tRESP:{}".format(atom.GetProp("EP"), atom.GetProp("RESP"))) 
    # C ESP:-0.48689845769479134    RESP:-0.29156449271229395
    # C ESP:0.8956952830247616  RESP:0.8167697222722353
    # O ESP:-0.6242623921340723 RESP:-0.6064825335690213
    # O ESP:-0.6781494526425089 RESP:-0.6674381342538175
    # H ESP:0.1338035637327267  RESP:0.08243826949989345
    # H ESP:0.15316056312500753 RESP:0.10522698717266236
    # H ESP:0.15374240384517748 RESP:0.1056556503974483
    # H ESP:0.4529084887436992  RESP:0.4553945311928933

### Jupyter notebook

- [RESP charge](examples/Rendering_RESP_charge/RESP%20charge%20of%20the%20tetrazole.ipynb)
- [Torsion Scan](examples/Torsion_scan/torsional_scan.ipynb)
- [Rendering Orbital](examples/Rendering_Orbital/Render_orbital.ipynb)