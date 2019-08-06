# Psikit: a thin wrapper library for Psi4 and RDKit

Inspired from the entry:[Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

## Install RDKit and Psi4 from Conda

    conda install -c psi4 psi4
    conda install -c rdkit rdkit
    conda install -c conda-forge debtcollector
    conda install -c psi4 resp # optional

## Install Psikit

Psikit is under development but you can install the current version of Psikit from pypi or conda.

**via conda**

    conda install -c iwatobipen psikit

**via pip**

    pip install psikit

**via pip from github**

    pip install git+https://github.com/Mishima-syk/psikit

## Testing Psikit

    pytest --disable-warnings -v

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
    pk.calc_resp_charges()
    # array([-0.32506898,  0.83672649, -0.61924915, -0.66135715,  0.10450057,
    #    0.10478188,  0.10780051,  0.45186584])

    for atom in pk.mol.GetAtoms(): 
        print(atom.GetSymbol(), "ESP:{}\tRESP:{}".format(atom.GetProp("EP"), atom.GetProp("RESP"))) 

    # C ESP:-0.49662019588648315	RESP:-0.3250689814483399
    # C ESP:0.91473263536048643		RESP:0.83672648554100837
    # O ESP:-0.63823808477114718	RESP:-0.61924915363703359
    # O ESP:-0.6763331997116846		RESP:-0.66135714989354499
    # H ESP:0.14625849864628995		RESP:0.10450056830656008
    # H ESP:0.14578513969681847		RESP:0.10478187811883517
    # H ESP:0.1530843954112609		RESP:0.1078005104750676
    # H ESP:0.45133081125445906		RESP:0.45186584253744722

    ### Compute Mulliken charges and Lowdin charges

    pk = Psikit()
    pk.read_from_smiles("CC(=O)O")
    pk.optimize() # or pk.energy()

    pk.calc_mulliken_charges()
    # array([-0.42203029,  0.72794785, -0.55419051, -0.59333358,  0.16369722,
    #    0.1636994 ,  0.15462075,  0.35958916])

    pk.calc_lowdin_charges()
    #array([-0.30006577,  0.33900448, -0.35983788, -0.28463832,  0.12439944,
    #    0.12810672,  0.11935266,  0.23367866])

### Rendering Molecular Orbitals
    
    # launch pymol as a RPC server, "pymol -R"
    from psikit import Psikit
    pk = Psikit()
    pk.read_from_smiles("c1ccccc1")
    pk.optimize(basis_sets="scf/sto-3g")
    pk.create_cube_files()
    pk.view_on_pymol()

![HOMO of benzene](images/homo.png)

### Adding RDKit mol object to Psikit object directly

    from psikit import Psikit
    pk = Psikit()
    pk.mol = your_mol_object

### Jupyter notebook

- [RESP charge](https://github.com/Mishima-syk/psikit/blob/master/examples/Rendering_RESP_charge/RESP%20charge%20of%20the%20tetrazole.ipynb)
- [Torsion Scan](https://github.com/Mishima-syk/psikit/blob/master/examples/Torsion_scan/torsional_scan.ipynb)
- [Rendering Orbital with VMD](https://github.com/Mishima-syk/psikit/blob/master/examples/Rendering_Orbital/Render_orbital.ipynb)
- [Rendering Orbital in PyMol](https://github.com/Mishima-syk/psikit/blob/master/examples/Rendering_Orbital_in_PyMol/Visualize_MO_in_PyMol.ipynb)
- [Charge Comparison](https://github.com/Mishima-syk/psikit/blob/master/examples/CHARGE_COMPARISON/charge_comparison.ipynb)
- [SAPT](https://github.com/Mishima-syk/psikit/blob/master/examples/example_sapt/sapt_ex.ipynb)
- [FSAPT](https://github.com/Mishima-syk/psikit/blob/master/examples/example_sapt/fsapt_ex.ipynb)

### License

Code released under the [BSD license](LICENSE).
