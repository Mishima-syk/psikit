# Psikit: a thin wrapper library for Psi4 and RDKit

Inspired from this blog entry [Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

## Install RDKit and Psi4 from Conda

    conda install -c psi4 psi4
    conda install -c rdkit rdkit

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
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    x, y, z, total = pk.dipolemoment
    print("SCF Total Dipole Moment: {}".format(total))

    # Optimizer: Optimization complete!
    # Optimized SCF Energy:  -230.71352354223438
    # HOMO:  -0.3306834775917495
    # LUMO:  0.14908631857977886
    # SCF Total Dipole Moment: 2.527398024898661e-05