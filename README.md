from [Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

# Usage

    from psikit import Psikit
    from rdkit import Chem
    from rdkit.Chem import AllChem
        
    pk = Psikit()
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = pk.rdkit_optimize(mol)
    pk.geometry(mol)
    print("SCF Energy: ", pk.energy())
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    x, y, z, total = pk.dipolemoment
    print("SCF Total Dipole Moment: {}".format(total))
    print("Optimized SCF Energy: ", pk.optimize())
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    x, y, z, total = pk.dipolemoment
    print("SCF Total Dipole Moment: {}".format(total))

    # SCF Energy:  -230.71227964886188
    # HOMO:  -0.3284856200909428
    # LUMO:  0.14565152225064903
    # SCF Total Dipole Moment: 3.292464935735843e-05
    # Optimizer: Optimization complete!
    # Optimized SCF Energy:  -230.71352354208835
    # HOMO:  -0.33068347419415384
    # LUMO:  0.1490863128964776
    # SCF Total Dipole Moment: 2.5269468477613993e-05