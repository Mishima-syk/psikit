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
    print("Optimized SCF Energy: ", pk.optimize())
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    x, y, z, total = pk.get_dipolemoment()
    print("SCF Total Dipole Moment: {}".format(total))

    # SCF Energy:  -230.71227964886216
    # HOMO:  -0.3284856200909327
    # LUMO:  0.14565152225066552
    # Optimizer: Optimization complete!
    # Optimized SCF Energy:  -230.7135235420281
    # HOMO:  -0.3306834797644457
    # LUMO:  0.14908632271767028
    # SCF Total Dipole Moment: 3.292465110164062e-05
