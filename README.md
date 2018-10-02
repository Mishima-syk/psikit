from [Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

# Usage

    from psikit import Psikit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    pk = Psikit()
    mol = Chem.MolFromSmiles("c1ccccc1")
    pk.geometry(mol, addHs=True, optimize=True)
    print("SCF Energy: ", pk.energy())
    print("HOMO: ", pk.HOMO)
    print("LUMO: ", pk.LUMO)
    
    # SCF Energy:  -232.26253075623467
    # HOMO:  -0.2529008246007736
    # LUMO:  -0.006506665366600545
