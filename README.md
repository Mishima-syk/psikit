from [Calculate HOMO and LUMO with Psi4](https://iwatobipen.wordpress.com/2018/08/24/calculate-homo-and-lumo-with-psi4-rdkit-psi4/)

# psikit as a wrapper library?

    import psikit
    import numpy as np
    from rdkit import Chem
    
    psikit.core.set_output_file("out.dat", True)
    psikit.set_memory('4 GB')
    psikit.set_num_threads(4)
    
    mol = Chem.MolFromSmiles("c1ccccc1")
    benz = psikit.geometry(mol)

# or psikit as a set of util functions library?

    import psi4
    from psikit import mol2xyz
    import numpy as np
    from rdkit import Chem
    psi4.core.set_output_file("out.dat", True)
    mol = Chem.MolFromSmiles("c1ccccc1")
    xyz, mol = mol2xyz(mol)
    psi4.set_memory('4 GB')
    psi4.set_num_threads(4)
    benz = psi4.geometry(xyz)