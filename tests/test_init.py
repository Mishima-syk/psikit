from psikit import Psikit
import rdkit

class TestPsikit(object):
    def test_read_from_smiles(self):
        pk = Psikit()
        pk.read_from_smiles("C")
        assert type(pk.mol) is rdkit.Chem.rdchem.Mol