from psikit import Psikit
import rdkit
import pytest


def test_read_from_smiles():
    pk = Psikit()
    pk.read_from_smiles("C")
    assert type(pk.mol) is rdkit.Chem.rdchem.Mol

def test_energy():
    pk = Psikit()
    pk.read_from_smiles("C")
    energy = pk.energy()
    assert pytest.approx(-40.1999631329161, 0.00000000000005) == energy

def test_energy_sto3g():
    pk = Psikit()
    pk.read_from_smiles("C")
    energy = pk.energy(basis_sets="scf/sto-3g")
    assert pytest.approx(-39.72474793261219, 0.00000000000005) == energy