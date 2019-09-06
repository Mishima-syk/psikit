import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from psikit import Psikit
from psikit import sapt
import rdkit
import pytest


def test_read_from_smiles():
    pk = Psikit()
    pk.read_from_smiles("C")
    assert type(pk.mol) is rdkit.Chem.rdchem.Mol

def test_read_from_molfiles():
    pk = Psikit()
    pk.read_from_molfile("tests/test.mol")
    assert type(pk.mol) is rdkit.Chem.rdchem.Mol

def test_energy():
    pk = Psikit()
    pk.read_from_smiles("C")
    energy = pk.energy()
    assert pytest.approx(-40.19996313, 0.000000005) == energy

def test_energy_sto3g():
    pk = Psikit()
    pk.read_from_smiles("C")
    energy = pk.energy(basis_sets="scf/sto-3g")
    assert pytest.approx(-39.724747932, 0.000000005) == energy

def test_optimize():
    pk = Psikit()
    pk.read_from_smiles("C")
    energy = pk.optimize()
    assert pytest.approx(-40.20171733, 0.000000005) == energy
