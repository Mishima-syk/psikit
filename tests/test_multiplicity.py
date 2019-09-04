import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from psikit import Psikit
from psikit import Sapt
import rdkit
import pytest


def test_mul1_energy():
    pk = Psikit()
    pk.read_from_smiles("CO")
    energy = pk.energy(multiplicity=1)
    assert pytest.approx(-115.0421871159, 0.000000005) == energy

def test_mul3_energy():
    pk = Psikit()
    pk.read_from_smiles("CO")
    energy = pk.energy(multiplicity=4)
    assert pytest.approx(-40.19996313, 0.000000005) == energy

