import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from psikit import Psikit
from psikit.sapt import Sapt
import rdkit
import pytest

def test_sapt():
    sapt = Sapt()
    sapt.psi4.core.clean()
    sapt.monomer1_from_molfile('tests/saptex/water1.mol')
    sapt.monomer2_from_molfile('tests/saptex/water2.mol')
    sapt.make_dimer()
    sapt0, Exch100, Elst10, Disp200, ExchDisp20, Ind20r, ExchInd20r = sapt.run_sapt()
    assert pytest.approx(-0.007335250653651525, 0.000000005) == sapt0

def test_fisapt():
    sapt = Sapt()
    sapt.psi4.core.clean()
    sapt.monomer1_from_molfile('tests/saptex/water1.mol')
    sapt.monomer2_from_molfile('tests/saptex/water2.mol')
    sapt.make_dimer()
    e = sapt.run_fisapt()
    assert pytest.approx(-0.00822262799944366, 0.000000005) == e
