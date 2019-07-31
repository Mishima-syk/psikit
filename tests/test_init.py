from psikit import Psikit
from psikit import Sapt
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
    sapt.monomer1_from_molfile('tests/saptex/phenol1.mol')
    sapt.monomer2_from_molfile('tests/saptex/phenol2.mol')
    sapt.make_dimer()
    e = sapt.run_fisapt()
    assert pytest.approx(-0.011385703498804293, 0.000000005) == e
