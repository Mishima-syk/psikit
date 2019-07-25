from psikit import Sapt
sapt=Sapt()
sapt.monomer1_from_molfile('water1.mol')
sapt.monomer2_from_molfile('water2.mol')
sapt.make_dimer()
print(sapt.dimer)
sapt0, Exch100, Elst10, Disp200, ExchDisp20, Ind20r, ExchInd20r = sapt.run_sapt()
print('total sapt0', sapt0)
print('Exch10 (S^2)', Exch100)
print('Elst10', Elst10)
print('Disp20', Disp200)
print('Exch-Disp20', ExchDisp20)
print('Ind20,r', Ind20r)
print('Exch-Ind20,r', ExchInd20r)
