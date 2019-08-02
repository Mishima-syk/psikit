import psikit
sapt = psikit.Sapt()
sapt.monomer1_from_molfile('phenol1.mol')
sapt.monomer2_from_molfile('phenol2.mol')
sapt.make_dimer()
p = sapt.dimer
sapt.run_fisapt()
