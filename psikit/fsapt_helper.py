from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Recap
from rdkit.Chem import rdChemReactions
import os
    
Recap.reactions += tuple([rdChemReactions.ReactionFromSmarts('[c:1]-[X:2]>>[c:1]*.*[X:2]'), rdChemReactions.ReactionFromSmarts('[c:1]-[OH1:2]>>[c:1]*.*[OH1:2]')])

def get_neighbor_h(atom_idx, mol):
    atom = mol.GetAtomWithIdx(atom_idx)
    neis = atom.GetNeighbors()
    res = []
    for nei in neis:
        if nei.GetSymbol() == 'H':
            res.append(nei.GetIdx())
    return res

def make_feat_data(mol, offset=1):
    res = []
    check_atom = set()
    nohmol = Chem.RemoveHs(mol)
    recap_res = Recap.RecapDecompose(nohmol)
    leaves = [key.replace('*','').replace('()','') for key in recap_res.GetLeaves().keys()]
    leaves = [leave.replace('[H]', '') for leave in leaves if leave != '[H]']
    leaves = sorted(leaves, key=lambda x: Chem.MolFromSmarts(x).GetNumAtoms(), reverse=True)
    if len(leaves) == 0:

        line = [i for i in range(mol.GetNumAtoms())]
        line = [str(n + offset) for n in line]
        line = [Chem.MolToSmiles(mol)] + line
        return [line]
    for leavsmi in leaves:
        leav = Chem.MolFromSmarts(leavsmi)
        matches = mol.GetSubstructMatches(leav)
        for i, match in enumerate(matches):
            line = list(match)
            if len(check_atom & set(line)) > 0:
                continue
            check_atom = check_atom|set(line)
            for idx in match:
                nei = get_neighbor_h(idx, mol)
                line += nei
            line = [str(j + offset) for j in line]
            line = [leavsmi + '_' + str(i)] + line
            res.append(line)
    return res
