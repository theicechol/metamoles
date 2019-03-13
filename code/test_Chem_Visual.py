# this will be test_Chem_Visual.py

from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from Chem_Visual import Draw_mol
from Chem_Visual import Draw_Grid
from Chem_Visual import Draw_Rxn_1to1
from Chem_Visual import Draw_Rxn_2to1

rapamycin ='C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC'
smi_succinate = 'C(CC(=O)O)C(=O)O'
smi_butanediol = 'C(CCO)CO'
pyruvate = 'CC(=O)C(=O)[O-]'
benzaldehyde = 'C1=CC=C(C=C1)C=O'
PAC = 'CC(=O)[C@@H](C1=CC=CC=C1)O'

def test_Draw_mol():
    test = Draw_mol(rapamycin)
    return test

def test_Draw_Grid():
    test = Draw_Grid([smi_succinate, smi_butanediol])
    return test

def test_Draw_Rxn_1to1():
    test = Draw_Rxn_1to1(smi_succinate, smi_butanediol)
    return test

def test_Draw_Rxn_2to1():
    test = Draw_Rxn_2to1(pyruvate, benzaldehyde, PAC)
    return test
