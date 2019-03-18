# This cell is test_cpd_info.py

from rdkit import Chem
from rdkit.Chem import AllChem
from cpd_info import cpd_inform
from cpd_info import create_cpd_info

def test_cpd_inform():
    rapamycin = 'C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC'
    test = cpd_inform(rapamycin)
    
    assert test[0] == 51, "Carbon count is incorrect"
    assert test[1] == 79, "Hydrogen count is incorrect"
    assert type(test[-1]) == type(1.0), "TypeError: Molecular Weight should be float"
    
    return 'Test pass, yayyyyyy'
import pandas as pd

def test_create_cpd_info():
    
    """A unit test for create compound info"""
    
    df_master = pd.DataFrame(['C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O',
         'C([C@@H]1[C@@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O',
         'C([C@H]([C@H]([C@@H](C(=O)CO)O)O)O)O', 
         
'C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC'] 
, columns=['SMILES'])
    test = create_cpd_info(df_master)
    
    assert test['n_C'][0] == 6, "ValueError: Carbon count is incorrect"
    assert test['DoU'][3] == 13, "ValueError: Degree of Unsaturation in inaccurate"
    assert type(test['MW'][2]) == type(test['n_C'][0]), "TypeError: MW should be float"
    assert type(test['n_H'][3]) == type(test['n_C'][0]), "TypeError: All data should be float"
    
    return 'Test pass, you can use it to create compound info columns'
