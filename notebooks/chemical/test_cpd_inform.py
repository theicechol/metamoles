from rdkit import Chem
from rdkit.Chem import AllChem
from cpd_inform import cpd_inform
def test_cpd_inform():
    rapamycin = 'C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC'
    test = cpd_inform(rapamycin)
    
    assert test['n_C'] == 51, "Carbon count is incorrect"
    assert test['n_H'] == 79, "Hydrogen count is incorrect"
    assert type(test['DoU']) == type(5), "TypeError: Degree of Saturation should be an int"
    assert type(test['MW']) == type(1.0), "TypeError: Molecular Weight should be float"
    
    return 'Test pass, yayyyyyy'
