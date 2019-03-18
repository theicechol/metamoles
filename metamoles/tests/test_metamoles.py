import metamoles
from metamoles import *
#from metamoles import create_cpd_info

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

#Tests for the RDKit molecular similarity functions
#Requires playground_df_cleaned_kegg_with_smiles.csv to be in the same directory for tests to pass.

def test_input_data():
    '''Tests input_data function in metamoles.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = metamoles.input_data(input_df)

    assert isinstance(test_df, pd.DataFrame) == True, """TypeError,
    function should return a pandas dataframe"""
    #assert

    return '1/1 tests successful'

def test_fingerprint_products():
    '''Tests fingerprint_products function in metamoles.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = metamoles.input_data(input_df)

    assert isinstance(metamoles.fingerprint_products(test_df), pd.DataFrame) == True, """TypeError,
    function should return a pandas dataframe"""
    #assert

    return '1/1 tests successful'

def test_sim_i_j():
    '''Tests sim_i_j function in metamoles.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = metamoles.fingerprint_products(metamoles.input_data(input_df))

    A = test_df.iloc[0]
    #B = test_df.iloc[1]
    #C = test_df.iloc[2]

    assert metamoles.sim_i_j(A, A) == 1, "Self correlation is broken"
    #assert metamoles.sim_i_j(A, B) == -1, "Standard correlation is broken"
    #assert metamoles.sim_i_j(A, C) == 0, "Standard correlation is broken"

    return '1/1 tests successful'

def test_sim_i_all():
    '''Test sim_i_all function in metamoles.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = metamoles.fingerprint_products(metamoles.input_data(input_df))
    metric = pd.DataFrame()

    assert metric.empty == True, """ShapeError, input metric dataframe
    should be initialized as empty"""

    for index, row in test_df.iterrows():
        assert metamoles.sim_i_all(test_df, index, row, metric) == None, """OutputError, function
        shouldn't return anything"""
        assert metric[index].all() >= 0 and metric[index].all() <= 1.0, """ValueError,
        metric should be between 0 and 1"""

    return "3/3 Tests successful"

def test_sim_metric():
    '''Test sim_i_all function in metamoles.py'''
    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = metamoles.fingerprint_products(metamoles.input_data(input_df))
    assert isinstance(metamoles.sim_metric(test_df), pd.DataFrame) == True, """TypeError,
    function should return a dataframe"""
    assert metamoles.sim_metric(test_df).isnull().values.any() == False, """ValueError,
    function-generated dataframe should not contain null values"""
    #assert test_df.size == metamoles.sim_metric(test_df).size, """ShapeError,
    #function-generated dataframe should be the same size as input dataframe"""

    return "2/2 Tests successful"

def test_calculate_dist():
	'''Test calculate_dist function in metamoles.py'''
	
	df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    
	test_df = metamoles.calculate_dist(df)
	assert isinstance(test_df, pd.DataFrame) == True, """TypeError,
	function should return a dataframe"""
    #assert len(test_df.columns) == 3+len(df.columns), """ShapeError, 
    #function should add 3 columns to dataframe"""


	return "1/1 Tests successful"

#Tests for the RDKit compound inform functions

def test_cpd_inform():
	'''Test cpd_inform function in metamoles.py'''	
	rapamycin = 'C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC'
	test = cpd_inform(rapamycin)
	
	assert test[0] == 51, "Carbon count is incorrect"
	assert test[1] == 79, "Hydrogen count is incorrect"
	assert type(test[-1]) == type(1.0), "TypeError: Molecular Weight should be float"
	
	return '3/3 Tests successful'

def test_create_cpd_info():
	'''Test create_cpd_info function in metamoles.py'''
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
	
	return '3/3 Tests successful'
	
def test_count_C():
	'''Test count_C function in metamoles.py'''
	mol='CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_C(mol) == 43, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_O():
	'''Test count_O function in metamoles.py'''
	mol='CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_O(mol) == 12, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_N():
	'''Test count_N function in metamoles.py'''
	mol='CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_N(mol) == 12, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_P():
	'''Test count_P function in metamoles.py'''
	mol='ClP(Cl)Cl'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_P(mol) == 1, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_S():
	'''Test count_S function in metamoles.py'''
	mol='CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_S(mol) == 2, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_X():
	'''Test count_X function in metamoles.py'''
	mol='ClP(Cl)Cl'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_X(mol) == 3, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
	
def test_count_H():
	'''Test count_H function in metamoles.py'''
	mol='CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O'
	mol=Chem.rdmolfiles.MolFromSmiles(mol)
	assert count_H(mol) == 66, "ValueError: Count is incorrect"
	return '1/1 Tests successful'
