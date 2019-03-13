import pandas as pd
import mol_sim

def test_input_data():
    '''Tests input_data function in mol_sim.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.input_data(input_df)

    assert isinstance(test_df, pd.DataFrame) == True, """TypeError,
    function should return a pandas dataframe"""
    #assert

    return '1/1 tests successful'

def test_fingerprint_products():
    '''Tests fingerprint_products function in mol_sim.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.input_data(input_df)

    assert isinstance(mol_sim.fingerprint_products(test_df), pd.DataFrame) == True, """TypeError,
    function should return a pandas dataframe"""
    #assert

    return '1/1 tests successful'

def test_split_by_enzyme():
    '''Tests split_by_enzyme function in mol_sim.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))

    assert isinstance(mol_sim.split_by_enzyme(test_df), list) == True, """TypeError,
    function should return a pandas dataframe"""
    #assert

    return '1/1 tests successful'

def test_sim_i_j():
    '''Tests sim_i_j function in mol_sim.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))

    A = test_df.iloc[0]
    #B = test_df.iloc[1]
    #C = test_df.iloc[2]

    assert mol_sim.sim_i_j(A, A) == 1, "Self correlation is broken"
    #assert mol_sim.sim_i_j(A, B) == -1, "Standard correlation is broken"
    #assert mol_sim.sim_i_j(A, C) == 0, "Standard correlation is broken"

    return '1/1 tests successful'

def test_sim_i_all():
    '''Test sim_i_all functionin mol_sim.py'''

    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))
    metric = pd.DataFrame()

    assert metric.empty == True, """ShapeError, input metric dataframe
    should be initialized as empty"""

    for index, row in test_df.iterrows():
        assert mol_sim.sim_i_all(test_df, index, row, metric) == None, """OutputError, function
        shouldn't return anything"""
        assert metric[index].all() >= 0 and metric[index].all() <= 1.0, """ValueError,
        metric should be between 0 and 1"""

    return "3/3 Tests successful"

def test_sim_metric():
    '''Test sim_i_all functionin mol_sim.py'''
    input_df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))
    assert isinstance(mol_sim.sim_metric(test_df), pd.DataFrame) == True, """TypeError,
    function should return a dataframe"""
    assert mol_sim.sim_metric(test_df).isnull().values.any() == False, """ValueError,
    function-generated dataframe should not contain null values"""
    #assert test_df.size == mol_sim.sim_metric(test_df).size, """ShapeError,
    #function-generated dataframe should be the same size as input dataframe"""

    return "2/2 Tests successful"

def test_calculate_dist():

    df = pd.read_csv('playground_df_cleaned_kegg_with_smiles.csv')
    test_df = mol_sim.calculate_dist(df)

    assert isinstance(test_df, pd.DataFrame) == True, """TypeError,
    function should return a dataframe"""
    #assert len(test_df.columns) == 3+len(df.columns), """ShapeError, 
    #function should add 3 columns to dataframe"""


    return "1/1 Tests successful"
