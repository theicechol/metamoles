import metamoles
from metamoles import *
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from pandas.util.testing import assert_frame_equal

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

def test_create_kegg_df():
    """Unit tests for create_kegg_df()
    """
    df1 = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    df2 = metamoles.create_kegg_df("test_kegg_compound_records.txt.gz","compound")
    # for enzyme database
    assert df1.shape == (3, 16)
    assert df1.columns.tolist() == ['classname', 'cofactor', 'comment', 'dblinks',
        'disease', 'effector', 'entry', 'genes', 'inhibitor', 'name',
        'pathway', 'product', 'reaction', 'structures', 'substrate', 'sysname']
    assert df1['entry'].tolist() == ['1.1.1.1', '1.1.1.2', '1.1.1.3']
    # for compound database
    assert df2.shape == (55, 8)
    assert df2.columns.tolist() == ['dblinks', 'entry', 'enzyme', 'formula', 'mass',
            'name', 'pathway', 'structures']
    assert df2['entry'].tolist() == ['C00001', 'C00002', 'C00003', 'C00004',
            'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00011',
            'C00012', 'C00013', 'C00014', 'C00015', 'C00016', 'C00017', 'C00018',
            'C00019', 'C00020', 'C00021', 'C00022', 'C00023', 'C00024', 'C00025',
            'C00026', 'C00027', 'C00028', 'C00029', 'C00030', 'C00031', 'C00032',
            'C00033', 'C00034', 'C00035', 'C00036', 'C00037', 'C00038', 'C00039',
            'C00040', 'C00041', 'C00042', 'C00043', 'C00044', 'C00045', 'C00046',
            'C00047', 'C00048', 'C00049', 'C00050', 'C00051', 'C00052', 'C00053',
            'C00054', 'C00055']
    return

def test_select_promiscuous_enzymes():
    """Unit tests for select_promiscuous_enzymes()
    """
    df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    test_prom_df = metamoles.select_promiscuous_enzymes(df)
    expected_column = ['entry', 'reaction', 'product', 'substrate']
    actual_column = test_prom_df.columns.tolist()
    assert test_prom_df.shape == (1,4)
    assert expected_column == actual_column, "expected column names and actual column names do not match"

    return

def test_parse_compound_ids():
    """Unit tests for parse_compound_ids()
    """
    df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    expected = ['C00071','C00004','C00080','C01450','C00071','C00005',
    'C00080','C00441','C00004','C00005','C00080']
    actual = metamoles.parse_compound_ids(df["product"])
    assert expected == actual, "expected result of parse_compound_ids does not match the actual result"

    return

def test_parse_pubchem_ids():
    """Unit tests for parse_pubchem_ids()
    """
    compound_df = metamoles.create_kegg_df("test_kegg_compound_records.txt.gz","compound")
    PubChemID_list = []
    for _, row in compound_df.iterrows():
        pubchem_id = metamoles.parse_pubchem_ids(row['dblinks'])
        PubChemID_list.append(pubchem_id)
    assert PubChemID_list == ['3303', '3304', '3305', '3306', '3307', '3308',
        '3309', '3310', '3311', '3312', '3313', '3314', '3315', '3316', '3317',
        '3318', '3319', '3320', '3321', '3322', '3323', '3324', '3325', '3326',
        '3327', '3328', '3329', '3330', '3331', '3332', '3333', '3334', '3335',
        '3336', '3337', '3338', '3339', '3340', '3341', '3342', '3343', '3344',
        '3345', '3346', '3347', '3348', '3349', '3350', '3351', '3352', '3353',
        '3354', '3355', '3356', '3357']
    return

def test_explode_dataframe():
    """Unit tests for explode_dataframe()
    """
    df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    exploded_df = metamoles.explode_dataframe(df, metamoles.parse_compound_ids,
                                    'product', ['entry'])
    assert exploded_df.shape == (11, 2)
    assert exploded_df['product'].tolist() == ['C00071','C00004','C00080','C01450',
        'C00071','C00005','C00080','C00441','C00004','C00005', 'C00080']
    return

def test_create_negative_matches():
    """Unit tests for create_negative_matches()
    """
    df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    exploded_df = metamoles.explode_dataframe(df, metamoles.parse_compound_ids,
                                    'product', ['entry'])
    pos_df, neg_df = metamoles.create_negative_matches(exploded_df, 'entry', 'product')
    assert pos_df.shape == (11, 3)
    assert neg_df.shape == (7, 3)
#    assert pos_df['enzyme'].tolist() == ['1.1.1.2', '1.1.1.2', '1.1.1.2', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.3', '1.1.1.3', '1.1.1.3', '1.1.1.3'], "ValueError enzyme listing is wrong"
#    assert neg_df['product'].tolist() == ['C00004', 'C01450', 'C00441', 'C00441', 'C00005', 'C01450', 'C00071'], "ValueError compound listing is wrong"
    return

def test_parse_reaction_ids():
    """Unit tests for parse_reaction_ids()
    """
    df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    assert metamoles.parse_reaction_ids(df)==['R00623', 'R00624', 'R07328', 'R01773', 'R01775']
    return

def test_parse_reversible_reactions():
     """Unit tests for parse_reversible_reactions()
     """
     df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
     reaction_list = metamoles.parse_reaction_ids(df)
     assert metamoles.parse_reversible_reactions(reaction_list) == reaction_list
     return


def test_binarize_enzyme_class():
     """Unit tests for binarize_enzyme_class()
     """
     df = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
     actual = metamoles.binarize_enzyme_class(df,column="entry").enzyme_class_1.tolist()
     assert actual == [1, 1, 1]
     return

def test_remove_single_cpd_rows():
     """Unit tests for remove_single_cpd_rows()
     """
     smilesdf = pd.read_csv("playground_df_cleaned_kegg_with_smiles.csv")
     assert metamoles.remove_single_cpd_rows(smilesdf,"entry","SMILES").shape == (50, 6)
     return

def test_join_pubchem_ids():
     """Unit tests for join_pubchem_ids()
     """
     pubchemdf = metamoles.create_kegg_df("test_kegg_compound_records.txt.gz","compound")
     masterdf =  metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
     return

def test_sid_to_smiles():
    """Unit test for pubchem_client.py sid_to_smiles."""
    sids = ['3489', '3990']
    expected = ['C(CO)N', 'C1CSSC1CCCCC(=O)O']
    actual = []
    for sid in sids:
        result_smile = sid_to_smiles(sid)
        assert len(
            result_smile) >= 1, 'SMILES string is very short. Check SMILES.'
        isinstance(result_smile, str), 'SMILES not returned as string.'
        actual.append(result_smile[0])
    assert expected == actual, 'Actual SMILES are not the expected SMILES.'
    return

def test_kegg_df_to_smiles():
    """Unit test for pubchem_client.py kegg_df_to_smiles."""
    test_frame = pd.DataFrame([['space fill', 'ethanolamine', '3489'], [
                              'space fill', 'pyruvate', '3324']], columns=['Filler', 'Compound Name', 'SID'])
    expected_frame = pd.DataFrame([[int(700),
                                    'C(CO)N',
                                    'space fill',
                                    'ethanolamine',
                                    '3489'],
                                   [int(1060),
                                    'CC(=O)C(=O)O',
                                    'space fill',
                                    'pyruvate',
                                    '3324']],
                                  columns=['CID',
                                           'SMILES',
                                           'Filler',
                                           'Compound Name',
                                           'SID'])
    result_frame = kegg_df_to_smiles(test_frame, 'SID')
    assert_frame_equal(
        result_frame[0], expected_frame), 'Did not generate expected df.'
    return
