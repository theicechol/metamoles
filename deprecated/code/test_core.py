# imports
import metamoles
from metamoles import *

def test_create_kegg_df():
    """Unit tests for create_kegg_df()
    """

    df1 = metamoles.create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    df2 = metamoles.create_kegg_df("test_kegg_compound_records.txt.gz","compound")

    # for enzyme database
    assert df1.shape == (3, 16)
    assert df1.columns == ['classname', 'cofactor', 'comment', 'dblinks',
        'disease', 'effector', 'entry', 'genes', 'inhibitor', 'name',
        'pathway', 'product', 'reaction', 'structures', 'substrate', 'sysname']
    assert df1['entry'].tolist() == ['1.1.1.1', '1.1.1.2', '1.1.1.3']

    # for compound database

    assert df2.shape == (55, 8)
    assert df2.columns == ['dblinks', 'entry', 'enzyme', 'formula', 'mass',
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

    df = metamoles.test_create_kegg_df("test_kegg_enzyme_records.txt.gz","enzyme")
    test_prom_df = metamoles.test_select_promiscuous_enzymes(df)
    expected_column = ['reaction', 'entry', 'product', 'substrate']
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
    assert pos_df['enzyme'].tolist() == ['1.1.1.2', '1.1.1.2', '1.1.1.2',
        '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.3',
        '1.1.1.3', '1.1.1.3', '1.1.1.3']
    assert neg_df['product'].tolist() == ['C00441', 'C01450', 'C00004', 'C00441',
     'C00005', 'C01450', 'C00071']
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


binarize_enzyme_class(df,column="entry").enzyme_class_1.tolist()




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
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
# def test_func():
#     """Unit tests for func()
#     """
#
#     return
#
