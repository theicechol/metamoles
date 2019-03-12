# imports

def test_enzyme_records_to_df():
    """Unit tests for enzyme_records_to_df()
    """
    df = enzyme_records_to_df(test_kegg_enzyme_records.txt.gz)
    assert df.shape == (3, 16)
    assert df.columns == ['classname', 'cofactor', 'comment', 'dblinks',
        'disease', 'effector', 'entry', 'genes', 'inhibitor', 'name',
        'pathway', 'product', 'reaction', 'structures', 'substrate', 'sysname']
    assert df['entry'].tolist() == ['1.1.1.1', '1.1.1.2', '1.1.1.3']
    return

def test_compound_records_to_df():
    """Unit tests for compound_records_to_df()
    """
    df = test_kegg_compound_records.txt.gz
    assert df.shape == (55, 8)
    assert df.columns == ['dblinks', 'entry', 'enzyme', 'formula', 'mass',
        'name', 'pathway', 'structures']
    assert df['entry'].tolist() == ['C00001', 'C00002', 'C00003', 'C00004',
        'C00005', 'C00006', 'C00007', 'C00008', 'C00009', 'C00010', 'C00011',
        'C00012', 'C00013', 'C00014', 'C00015', 'C00016', 'C00017', 'C00018',
        'C00019', 'C00020', 'C00021', 'C00022', 'C00023', 'C00024', 'C00025',
        'C00026', 'C00027', 'C00028', 'C00029', 'C00030', 'C00031', 'C00032',
        'C00033', 'C00034', 'C00035', 'C00036', 'C00037', 'C00038', 'C00039',
        'C00040', 'C00041', 'C00042', 'C00043', 'C00044', 'C00045', 'C00046',
        'C00047', 'C00048', 'C00049', 'C00050', 'C00051', 'C00052', 'C00053',
        'C00054', 'C00055']
    return

def test_extract_KEGG_compound_IDs():
    """Unit tests for extract_KEGG_compound_IDs()
    """
    df = enzyme_records_to_df(test_kegg_enzyme_records.txt.gz)
    assert extract_KEGG_compound_IDs(df['product']) == ['C00071', 'C00004',
        'C00080', 'C01450', 'c00071', 'C00071', 'C00005', 'C00080',
        'C00441', 'C00004', 'C00005', 'C00080']
    return

def test_extract_PubChem_id():
    """Unit tests for extract_PubChem_id()
    """
    df = enzyme_records_to_df(test_kegg_enzyme_records.txt.gz)
    PubChemID_list = []
    for _, row in compound_df.iterrows():
        pubchem_id = extract_PubChem_id(row['dblinks'])
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
    df = enzyme_records_to_df(test_kegg_enzyme_records.txt.gz)
    exploded_df = explode_dataframe(df, extract_KEGG_compound_IDs,
                                    'product', ['entry'])
    assert exploded_df.shape == (12, 2)
    assert exploded_df['product'].tolist() == ['C00071', 'C00004', 'C00080',
        'C01450', 'c00071', 'C00071', 'C00005', 'C00080', 'C00441', 'C00004',
        'C00005', 'C00080']
    return

def test_neg_data_matchmaker():
    """Unit tests for neg_data_matchmaker()
    """
    df = enzyme_records_to_df(test_kegg_enzyme_records.txt.gz)
    exploded_df = explode_dataframe(df, extract_KEGG_compound_IDs,
                                    'product', ['entry'])
    pos_df, neg_df = neg_data_matchmaker(exploded_df, 'entry', 'product')
    assert pos_df.shape == (12, 3)
    assert neg_df.shape == (9, 3)
    assert pos_df['enzyme'].tolist() == ['1.1.1.2', '1.1.1.2', '1.1.1.2',
        '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.3',
        '1.1.1.3', '1.1.1.3', '1.1.1.3']
    assert neg_df['product'].tolist() == ['C01450', 'C00441', 'C00004',
        'c00071', 'C00441', 'C00005', 'C01450', 'C00071', 'c00071']
    return

def test_func():
    """Unit tests for func()
    """

    return

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
