
import pandas as pd
import pubchempy as pc

import kegg_data

def test_sid_to_smiles():
    
    sids = ['3489', '3990']
    expected = ['C(CO)N', 'C1CSSC1CCCCC(=O)O']
    actual = []
    
    for sid in sids:
        result_smile = kegg_data.sid_to_smiles(sid)
        
        assert len(result_smile) >= 1, 'SMILES string is very short. Check SMILES.'
        isinstance(result_smile, str), 'SMILES not returned as string.'
        
        actual.append(result_smile)
    
    assert expected == actual, 'Actual SMILES are not the expected SMILES.'
    
    return
   
    
    
def test_kegg_df_to_smiles():
    
    test_frame = pd.DataFrame([['ethanolamine', '3489'], ['pyruvate', '3324']], columns=['Compound Name', 'SID'])
    
    expected_frame = pd.DataFrame([['ethanolamine', '3489', 'C(CO)N'], ['pyruvate', '3324', 'CC(=O)C(=O)O']], columns=['Compound Name', 'SID', 'SMILES'])
    
    result_frame = kegg_data.kegg_df_to_smiles(test_frame)
    
    assert result_frame.equals(expected_frame), 'Did not generate expected df.'
    
    return