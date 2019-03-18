# function to query the SMILES string and append new pairs to the master dataframe

def pair_query_compound(master_df, enzyme_col, pubchem_col, smiles_col, pubchem_sid):
    """
    pair_query_compound_with_enzymes() queries pubchem to get a SMILES string from an input pubchem_sid,
        then pairs that query compound with each unique enzyme id in the master dataframe
        
    Args:
        master_df (pandas.DataFrame): master dataframe containing enzyme ids
        enzyme_col (str): column containing enzyme id
        pubchem_col (str): column containing pubchem sid
        smiles_col (str): column containing SMILES string
        pubchem_sid (str): query PubChem sid
        
    Returns:
        pandas.DataFrame: with rows added to include query compound
    """
    master_df = master_df[[enzyme_col, pubchem_col, smiles_col]]
    new_pairs = []
    smiles, _ = sid_to_smiles(pubchem_sid)
    if len(smiles) == 0:
        raise 'query compound SMILES string could not be retrieved'
    else:
        pass
    unique_enzymes = master_df[enzyme_col].unique().tolist()
    for enzyme in unique_enzymes:
        pair = {enzyme_col:enzyme, pubchem_col:pubchem_sid, smiles_col:smiles}
        new_pairs.append(pair)
    new_paris_df = pd.DataFrame(new_pairs)
    output_df = pd.concat((master_df, new_paris_df), axis=0, sort=False)
    return output_df