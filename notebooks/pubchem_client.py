
import numpy as np
import pandas as pd
import pubchempy as pc


def sid_to_smiles(sid):
    """Takes a PubChem SID. Returns the associated isomeric SMILES string and PubChem CID.

    Args:
        sid : The PubChem SID number.

    Returns:
        str: isomeric smiles.
        int: Pubchem CID number.

    """

    substance = pc.Substance.from_sid(sid)
    cid = substance.standardized_cid
    compound = pc.get_compounds(cid)[0]

    return compound.isomeric_smiles, cid


def kegg_df_to_smiles(kegg_df):
    """
    Args:
        kegg_df : pandas dataframe with SID numbers in the third column

    Returns:
        kegg_df : modified with a fourth column containing CID and fifth column containing SMILES
        unsuccessful_list : list of SIDs for which no CID or SMILES were found

    """

    res = []
    cid_list = []
    unsuccessful_list = []

    for i in range(len(kegg_df)):
        # cell index of desired SID
        sid = kegg_df.iloc[i, 2]
        try:
            smile_result = sid_to_smiles(sid)[0]
            res.append(smile_result)
            cid_result = sid_to_smiles(sid)[1]
            cid_list.append(cid_result)
        except BaseException:
            res.append('none')
            cid_list.append('none')
            unsuccessful_list.append(sid)
            pass

    kegg_df.insert(3, column='CID', value=cid_list)
    # Change this 2 to the number where the smiles column should be
    kegg_df.insert(4, column='SMILES', value=res)
    # kegg_df.to_csv(r'../datasets/df_cleaned_kegg_with_smiles.csv')

    return kegg_df, unsuccessful_list