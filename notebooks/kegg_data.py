
import pandas as pd
import pubchempy as pc

def sid_to_smiles(sid):
    """Takes an SID and prints the associated SMILES string."""
    
    substance = pc.Substance.from_sid(sid)
    cid = substance.standardized_cid
    compound = pc.get_compounds(cid)[0]
    return compound.isomeric_smiles


def kegg_df_to_smiles(kegg_df):
    """Takes a pandas dataframe that includes a column of SIDs, gets the isomeric SMILES for each SID, stores them as a list, then adds a SMILES column."""

    res = [] 
    
    for i in range(len(kegg_df)):
        sid = kegg_df.iloc[i, 1] #CHANGE THIS 1 TO THE PROPER COLUMN NUMBER FOR SID 
        result = sid_to_smiles(sid)
        res.append(result)
        
    
    kegg_df.insert(2, column='SMILES', value=res) #Change this 2 to the number where the smiles column should be
    
    return kegg_df