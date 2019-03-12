# Create a function that create a new column of chemical information

import pandas as pd
from cpd_inform import cpd_inform

def create_cpd_info(df, col_name='SMILES'):
    
    """Receive a DataFrame and return a dataframe with additional columns named n_C, n_H, ..., DoU, and MW"""
    
    df['n_C'] = 'Empty'
    df['n_H'] = 'Empty'
    df['n_O'] = 'Empty'
    df['n_N'] = 'Empty'
    df['n_P'] = 'Empty'
    df['n_S'] = 'Empty'
    df['n_X'] = 'Empty'
    df['DoU'] = 'Empty'
    df['MW'] = 'Empty'
    
    for index in range(df.shape[0]):
        mol = df[col_name][index]
        info = cpd_inform(mol)
        df['n_C'][index] = info['n_C']
        df['n_H'][index] = info['n_H']
        df['n_O'][index] = info['n_O']
        df['n_N'][index] = info['n_N']
        df['n_P'][index] = info['n_P']
        df['n_S'][index] = info['n_S']
        df['n_X'][index] = info['n_X']
        df['DoU'][index] = info['DoU']
        df['MW'][index] = info['MW']
        
    return df
