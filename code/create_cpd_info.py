# Create a function that create a new column of chemical information

import pandas as pd
import numpy as np
from cpd_inform import cpd_inform

def create_cpd_info(input_df, col_name='SMILES'):

    """Receive a DataFrame and return a dataframe with additional columns named n_C, n_H, ..., DoU, and MW"""

    # create an empty funciton with either empty strint '' or NaN by np.nan

#     #df['n_C'] = ''
#     #df['n_H'] = ''
#     df['n_O'] = ''
#     df['n_N'] = ''
#     df['n_P'] = ''
#     df['n_S'] = ''
#     df['n_X'] = ''
#     df['DoU'] = ''
#     df['MW'] = ''
    n_C = []
    n_H = []
    n_O = []
    n_N = []
    n_P = []
    n_S = []
    n_X = []
    DoU = []
    MW = []
    # see df.iterrows()

    #for _, row in input_df.iterrows():
    for row in range(input_df.shape[0]):
        mol = input_df[col_name][row]
        info = cpd_inform(mol)
        n_C.append(info[0])
        n_H.append(info[1])
        n_O.append(info[2])
        n_N.append(info[3])
        n_P.append(info[4])
        n_S.append(info[5])
        n_X.append(info[6])
        DoU.append(info[7])
        MW.append(info[8])
#         df['n_C'][index] = info['n_C']
#         df['n_H'][index] = info['n_H']
#         df['n_O'][index] = info['n_O']
#         df['n_N'][index] = info['n_N']
#         df['n_P'][index] = info['n_P']
#         df['n_S'][index] = info['n_S']
#         df['n_X'][index] = info['n_X']
#         df['DoU'][index] = info['DoU']
#         df['MW'][index] = info['MW']

    input_df['n_C'] = pd.DataFrame(n_C)
    input_df['n_H'] = pd.DataFrame(n_H)
    input_df['n_O'] = pd.DataFrame(n_O)
    input_df['n_N'] = pd.DataFrame(n_N)
    input_df['n_P'] = pd.DataFrame(n_P)
    input_df['n_S'] = pd.DataFrame(n_S)
    input_df['n_X'] = pd.DataFrame(n_X)
    input_df['DoU'] = pd.DataFrame(DoU)
    input_df['MW'] = pd.DataFrame(MW)

    return input_df
