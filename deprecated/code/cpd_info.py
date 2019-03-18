# create a function that return the general information of that compound e.g., number of each atom, MW, unsaturation
# this will be cpd_info.py
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt

def count_C(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
def count_O(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
def count_N(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
def count_P(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
def count_S(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
def count_X(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9 or atom.GetAtomicNum() == 17 or
        atom.GetAtomicNum() == 35 or atom.GetAtomicNum() == 53)
def count_H(mol):
    H = 0
    for i in range(mol.GetNumAtoms()):
        H += mol.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True)
    return H

def cpd_inform(SMILES):
    """
    A function for getting compound information from SMILES string it received
        a SMILES string and return a dictionary of information consisted of
        number of C, H, O , N, P, S, X, Degree of Unsaturation
        and Molecular Weight
    """
    info = []
    mol = Chem.rdmolfiles.MolFromSmiles(SMILES)
    info.append(float(count_C(mol)))
    info.append(float(count_H(mol)))
    info.append(float(count_O(mol)))
    info.append(float(count_N(mol)))
    info.append(float(count_P(mol)))
    info.append(float(count_S(mol)))
    info.append(float(count_X(mol)))
    info.append((2*info[0] + 2 + info[3] + info[4] - info[6] - info[1])/2) # it is (2*C + 2 + N + P - X - H)/2
    info.append(MolWt(mol))
    return info

# Create a function that create a new column of chemical information

def create_cpd_info(input_df, col_name='SMILES'):

    """
    Receive a DataFrame and return a dataframe with additional columns
        named n_C, n_H, ..., DoU, and MW
    """
    n_C = []
    n_H = []
    n_O = []
    n_N = []
    n_P = []
    n_S = []
    n_X = []
    DoU = []
    MW = []

    for _, row in input_df.iterrows():
        mol = row[col_name]
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

    input_df['n_C'] = n_C
    input_df['n_H'] = n_H
    input_df['n_O'] = n_O
    input_df['n_N'] = n_N
    input_df['n_P'] = n_P
    input_df['n_S'] = n_S
    input_df['n_X'] = n_X
    input_df['DoU'] = DoU
    input_df['MW'] = MW

    return input_df
