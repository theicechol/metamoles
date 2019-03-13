# create a function that return the general information of that compound e.g., number of each atom, MW,
# unsaturation

from rdkit import Chem
from rdkit.Chem import AllChem

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

# DoU funciton is dropped out because of compoutational workload
def DoU(mol):
    """Calculation of Degree of Unsaturation"""
    C = count_C(mol)
    H = count_H(mol)
    O = count_O(mol)
    N = count_N(mol)
    P = count_P(mol)
    S = count_S(mol)
    X = count_X(mol)
    return float((2*C+2+N+P-X-H)/2)

from rdkit.Chem.Descriptors import MolWt

def cpd_inform(SMILES):

    """A function for getting compound information from SMILES string
    it received a SMILES string and return a dictionary of information consisted of number of C, H, O , N, P, S,
X, Degree of Unsaturation and Molecula) Weight"""
#     info = {}
    info = []
    mol = Chem.MolFromSmiles(SMILES)
#     info['n_C'] = float(count_C(mol))
#     info['n_H'] = float(count_H(mol))
#     info['n_O'] = float(count_O(mol))
#     info['n_N'] = float(count_N(mol))
#     info['n_P'] = float(count_P(mol))
#     info['n_S'] = float(count_S(mol))
#     info['n_X'] = float(count_X(mol))
#     info['DoU'] = (2*info['n_C']+2+info['n_N']+info['n_P']-info['n_X']-info['n_H'])/2
#     #info['DoU'] = DoU(mol)
#     info['MW'] = MolWt(mol)
    info.append(float(count_C(mol)))
    info.append(float(count_H(mol)))
    info.append(float(count_O(mol)))
    info.append(float(count_N(mol)))
    info.append(float(count_P(mol)))
    info.append(float(count_S(mol)))
    info.append(float(count_X(mol)))
    info.append((2*info[0] + info[3] + info[4] - info[6] - info[1])/2)
#    info['DoU'] = (2*info['n_C']+2+info['n_N']+info['n_P']-info['n_X']-info['n_H'])/2
#     #info['DoU'] = DoU(mol)
    info.append(MolWt(mol))
#     info['MW'] = MolWt(mol)
    return info
