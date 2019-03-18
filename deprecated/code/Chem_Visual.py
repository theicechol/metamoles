# This cell will be Chem_Visual.py

from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

# Create mol object from smiles string or rdkit mol object
def Draw_mol(mol):
    """Input object is mol which can be either SMILES string or RDKit mol object
    They will be visualize in IPython Display"""
    if type(mol) == type('SMILES'):
        mol = Chem.MolFromSmiles(mol)
    else:
        pass

    molSize=(450,150)
    mc = Chem.Mol(mol.ToBinary())
        
    rdDepictor.Compute2DCoords(mc)

    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])

    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()
        
    return SVG(svg)

# To show set of molecules in grid or side-by-side

def Draw_Grid(mol, n_row = 0):
    """ Input a list of molecule to be written"""
    if n_row == 0:
        n_row = len(mol)
    grid = []
    for molecule in range(len(mol)):
        if type(mol[molecule]) == type('SMILES'):
            grid.append(Chem.MolFromSmiles(mol[molecule]))
        else:
            grid.append(mol[molecule])
    fig = Draw.MolsToGridImage( grid, molsPerRow=n_row)
    return fig

# create a function that make a retrosynthetic equation from two molecules
def Draw_Rxn_1to1(sm, pdt):
    """The input compounds are product (pdt) and starting material (sm) in SMILES string format"""
    if type(sm) == type('SMILES'):
        mol1 = Chem.MolFromSmiles(sm)
    else:
        mol1 = sm
    if type(pdt) == type('SMILES'):
        mol2 = Chem.MolFromSmiles(pdt)
    else:
        mol2 = pdt
    smart1 = Chem.MolToSmarts(mol1, isomericSmiles=True)
    smart2 = Chem.MolToSmarts(mol2, isomericSmiles=True)
    rxn_str = smart1 + '>>' + smart2
    rxn = AllChem.ReactionFromSmarts(rxn_str)
    return Draw.ReactionToImage(rxn)

def Draw_Rxn_2to1(sm1, sm2, pdt):
    """The input compounds are starting material 1 and 2 (sm), and product (pdt)"""
    if type(sm1) == type('SMILES'):
        mol1 = Chem.MolFromSmiles(sm1)
    else:
        mol1 = sm1
    
    if type(sm2) == type('SMILES'):
        mol2 = Chem.MolFromSmiles(sm2)
    else:
        mol2 = sm2

    if type(pdt) == type('SMILES'):
        mol3 = Chem.MolFromSmiles(pdt)
    else:
        mol3 = pdt
    
    smart1 = Chem.MolToSmarts(mol1, isomericSmiles=True)
    smart2 = Chem.MolToSmarts(mol2, isomericSmiles=True)
    smart3 = Chem.MolToSmarts(mol3, isomericSmiles=True)
    
    rxn_str = smart1 + '.' + smart2 + '>>' + smart3
    rxn = AllChem.ReactionFromSmarts(rxn_str)
    
    return Draw.ReactionToImage(rxn)
