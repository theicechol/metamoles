#rdkit imports
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

#housekeeping imports
import pandas as pd
import matplotlib
import numpy as np
import scipy as sp


def input_data(input_df): #cleans input df and returns neccessary elements
    '''DocString'''
    
    for index, row in input_df.iterrows():
        
        if row['SMILES'] == 'none':
            
            input_df.drop(index, inplace=True)            
        
    return input_df
	
def fingerprint_products(input_df): #fingerprints all products in a given df
    '''DocString'''
    
    mol_list = []
    fp_list = []
    
    for index, row in input_df.iterrows():
        mol_list.append(Chem.rdmolfiles.MolFromSmiles(row['SMILES'])) #get mols from SMILES and add mols to list
        fp_list.append(FingerprintMols.FingerprintMol(Chem.rdmolfiles.MolFromSmiles(row['SMILES']))) #get fingerprints from mols and and fingerprints to list
        
    input_df.insert(6, column='Mol', value=mol_list)
    input_df.insert(7, column='Fingerprint', value= fp_list)
            
    return input_df
	
def split_by_enzyme(input_df):
    '''DocString'''
    
    unique_enzymes = set(input_df['entry'].unique())
    
    enzyme_df_list = []
    
    for entry in unique_enzymes: #for each unique enzyme in the input dataframe...
        
        enzyme_df = pd.DataFrame(columns=input_df.columns) #...initialize a new dataframe with the same columns as the input dataframe...
        
        for index, row in input_df.iterrows(): #...iterate through the input dataframe...
            
            if row['entry'] == entry: #... and add product rows that correspond to the unique enzyme entry...
                enzyme_df.loc[index] = row
                
        enzyme_df_list.append(enzyme_df) #...then add the completed dataframe of unique enzyme products to a list
           
    return enzyme_df_list #return list of dataframes
	
def sim_i_j(row_i, row_j):
    """DocString"""
    return DataStructs.FingerprintSimilarity(row_i['Fingerprint'], row_j['Fingerprint'], metric=DataStructs.TanimotoSimilarity)
	
def sim_i_all(input_df, index_i, row_i, metric):
    """DocString"""
    for index_j, row_j in input_df.iterrows():
        if index_j < index_i: #skip redundant rows
            continue
        elif index_i == index_j: #autocorrelate rows
            metric.loc[index_i, index_j] = 1
        else:
            metric.loc[index_i, index_j] = sim_i_j(row_i, row_j) #fill matrix with calculated similarity at two positions at once 
            metric.loc[index_j, index_i] = metric.loc[index_i, index_j]
    return 
	
def sim_metric(input_df):
    """DocString"""
    metric = pd.DataFrame()
    for index_i, row_i in input_df.iterrows():
        sim_i_all(input_df, index_i, row_i, metric)
    return metric
	
def main(input_df):
    '''DocString'''
        
    input_df = fingerprint_products(input_data(input_df))    #expand input df: generate mols from SMILES then generate fingerprints from mols, adding columns for each
    
    enzyme_df_list = split_by_enzyme(input_df)    #split expanded df by rows, grouped by enzyme entry (1.1.1.110 etc), into a list of dataframes
    
    for enzyme_df in enzyme_df_list:    #loop through list of enzyme dataframes
        
        enzyme_df['Dist'] = '' #initialize distance column
        
        metric = sim_metric(enzyme_df) #get similarity matrix dataframe
        
        vals = metric.values #use np array of similarity matrix
        
        start_at = 1 #skip autocorrelation
        
        dist_list =[] #initialize list
        
        for i in range(len(vals)-1): #row of matrix except for last row
            
            for j in range(start_at, len(vals)): #col of matrix skipping first column
                
                dist_list.append(vals[i][j]) #add distance value to list
            
            start_at += 1 #start at higher index to skip redundancy
        
        avg_dist = sum(dist_list)/len(dist_list) #compute average distance
        
        for index, row in enzyme_df.iterrows():    #loop through enzyme dataframe 
            enzyme_df['Dist'].loc[index] = avg_dist #add averaged distance to each product row of enzyme dataframe
    
    master_df = pd.concat(enzyme_df_list) #concatenate enzyme dataframes into master_df
    
    return master_df