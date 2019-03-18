#!/usr/bin/env python
import Bio
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
import re
from Bio.KEGG import Compound

import gzip
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def create_enzyme_df(path_to_file):
    """
    input:path_to_file. file.gz format
    output:enzyme dataframe
    """

    enzyme_fields = [method for method in dir(Enzyme.Record()) if not method.startswith('_')]
    data_matrix = []

    with gzip.open(path_to_file, 'rt') as file:
        for record in enzyme.parse(file):
            data_matrix.append([getattr(record, field) for field in enzyme_fields])

    enzyme_df = pd.DataFrame(data_matrix, columns=enzyme_fields)
    return enzyme_df



def get_compact_promiscuous_df(enzyme_df):
    """
    input:enzyme dataframe (dataframe)
    output:promiscuous enzyme dataframe (dataframe)
    """

    promiscuous_df = enzyme_df[[True if len(rxn) > 1 else False for rxn in enzyme_df['reaction']]]
    compact_promiscuous_df = promiscuous_df[['entry','reaction','product','substrate']]

    return compact_promiscuous_df



def get_reaction_list(df):
    """
    get the list of reaction from a dataframe that contains reaction column
    input:dataframe with reaction column (df)
    output: list of reaction (list)
    """
    reaction_list = []
    for index,row in df.iterrows():
        for reaction in row['reaction']:
            reaction_split = reaction.split("[RN:")[-1]
            if reaction_split.startswith("R") and not reaction_split.startswith("RN"):
                for i in reaction_split[:-1].split(" "):
                    reaction_list.append(i)
    return reaction_list



def query_reversible_reaction(reaction_list):
    """
    get the list of reversible reaction
    input:list of reactions(list) eg)["R00709"]
    output:list of reversible reactions(list) 
    """

    reversible_reaction = []
    for reaction in reaction_list:
        reaction_file = REST.kegg_get(reaction).read()
        for i in reaction_file.rstrip().split("\n"):
            if i.startswith("EQUATION") and "<=>" in i:
                reversible_reaction.append(reaction)
    return reversible_reaction



def combine_substrate_product(df):
    """
    append substrates to product column.
    should not be run multiple times. 
    it will append substrates multiple times
    input:dataframe with substrate and product(df)
    output:dataframe with combined substrate and product. named under product column(df)
    """

    rowindex = np.arange(0,len(df))
    df_with_ordered_index = df.set_index(rowindex)

    newdf = df_with_ordered_index
    
    for index,row in df_with_ordered_index.iterrows():
        productlist = row['product']
        substratelist = row['substrate']
        newdf.iloc[index,2] = productlist + substratelist 

    return newdf[["entry","product"]]



def get_cofactor_list(cofactor_df,CPDcolumnname):
    """
    <input>
    cofactor_df : cofactor dataframe(df)
    CPDcolumnname : name of CPD columnname from cofactor dataframe(str) 
    <output>
    cofactor_list : list of cofactors from cofactor dataframe (list)
    """

    cofactor_list = [cofactor[4:10] for cofactor in cofactor_df[CPDcolumnname]]
    return cofactor_list


def get_cpd_id(compound_full):
    """
    input:compound_full = compound full name (str) eg) 'oxalureate [CPD:C00802]'
    output: cpd =  cpd id (str) eg) 'C01007'
    """
    cpd = compound_full[-7:-1]
    return cpd 



def rm_cofactor_only_cpd(enzyme_df,cofactor_list,compound_columnname="product",keepNA=True):
    """
    <input>
    enzyme_df : dataframe with enzyme information. should have substrate and product combined(df)
    compound_columnname : name of the column with compounds (str)
    cofactor_list : list of cofactors to be removed (list)
    keepNA : if false, will drop the row with no compounds (boolean, default:True) 
    <output>
    clean dataframe (df) 
    """
    newdf = enzyme_df.drop(["product"],axis=1)
    cleaned_compound_column = []
    for index,row in enzyme_df.iterrows():
        cpd_compound_list =[]
        for compound in row[compound_columnname]:
            if "CPD" in compound:
                onlycpd = get_cpd(compound)
                if onlycpd not in cofactor_list:
                    cpd_compound_list.append(onlycpd)
                else:
                    pass
        if len(cpd_compound_list)==0:
            cleaned_compound_column.append("NA")
        else: 
            cleaned_compound_column.append(cpd_compound_list)
    newdf['product'] = cleaned_compound_column

    if keepNA==False:
        newdf = newdf.loc[cleaned_df_productinList['product']!='NA']
    
    return newdf



def itemlist_eachrow(df,oldcolumnname,newcolumnname,sorting_column):
    """
    <input>
    df: dataframe with list items in one column (dataframe)
    oldcolumnname : name of the old column to be replaced (str) eg)"products"
    newcolumnname : name of the new column to replace (str) eg)"product"
    sorting_column : name of the column to be sorted by (str) eg)"entry"

    <output>
    dataframe with each item in each row. 
 
    """
    newdf = df[oldcolumnname].\
    apply(pd.Series).\
    merge(df, left_index=True, right_index=True).\
    drop([oldcolumnname],axis=1).\
    melt(id_vars=[enzymecolumn],value_name=newcolumnname).\
    sort_values(by=[sorting_column]).\
    dropna().\
    drop(columns=["variable"])
    return newdf


def compound_records_to_df(file_path):
    """
    Function parses all records using Biopython.Bio.KEGG.Compound parser, and returns a pandas dataframe.
    <Input>
    filepath = file path to a gzipped text file of KEGG enzyme records (str) 
    <output>
    compound dataframe 
    """
    compound_fields = [method for method in dir(Compound.Record()) if not method.startswith('_')]
    data_matrix = []

    with gzip.open(file_path, 'rt') as file:
        for record in Compound.parse(file):
            data_matrix.append([getattr(record, field) for field in compound_fields])
    
    compound_df = pd.DataFrame(data_matrix, columns=compound_fields)
    return compound_df



def extract_PubChem_id(field):
    """
    This function uses regular expressions to extract the PubChem compound IDs from a field in a record
    input : field 
    output : pubchem_id 
    """

    regex = "'PubChem', \[\'(\d+)\'\]\)" # matches "'PubChem', ['" characters exactly, then captures any number of digits (\d+), before another literal "']" character match
    ids = re.findall(regex, str(field), re.IGNORECASE)
    if len(ids) > 0:
        pubchem_id = ids[0]
    else:
        pubchem_id = ''
    
    return pubchem_id




