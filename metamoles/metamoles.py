#!/usr/bin/env python

import scipy as sp
import gzip
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
import Bio
import rdkit
import re
from Bio.KEGG import Compound
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from sklearn import linear_model
from sklearn.model_selection import train_test_split


def create_kegg_df(file_path: str, kegg_db: str):
    """
    create_kegg_df() parses a gzipped text file of KEGG records using the
        Biopython.Bio.KEGG package, and returns the data as a pandas dataframe.

        NOTE: only 'enzyme' and 'compound' KEGG records are supported

    Args:
        file_path (str): filepath string pointing to gzipped text file
            of KEGG records
        kegg_db (str): either 'enzyme' or 'compound' keyword argument

    Returns:
        pandas.DataFrame: containing all Enzyme.Record() fields in columns
    """
    supported_dbs = ['enzyme', 'compound']

    if kegg_db == 'enzyme':
        parser = Enzyme
    elif kegg_db == 'compound':
        parser = Compound
    else:
        raise ValueError('supported kegg_db values include: {}'.format(supported_dbs))

    field_list = [method for method in dir(parser.Record()) if not method.startswith('_')]
    data_matrix = []

    with gzip.open(file_path, 'rt') as file:
        for record in parser.parse(file):
            data_matrix.append([getattr(record, field) for field in field_list])

    kegg_df = pd.DataFrame(data_matrix, columns=field_list)
    return kegg_df


def select_promiscuous_enzymes(enzyme_df: pd.DataFrame):
    """
    select_promiscuous_enzymes() down selects the enzymes from an input
        dataframe to include only those in which 2 or more reactions are
        associated with their KEGG record.

    Args:
        enzyme_df (pandas.DataFrame): must contain at least fields
            ['reaction', 'entry', 'product', 'substrate']

    Returns:
        pandas.DataFrame: pandas dataframe containing only the fields
            ['reaction', 'entry', 'product', 'substrate'] and containing only
            rows of enzymes with two or more reactions in their KEGG record
    """

    promiscuous_df = enzyme_df[[True if len(rxn) > 1 else False for rxn in enzyme_df['reaction']]]
    compact_promiscuous_df = promiscuous_df[['entry','reaction','product','substrate']]

    return compact_promiscuous_df


def parse_compound_ids(field: str):
    """
    parse_compound_ids() uses regular expressions to extract the KEGG compound
        IDs from a product or substrate field in a KEGG record field

    Args:
        field (str): name of field that contains KEGG compound IDs in a string

    Returns:
        list: contains parsed KEGG compound IDs
    """

    cpd_list = []
    regex = 'CPD:(C\d+)'
    # matches 'CPD:' chars exactly and captures 'C' + any following digits (\d+)
    for entry in field:
        ids = re.findall(regex, str(entry), re.IGNORECASE)
        for i in ids:
            cpd_list.append(i)

    return cpd_list


def parse_pubchem_ids(field: str):
    """
    parse_pubchem_ids() uses regular expressions to extract the PubChem
        compound IDs from a field in a record

    Args:
        field (str): name of a pandas.DataFrame field containing PubChem
            compound IDs in a string

    Returns:
        str: extracted pubchem_id
    """

    regex = "'PubChem', \[\'(\d+)\'\]\)"
    # matches "'PubChem', ['" characters exactly, then captures any following
    # digits (\d+), before another literal "']" character match

    ids = re.findall(regex, str(field), re.IGNORECASE)
    if len(ids) > 0:
        pubchem_id = ids[0]
    else:
        pubchem_id = ''

    return pubchem_id


def parse_reaction_ids(df: pd.DataFrame):
    """
    parse_reaction_ids() parses the list of reaction numbers from a
        dataframe containing a column labeled 'reaction'

    Args:
        df (pandas.DataFrame): must contain column 'reaction'

    Returns:
        list: contains parsed KEGG reaction IDs
    """
    reaction_list = []
    for index,row in df.iterrows():
        for reaction in row['reaction']:
            reaction_split = reaction.split("[RN:")[-1]
            if reaction_split.startswith("R") and not reaction_split.startswith("RN"):
                for i in reaction_split[:-1].split(" "):
                    reaction_list.append(i)
    # note - this approach is case sensitive and will miss reaction numbers
    # labeled with a lowercase r, or with incorrect spacing. A regex extraction
    # may be one option for improved robustness, though likely sacrifices speed.
    return reaction_list


def parse_reversible_reactions(reaction_list: list):
    """
    parse_reversible_reactions() queries the KEGG database with the input
        reaction list, and parses the results for all reactions that have been
        annotated with "<=>" in the reaction equation, which suggests that the
        catalyzed reaction is reversible

    Args:
        reaction_list (list): contains KEGG reaction IDs (e.g. 'R00709')

    Returns:
        list: contains KEGG IDs of reversible reactions
    """

    reversible_reaction = []
    for reaction in reaction_list:
        reaction_file = REST.kegg_get(reaction).read()
        for i in reaction_file.rstrip().split("\n"):
            if i.startswith("EQUATION") and "<=>" in i:
                reversible_reaction.append(reaction)
    return reversible_reaction


def combine_substrates_products(df: pd.DataFrame):
    """
    combine_substrates_products() is for use with a collection of enzymes
        in which it is understood that they are capable of catalyzing both the
        forward and reverser reactions. In this case, both the substrates and
        the products should be considered as bioreachable products.
        This function parses the list of substrates and products from their
        respective fields in the input dataframe, and returns a new dataframe
        with the combined substrates and products in a column labeled 'product'

    WARNING: combine_substrates_products() should not be run multiple times on
        the same dataframe becuase it will will append duplicate substrates

    Args:
        df (pandas.DataFrame): must contain fields
            ['entry', 'substrate', 'product']

    Returns:
        pandas.DataFrame: contains only fields ['entry', 'product']
    """

    rowindex = np.arange(0,len(df))
    df_with_ordered_index = df.set_index(rowindex)

    newdf = df_with_ordered_index
    # should this be a .copy()?

    for index,row in df_with_ordered_index.iterrows():
        productlist = row['product']
        substratelist = row['substrate']
        newdf.iloc[index,2] = productlist + substratelist

    return newdf[["entry","product"]]


def explode_dataframe(dataframe: pd.DataFrame, explosion_function,
                        explosion_target_field: str, fields_to_include: list):
    """
    explode_dataframe() applies the input explosion_function to the target
        field in each row of a dataframe. Each item in the output of the
        explosion_function is an anchor for a new row in the new dataframe. All
        of the supplied fields_to_include are added to the explosion item,
        and appended to the new dataframe row.

    Args:
        dataframe (pandas.DataFrame): input dataset
        explosion_function (function): function to be applied to target
            column in dataframe
        explosion_target_field (str): name of field in dataframe to which the
            explosion funciton will be applied
        fields_to_include (list): a list of strings that denote the columns of
            the input dataframe to be included in the output

    Returns:
        pandas.DataFrame: new exploded dataframe
    """
    new_rows = []
    for _, row in dataframe.iterrows():
        explosion_list = explosion_function(row[explosion_target_field])
        for item in explosion_list:
            row_data = [row[field] for field in fields_to_include]
            row_data.append(item)
            new_rows.append(row_data)

    fields_to_include.append(explosion_target_field)
    new_df = pd.DataFrame(new_rows, columns=fields_to_include)

    return new_df


def remove_cofactors(master_df: pd.DataFrame, master_cpd_field: str,
                     cofactor_df: pd.DataFrame, cofactor_field: str,
                     drop_na=True):
    """
    remove_cofactors() should be used to clean the dataset of cofactors. These
        will be included in the KEGG records as substrates and products, but
        are not actually products in the reaction

    Args:
        master_df (pandas.DataFrame): input dataset
        master_cpd_field (str): field that contains products
        cofactor_df (pandas.DataFrame): contains cofactors to be removed
        cofactor_field (str): field that contains cofactors
        drop_na (bool): default True

    Returns:
        pandas.DataFrame: cleaned data without cofactor entries
    """
    cofactor_list = parse_compound_ids(cofactor_df[cofactor_field])
    bool_mask = [False if cpd in cofactor_list else True for cpd in master_df[master_cpd_field]]
    clean_df = master_df[bool_mask]
    clean_df = clean_df.drop_duplicates()

    if drop_na:
        clean_df = clean_df[clean_df[master_cpd_field] != 'NA']
    else:
        pass

    return clean_df


def binarize_enzyme_class(dataframe, column):
    """
    binarize_enzyme_class() converts the enzyme class into binary dummy variables
        that are appended onto the input dataframe

    Args:
        dataframe (pandas.DataFrame): input dataset
        column (str): column name containing kegg enzyme id

    Returns:
        pandas.DataFrame: with seven columns appended for the seven enzyme classes
    """
    dataframe['enzyme_class'] = [row[column][0] for _, row in dataframe.iterrows()]
    dataframe = pd.get_dummies(dataframe, columns=['enzyme_class'])
    return dataframe


def create_negative_matches(dataframe: pd.DataFrame,
                            enzyme_field: str, compound_field: str):
    """
    create_negative_matches() returns two dataframes.
        One dataframe is positive data that contains all the enzyme-compound
        pairs that exist in the input dataset.
        The second data frame is negative data made from matching all
        enzyme-compound pairs that do not exist in the dataset.

    Args:
        dataframe (pandas.DataFrame): input dataset
        enzyme_field (str): column in dataframe that contains enzyme ids
        compound_field (str): column in dataframe that contains compound ids

    Returns:
        pandas.DataFrame: positive data
            (contains fields ['enzyme', 'product', 'reacts'])
        pandas.DataFrame: negative data
            (contains fields ['enzyme', 'product', 'reacts'])
    """
    unique_enzymes = set(dataframe[enzyme_field].unique())
    # set of all unique enzymes in provided dataframes
    unique_cpds = set(dataframe[compound_field].unique())
    # set of all unique compounds in provided dataframe

    positive_data = []
    negative_data = []
    # initialize empty lists

    for enzyme in unique_enzymes:
    # iterate through unique enzyme set
        working_prods = set(dataframe[dataframe[enzyme_field] == enzyme][compound_field].unique())
        # unique set of all products reported to reaction with this enzyme in provided dataset
        non_working_prods = (unique_cpds - working_prods)
        # set math of all remaining products in the dataset minus those reported to react

        reactions = [{'reacts':1.0, 'enzyme':enzyme, 'product':product} for product in working_prods]
        # create new entry for each positive reaction
        non_reactions = [{'reacts':0.0, 'enzyme':enzyme, 'product':product} for product in non_working_prods]
        # create new entry for each negative reaction

        positive_data.extend(reactions)
        # add positive reactions to master list
        negative_data.extend(non_reactions)
        # add negative reactions to master list

    positive_df = pd.DataFrame(positive_data)
    negative_df = pd.DataFrame(negative_data)

    return positive_df, negative_df


def remove_single_cpd_rows(dataframe, enzyme_col, smiles_col):
    """
    remove_single_cpd_rows() is meant to be a pre-processing function prior to passing a dataframe to the
        calculate_dist() function

    Args:
        dataframe (pandas.Dataframe): input dataset
        enzyme_col (str): name for column that contains kegg enzyme ids
        smiles_col (str): name for column that contains smiles string

    Returns:
        pandas.Dataframe: output dataframe with rows removed in which there was only one product paired with
            the enzyme entry, enzyme_col renamed 'entry', and smiles_col renamed 'SMILES'
    """
    dataframe = dataframe.rename(columns={enzyme_col:'entry', smiles_col:'SMILES'})
    counts_df = dataframe.groupby('entry').count()
    singles_df = counts_df[counts_df['SMILES'] == 1]
    singles = singles_df.index.tolist()
    bool_mask = [False if row['entry'] in singles else True for _, row in dataframe.iterrows()]
    clean_df = dataframe[bool_mask]
    return clean_df


def join_pubchem_ids(master_df, pubchem_df, master_join_key, pubchem_join_key,
                        pubchem_id_field):
    """
    join_pubchem_ids() takes an input dataframe containing a column of KEGG
        compound ids, and a second dataframe containing KEGG compound ids and
        their corresponding PubChem ids. The function parses the PubChem ids
        from the correct column, and joins these onto the input dataframe

    Args:
        master_df (pandas.DataFrame): input dataset
        pubchem_df (pandas.DataFrame): dataframe containing PubChem ids
        master_join_key (str): field in master_df with KEGG compound ids
        pubchem_join_key (str): field in pubchem_df with KEGG compound ids
        pubchem_id_field (str): field in pubchem_df with PubChem ids

    Returns:
         pandas.DataFrame:
    """
    pubchem_id_data = []

    for _, row in pubchem_df.iterrows():
        pubchem_id = parse_pubchem_ids(row[pubchem_id_field])
        join_key = row[pubchem_join_key]
        entry = {'pubchem_id': pubchem_id, pubchem_join_key: join_key}
        pubchem_id_data.append(entry)

    join_df = pd.DataFrame(pubchem_id_data)
    master_df = master_df.merge(join_df, left_on=master_join_key,
                                right_on=pubchem_join_key)
    master_df = master_df.drop(columns=pubchem_join_key)

    return master_df

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


def kegg_df_to_smiles(kegg_df, column_name):
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
        sid = kegg_df.loc[i, column_name]
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

    kegg_df.insert(0, column='CID', value=cid_list)
    # Change this 2 to the number where the smiles column should be
    kegg_df.insert(1, column='SMILES', value=res)
    # kegg_df.to_csv(r'../datasets/df_cleaned_kegg_with_smiles.csv')

    return kegg_df, unsuccessful_list

def input_data(input_df): #cleans input df and returns neccessary elements
    '''From the input dataframe, removes rows that do not contain product
    SMILES strings. Returns the cleaned dataframe'''
    for index, row in input_df.iterrows():

        if row['SMILES'] == 'none':

            input_df.drop(index, inplace=True)

    return input_df

def fingerprint_products(input_df): #fingerprints all products in a given df
    '''From the input dataframe, makes a list of rdkit Mol objects and makes a
    list of rdkit fingerprints generated from those Mol objects. Inserts both
    lists as new columns and returns the expanded dataframe.'''
    mol_list = []
    fp_list = []

    for index, row in input_df.iterrows():
        mol_list.append(Chem.rdmolfiles.MolFromSmiles(row['SMILES'])) #get mols from SMILES and add mols to list
        fp_list.append(FingerprintMols.FingerprintMol(Chem.rdmolfiles.MolFromSmiles(row['SMILES']))) #get fingerprints from mols and and fingerprints to list

    input_df['Mol'] = mol_list
    input_df['Fingerprint'] = fp_list

    return input_df

def sim_i_j(row_i, row_j):
    """For two given rows of a dataframe, use the rdkit fingerprints to compute
    TanimotoSimilarity and return the resulting float"""
    return DataStructs.FingerprintSimilarity(row_i['Fingerprint'], row_j['Fingerprint'], metric=DataStructs.TanimotoSimilarity)

def sim_i_all(input_df, index_i, row_i, metric):
    """From the input dataframe, check the passed indexes against the DataFrame,
    and construct a new dataframe which is the similarity matrix of all of the
    products contained in the dataframe."""
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
    """From an input_df, use sim_i_j and sim_i_all to build and return a
    similarity matrix dataframe."""
    metric = pd.DataFrame()
    for index_i, row_i in input_df.iterrows():
        sim_i_all(input_df, index_i, row_i, metric)
    return metric

def calculate_dist(input_df):
    '''Main method, takes an input dataframe and builds and returns a master
    dataframe which is the original dataframe, with three additional columns,
    an rdkit Mol column, an rdkit Fingerprint column, and a column which
    describes the average distance of a product row to all the products of the
    associated enzyme entry. Requires the KEGG enzyme entry column to be named 'entry'
	and the SMILES string column to be named 'SMILES' '''

    master_df = fingerprint_products(input_data(input_df))    #expand input df: generate mols from SMILES then generate fingerprints from mols, adding columns for each
    # enzyme_df_list = split_by_enzyme(input_df)    #split expanded df by rows, grouped by enzyme entry (1.1.1.110 etc), into a list of dataframes
    unique_enzymes = set(master_df['entry'].unique()) # create set of unique enzymes
    dist_lookup = {} # initialize master dist list
    for enzyme in unique_enzymes:    #loop through list of enzyme dataframes
        # enzyme_df['Dist'] = '' #initialize distance column
        enzyme_df = master_df[master_df['entry'] == enzyme]
        metric = sim_metric(enzyme_df) #get similarity matrix dataframe
        vals = metric.values #use np array of similarity matrix
        start_at = 1 #skip autocorrelation
        dist_list =[] #initialize list
        for i in range(len(vals)-1): #row of matrix except for last row
            for j in range(start_at, len(vals)): #col of matrix skipping first column
                dist_list.append(vals[i][j]) #add distance value to list
            start_at += 1 #start at higher index to skip redundancy
        avg_dist = sum(dist_list)/len(dist_list) #compute average distance
        dist_lookup[enzyme] = avg_dist
        # for _, row in enzyme_df.iterrows():    #loop through enzyme dataframe
        #     # enzyme_df['Dist'].loc[index] = avg_dist #add averaged distance to each product row of enzyme dataframe
    master_df['dist'] = [dist_lookup[row['entry']] for _, row in master_df.iterrows()]
    return master_df

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

def query_model(master_df, query_sid):
    """
    NOTE: Fields containing enzyme, compound PubChem sid, and SMILES string must be named
        ['entry', 'PubChem', 'SMILES'] respectively
    """
    # get query SMILES string & pair query compound with each unique enzyme in the master DataFrame
    updated_df = pair_query_compound(master_df, 'entry', 'PubChem', 'SMILES', query_sid)
    # calculate molecular distances between products of the same enzyme
    distance_df = calculate_dist(updated_df)
    # remove any rows that are not the query compound
    reduced_df = distance_df[distance_df['PubChem'] == query_sid]
    # get dummy variables to represent enzyme class
    query_df = binarize_enzyme_class(reduced_df, 'entry')
    query_df = query_df.reset_index(drop=True)
    # add in compound features with RDKit
    cpd_query_df = create_cpd_info(query_df)

    # re-instantiate log reg model
    feature_df = master_df[['dist', 'enzyme_class_1', 'enzyme_class_2', 'enzyme_class_3',
           'enzyme_class_4', 'enzyme_class_5', 'enzyme_class_6', 'enzyme_class_7',
           'n_O', 'n_N', 'n_S', 'n_X', 'DoU']]
    features = np.array(feature_df) #shape balance array for regression
    reactions = list(master_df['reacts'])
    feature_train, feature_test, reaction_train, reaction_test = train_test_split(features, reactions,
                                                      test_size=0.20, random_state=42)
    model_1 = linear_model.LogisticRegression(solver='liblinear', penalty='l1', random_state=1, class_weight='balanced')
    model_1.fit(feature_train, np.ravel(reaction_train))

    # select query features
    query_feat_df = query_df[['dist', 'enzyme_class_1', 'enzyme_class_2', 'enzyme_class_3',
           'enzyme_class_4', 'enzyme_class_5', 'enzyme_class_6', 'enzyme_class_7',
           'n_O', 'n_N', 'n_S', 'n_X', 'DoU']]
    # query reactive enzymes
    predictions = model_1.predict(query_feat_df)
    pred = model_1.predict_proba(query_feat_df)

    # write results to a DataFrame
    prediction_values = pd.DataFrame(pred)
    model_descriptive_df = pd.DataFrame()
#     model_descriptive_df['0']=prediction_values[0]
    model_descriptive_df['p_reacts']=prediction_values[1]
    prediction_df = pd.merge(model_descriptive_df, query_df, left_index=True, right_index=True)
    # sort DataFrame
    prediction_df = prediction_df.sort_values(by=['p_reacts'], ascending=False)
    # reset index in output dataframe
    prediction_df = prediction_df.reset_index(drop=True)
    # add rank to dataframe
    prediction_df['rank'] = prediction_df.index + 1
    # return DataFrame
    return prediction_df
