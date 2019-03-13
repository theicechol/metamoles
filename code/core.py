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
