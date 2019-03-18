# This program is designed to pull all data from a specified KEGG database
# and write the output to a text file

# imports
from Bio.KEGG import REST

import pandas as pd
import argparse

# currently supported KEGG databases
supported_KEGG_dbs = ['genome', 'compound', 'reaction', 'enzyme']

# argument parser
parser = argparse.ArgumentParser()
parser.add_argument('database',
                    type=str,
                    choices=supported_KEGG_dbs,
                    help='currently supported KEGG databases: {}'.format(supported_KEGG_dbs))
parser.add_argument('-e', '--entrylist',
                    help='output csv of all database entries',
                    action='store_true')
args = parser.parse_args()

# data runner that fetches a list of all entries in specified database
def fetch_entry_list(database):
    "connects to appropriate KEGG database and fetches list of all entries"
    all_entries_df = pd.read_csv(REST.kegg_list(database),
                                sep='\t',
                                header=None,
                                names=['id', 'description'])
    return all_entries_df

# main function
def main():
    "main code block"
    # parse databse argument as variable
    db = args.database

    # get list of all entries in database
    entries_df = fetch_entry_list(db)
    id_list = entries_df.id
    handle = 'KEGG_{}_db_entries'.format(db)

    # if command line option selected, output csv of all database entries
    if args.entrylist:
        print('\nwriting list of {} database entries to csv file'.format(db))
        entries_df.to_csv(handle + '.csv')

    # fetch all records from specified KEGG database
    total_entries = len(id_list)
    print('\nbeginning retrieval of {} records from KEGG {} database\n'.format(total_entries, db))

    # write records to output txt file as being fetched
    counter = total_entries
    with open(handle + '.txt', 'a') as f:
        for id in id_list:
            print('retrieving record {} of {}: {}'.format(total_entries - counter + 1, total_entries, id))
            f.writelines(REST.kegg_get(id))
            counter -= 1

# necessary underscore stuff to make everything work - don't understand this yet really
if __name__ == "__main__":
    main()
