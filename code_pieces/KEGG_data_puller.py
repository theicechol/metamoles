# This program is designed to pull all data from a specified KEGG database, and write the output
# to a csv file

# imports
from Bio.KEGG import REST
from Bio.KEGG import Enzyme

import pandas as pd
import argparse

# supported KEGG databases
KEGG_dbs = ['pathway', 'ko', 'genome', 'compound', 'glycan', 'reaction', 'enzyme', 'network']

# argument handling
parser = argparse.ArgumentParser()
parser.add_argument('database', help='database argument must belong to supported options: {}'.format(KEGG_dbs), type=str)

args = parser.parse_args()

print(args.database)
