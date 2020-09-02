## LAZY-- some basic imports and shared params
##
##
from glob import glob
from os import popen, system, chdir, remove, getcwd, mkdir
from os.path import exists, isdir, isfile
import math
from math import floor,sqrt
from sys import stderr,argv,exit
import random

from ..convert_svg_to_png import convert_svg_to_png, MONOSPACE_FONT_FAMILY

# combo dbs for tr and ig
db_file = 'combo_xcr.tsv'
#db_file = 'combo_db.tsv'
#db_file = 'alphabeta_db.tsv'

# should we use the crazy complicated tcr-dist scheme for picking a 'representative' gene for counting?
CLASSIC_COUNTREPS = False

## naming scheme for the gene segment types, occasionally useful for iterating
segtypes_uppercase = ['VA','JA','VB','JB']
segtypes_lowercase = ['va','ja','vb','jb']






# def Log(s): ## silly legacy helper function
#     stderr.write(s)
#     if s and not s.endswith('\n'):
#         stderr.write('\n')


