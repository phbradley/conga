import argparse
import sys
import os
#from os import system
from os.path import exists
import numpy as np

sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # so we can import conga
import conga
from conga.tcrdist.tcr_distances import TcrDistCalculator
import pandas as pd




if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--clones_file')
    parser.add_argument('--organism', choices=['mouse', 'human'], required=True)
    parser.add_argument('--n_components', type=int, default=50)

    args = parser.parse_args()
    assert exists(args.clones_file)

    compute_tcrdist_kpcs_from_clones_file(args.clones_file, args.organism, args.n_components)

