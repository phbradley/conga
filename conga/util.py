import numpy as np
import sys
from os import system
import os.path
from scipy.sparse import issparse
from collections import Counter
#import pandas as pd
#from . import preprocess as pp

# convenience paths
path_to_conga = os.path.dirname(os.path.realpath(__file__))
assert not path_to_conga.endswith('/')
path_to_conga += '/'

path_to_data = path_to_conga+'data/'
assert os.path.isdir( path_to_data )


# def get_mean_var_scanpy(X):
#     ''' this is based on---ie stolen from---scanpy's _get_mean_var function
#     wanted to be explicit about what our variance convention is
#     '''
#     #
#     mean = X.mean(axis=0)
#     if issparse(X):
#         mean_sq = X.multiply(X).mean(axis=0)
#         mean = mean.A1
#         mean_sq = mean_sq.A1
#     else:
#         mean_sq = np.multiply(X, X).mean(axis=0)
#     # do we want unbiased estimator?
#     #var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
#     var = (mean_sq - mean**2)
#     return mean, var


def run_command( cmd, verbose=False ):
    if verbose:
        print('util.run_command: cmd=', cmd)
    system(cmd)

def make_clones_file( tcrs, outfilename, subject = 'UNK', epitope = 'UNK_E' ):
    ''' This may not have all the standard fields
    Right now just adding the fields we need in order for make_tcr_logo.py to work...
    '''
    gene_fields = ['{}{}_{}'.format(x,y,z) for x in 'vj' for y in 'ab' for z in ['gene', 'genes']]
    outfields = 'clone_id subject epitope cdr3a cdr3a_nucseq cdr3b cdr3b_nucseq'.split() + gene_fields

    out = open(outfilename, 'w')
    out.write('\t'.join(outfields)+'\n')
    for ii,(atcr, btcr) in enumerate(tcrs):
        outl = { 'clone_id': 'clone_{}'.format(ii+1),
                 'subject': subject,
                 'epitope': epitope,
                 'va_gene': atcr[0],
                 'va_genes': atcr[0],
                 'ja_gene': atcr[1],
                 'ja_genes': atcr[1],
                 'cdr3a': atcr[2],
                 'cdr3a_nucseq': atcr[3],
                 'vb_gene': btcr[0],
                 'vb_genes': btcr[0],
                 'jb_gene': btcr[1],
                 'jb_genes': btcr[1],
                 'cdr3b': btcr[2],
                 'cdr3b_nucseq': btcr[3],
                 }
        out.write('\t'.join( outl[x] for x in outfields)+'\n')
    out.close()


