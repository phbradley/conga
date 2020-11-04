from os.path import exists
from pathlib import Path
import numpy as np
import pandas as pd
from . import util


amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

imhc_model_file = Path.joinpath( Path(util.path_to_data),'logreg_hobit_donor1_v4.tsv_0_model.tsv')
assert exists(imhc_model_file)
imhc_model_df = pd.read_csv(imhc_model_file, sep='\t')
imhc_model_df.set_index('feature', inplace=True)

# restrict to features with nonzero weight
mask = (np.abs(imhc_model_df['coef'])>1e-6)
imhc_model_df = imhc_model_df[mask]


def get_cdr3_aa_prop_length_fraction( cdr3, props ):
    if not cdr3:
        return 0.0 # probably should return mean(props) actually. But that would require refitting the score coefs
        # and a CDR3 of length 8 or less is actually pretty short...
    else:
        return sum( props[x] for x in cdr3 ) / float(len(cdr3))


def get_feature(tcr, ftag, aa_props_df):
    ''' taken from make_table2.py but note that tcr cdr3s here are not trimmed to [4:-4] yet
    '''
    trim_start, trim_stop = 4,-4
    if ftag.endswith('_AB'):
        cdr3s = [ tcr[0][2][trim_start:trim_stop], tcr[1][2][trim_start:trim_stop] ]
    elif ftag.endswith('_A'):
        cdr3s = [ tcr[0][2][trim_start:trim_stop] ]
    else:
        assert ftag.endswith('_B')
        cdr3s = [ tcr[1][2][trim_start:trim_stop] ]

    ftag = '_'.join( ftag.split('_')[:-1] )
    cdr3len = sum( len(x) for x in cdr3s )
    frac_norm = 1.0 / cdr3len if cdr3len else 1.0

    if ftag == 'len':
        return cdr3len
    elif ftag[0] in amino_acids and ftag[1:] == 'frac':
        return sum( x.count(ftag[0]) for x in cdr3s ) * frac_norm
    elif ftag == 'arofrac':
        return sum( x.count(y) for x in cdr3s for y in 'FYWH' ) * frac_norm
    else:
        return sum( get_cdr3_aa_prop_length_fraction(x, aa_props_df[ftag] ) for x in cdr3s)
    return None

def make_imhc_score_table_column(tcrs, aa_props_df):
    scores = np.zeros((len(tcrs),))
    for row in imhc_model_df.itertuples():
        col = np.array( [ get_feature(x, row.Index, aa_props_df) for x in tcrs] )
        col = row.coef * (col - row.mean)/np.sqrt(row.var)
        scores += col
        intercept = row.intercept # all the same
    scores += intercept
    return scores

def get_imhc_raw_score_terms_and_coefs(tcrs, aa_props_df):
    ''' Just for diagnostics/visualization
    '''
    cols = []
    features = []
    coefs = []
    for row in imhc_model_df.itertuples():
        col = np.array( [ get_feature(x, row.Index, aa_props_df) for x in tcrs] )
        cols.append(  (col - row.mean)/np.sqrt(row.var) )
        coefs.append(row.coef)
        features.append(row.Index)

    return np.array(cols).transpose(), features, coefs




