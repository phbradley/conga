#from . import tcr_scores
#from .tcr_scores import aa_props_df
import numpy as np
import pandas as pd

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

#aa_props_feature_set = frozenset(tcr_scores.aa_props_df.columns)


mhci_model_file ='/home/pbradley/gitrepos/conga/conga/data/logreg_hobit_donor1_v4.tsv_0_model.tsv'
mhci_model_df = pd.read_csv(mhci_model_file, sep='\t')
mhci_model_df.set_index('feature', inplace=True)

# restrict to features with nonzero weight
mask = (np.abs(mhci_model_df['coef'])>1e-6)
mhci_model_df = mhci_model_df[mask]


def get_cdr3_aa_prop_length_fraction( cdr3, props ):
    if not cdr3:
        return 0.0
    else:
        return sum( props[x] for x in cdr3 ) / float(len(cdr3))


def get_feature(tcr, ftag, aa_props_df):
    ''' taken from make_table2.py but note that tcr cdr3s here are not trimmed to [4:-4] yet
    '''
    trim_start, trim_stop = 4,-4
    if '_AB' in ftag:
        cdr3s = [ tcr[0][2][trim_start:trim_stop], tcr[1][2][trim_start:trim_stop] ]
    elif '_A' in ftag:
        cdr3s = [ tcr[0][2][trim_start:trim_stop] ]
    else:
        assert '_B' in ftag
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

def make_mhci_score_table_column(tcrs, aa_props_df):
    scores = np.zeros((len(tcrs),))
    for row in mhci_model_df.itertuples():
        col = np.array( [ get_feature(x, row.Index, aa_props_df) for x in tcrs] )
        col = row.coef * (col - row.mean)/np.sqrt(row.var)
        scores += col
        intercept = row.intercept # all the same
    scores += intercept
    return scores

def get_mhci_raw_score_terms_and_coefs(tcrs, aa_props_df):
    ''' Just for diagnostics/visualization
    '''
    cols = []
    features = []
    coefs = []
    for row in mhci_model_df.itertuples():
        col = np.array( [ get_feature(x, row.Index, aa_props_df) for x in tcrs] )
        cols.append(  (col - row.mean)/np.sqrt(row.var) )
        coefs.append(row.coef)
        features.append(row.Index)

    return np.array(cols).transpose(), features, coefs




