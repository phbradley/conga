from os.path import exists
from pathlib import Path
import numpy as np
from . import util

# comparison of the correlations between CD4/CD8 expression levels and
# either the old cd8 score or this new one on two 10x datasets
# Looks like this new one is a better predictor (which is should be!)
#
# based on this I am switching this to be the default
#
#                   R    Pearson P  Spearman  Kendall  dataset
# old corr: CD4   -0.171  1.10e-19  5.25e-20  1.07e-19 hs_pbmc3
# new corr: CD4   -0.184  1.60e-22  1.03e-23  2.27e-23 hs_pbmc3
# old corr: CD8A   0.360  3.56e-86  1.42e-91  4.98e-85 hs_pbmc3
# new corr: CD8A   0.384  1.57e-98 7.91e-109 6.17e-100 hs_pbmc3
# old corr: CD8B   0.312  5.94e-64  3.23e-74  1.25e-69 hs_pbmc3
# new corr: CD8B   0.329  3.76e-71  1.55e-91  1.28e-84 hs_pbmc3
# old corr: CD4   -0.136  8.38e-08  8.22e-08  1.02e-07 hs_pbmc
# new corr: CD4   -0.154  1.37e-09  4.40e-09  5.52e-09 hs_pbmc
# old corr: CD8A   0.348  3.83e-45  2.75e-46  2.76e-43 hs_pbmc
# new corr: CD8A   0.400  2.49e-60  1.46e-63  1.40e-57 hs_pbmc
# old corr: CD8B   0.258  8.65e-25  2.64e-30  1.75e-28 hs_pbmc
# new corr: CD8B   0.303  5.09e-34  5.89e-45  1.38e-41 hs_pbmc

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# read the model parameters
all_models = {}

for ab in 'AB':
    model_params = {
        'window_size':None,
        'min_lenbin':None,
        'max_lenbin':None,
        'NV':None,
        'NJ':None,
        'NL':None,
        'NC':None,
        }
    model_file = Path.joinpath( Path(util.path_to_data),f'cd8_logreg_params_{ab}.txt')
    tags, weights = [], []
    #print(f'loading cd8 model for chain {ab} {model_file}')
    with open(model_file, 'r') as data:
        for line in data:
            l = line.split()
            tag, value = l
            if tag in model_params:
                model_params[tag] = int(value)
            else:
                tags.append(tag)
                weights.append(float(value))
    NV, NJ, NL, NC = [model_params[x] for x in 'NV NJ NL NC'.split()]
    NTOT = NV + NJ + NL + NC
    assert tags[-1] == 'BIAS'
    bias = float(weights[-1])
    assert len(weights) == NTOT+1
    weights = np.array(weights)
    vgenes = tags[:NV-1] # last is UNK gene
    jgenes = tags[NV:NV+NJ-1]
    assert all(x.startswith(f'TR{ab}V') for x in vgenes)
    assert all(x.startswith(f'TR{ab}J') for x in jgenes)
    vgene_indexer = {x:i for i,x in enumerate(vgenes)}
    jgene_indexer = {x:i for i,x in enumerate(jgenes)}

    model_params['weights'] = weights
    model_params['vgene_indexer'] = vgene_indexer
    model_params['jgene_indexer'] = jgene_indexer

    all_models[ab] = model_params


def get_lenbin(L, min_lenbin, max_lenbin):
    if L <= min_lenbin:
        return min_lenbin
    elif L >= max_lenbin:
        return max_lenbin
    else:
        return L

def get_allele(g):
    assert g.count('*')==1
    return g[:g.index('*')]

def encode_single_chain_tcr(vgene, jgene, cdr3, model_params):
    ''' vgene and jgene still include '*01' (e.g.)
    '''
    window_size = model_params['window_size']
    min_lenbin = model_params['min_lenbin']
    max_lenbin = model_params['max_lenbin']
    vgene_indexer = model_params['vgene_indexer']
    jgene_indexer = model_params['jgene_indexer']
    NV = len(vgene_indexer)+1
    NJ = len(jgene_indexer)+1
    NL = max_lenbin-min_lenbin+1
    NC = 20 * (2*window_size+1)
    NTOT = NV + NJ + NL + NC

    # 1-hot encode the V/J genes and length
    x = np.zeros((NTOT+1,))
    iv = vgene_indexer.get(get_allele(vgene), NV-1)
    x[iv] = 1.0

    ij = jgene_indexer.get(get_allele(jgene), NJ-1)
    x[NV+ij] = 1.0

    il = get_lenbin(len(cdr3), min_lenbin, max_lenbin) - min_lenbin
    x[NV+NJ+il] = 1.0

    # 1-hot encode aas in the nterminal and cterminal cdr3 windows
    # k-hot encode the middle of the tcr
    nterm = min(window_size, len(cdr3)//2)
    cterm = min(window_size, len(cdr3)-nterm)
    #middle = len(cdr3)-nterm-cterm
    for i in range(window_size):
        if i<nterm:
            istart = NV+NJ+NL+20*i
            x[istart+amino_acids.index(cdr3[i])] += 1.0
        if i<cterm:
            istart = NV+NJ+NL+20*window_size+20*i
            x[istart+amino_acids.index(cdr3[-1-i])] += 1.0
    for aa in cdr3[nterm:-cterm]:
        assert len(cdr3)>2*window_size # sanity check
        istart = NV+NJ+NL+40*window_size
        x[istart+amino_acids.index(aa)] += 1.0

    x[NTOT] = 1. # for the bias term
    return x


def make_cd8_score_table_column(tcrs, use_sigmoid=False):
    ''' We fit a logistic regression model on single-chain bulk TCR data
    from sorted CD4/CD8 subsets. Here we just average the alpha-chain and beta-chain
    scores.
    '''
    global all_models
    score_totals = np.zeros((len(tcrs),))
    for iab, ab in enumerate('AB'):
        model = all_models[ab]
        tcr_vectors = [ encode_single_chain_tcr(x[iab][0], x[iab][1], x[iab][2], model)
                        for x in tcrs ]
        tcr_vectors = np.vstack(tcr_vectors)
        weights = model['weights']
        assert tcr_vectors.shape == (len(tcrs), weights.shape[0])
        ab_scores = np.dot( tcr_vectors, weights)
        if use_sigmoid:
            ab_scores = 1.0/(1.0+np.exp(-1*ab_scores))
        score_totals += 0.5 * ab_scores
    return score_totals





