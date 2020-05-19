import numpy as np
import sys
from os import system
from scipy.sparse import issparse
from collections import Counter
#import pandas as pd
from . import preprocess as pp

PYTHON2_EXE = '/home/pbradley/anaconda2/bin/python2.7'
TCRDIST_REPO = '/home/pbradley/gitrepos/tcr-dist/'

def get_mean_var_scanpy(X):
    ''' this is based on---ie stolen from---scanpy's _get_mean_var function
    wanted to be explicit about what our variance convention is
    '''
    #
    mean = X.mean(axis=0)
    if issparse(X):
        mean_sq = X.multiply(X).mean(axis=0)
        mean = mean.A1
        mean_sq = mean_sq.A1
    else:
        mean_sq = np.multiply(X, X).mean(axis=0)
    # do we want unbiased estimator?
    #var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    var = (mean_sq - mean**2)
    return mean, var


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


def get_vfam(vgene):
    assert vgene.startswith('TR') and vgene[3]=='V'
    pos = 4
    while pos<len(vgene) and vgene[pos].isdigit():
        pos += 1
    vno = int(vgene[4:pos])
    return '{}V{:d}'.format(vgene[2], vno)



def setup_tcr_cluster_names(adata):

    organism = adata.uns['organism']

    #clusters_gex = adata.obs['clusters_gex']
    clusters_tcr = adata.obs['clusters_tcr']

    tcrs = pp.retrieve_tcrs_from_adata(adata)

    num_clusters = np.max(clusters_tcr)+1
    names = []
    for c in range(num_clusters):
        cluster_size = np.sum(clusters_tcr==c)
        ctcrs = [ x for x,y in zip(tcrs,clusters_tcr) if y==c]
        counts = Counter( [ get_vfam(x[0][0]) for x in ctcrs]) +Counter([get_vfam(x[1][0]) for x in ctcrs])
        print( c, cluster_size, counts.most_common(3))
        top_vfam, top_count = counts.most_common(1)[0]
        eps=1e-3
        if top_count+eps >= 0.75*cluster_size:
            names.append('{}_{}'.format(c, top_vfam))
        elif top_count+eps >= 0.5*cluster_size:
            names.append('{}_{}'.format(c, top_vfam.lower()))
        else:
            if organism=='human': # special hack
                c2 = counts['AV14']+counts['AV38']
                c3 = counts['AV14']+counts['AV38']+counts['AV19']
                if c2 >= 0.75*cluster_size:
                    names.append('{}_AV14+'.format(c, top_vfam))
                elif c3 >= 0.75*cluster_size:
                    names.append('{}_AV14++'.format(c, top_vfam))
                elif c2 >= 0.5*cluster_size:
                    names.append('{}_av14+'.format(c, top_vfam))
                elif c3 >= 0.5*cluster_size:
                    names.append('{}_av14++'.format(c, top_vfam))
                else:
                    names.append(str(c))
            else:
                names.append(str(c))
    print('setup_tcr_cluster_names:', names)
    adata.uns['clusters_tcr_names'] = names
