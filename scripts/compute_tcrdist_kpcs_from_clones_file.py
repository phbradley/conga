import argparse
import sys
import os
#from os import system
from os.path import exists
import numpy as np
from sklearn.decomposition import KernelPCA

sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # so we can import conga
import conga
from conga.tcrdist.tcr_distances import TcrDistCalculator
import pandas as pd


def compute_tcrdist_kpcs_from_clones_file(
        clones_file,
        organism,
        n_components_in=50
):
    assert organism in ['mouse', 'human']

    tcrdist_calculator = TcrDistCalculator(organism)

    df = pd.read_csv(clones_file, sep='\t')

    # in conga we usually also have cdr3_nucseq but we don't need it for tcrdist
    tcrs = [ ( ( l.va_gene, l.ja_gene, l.cdr3a ), ( l.vb_gene, l.jb_gene, l.cdr3b ) ) for l in df.itertuples() ]
    ids = [ l.clone_id for l in df.itertuples() ]


    ## read distances
    print(f'compute tcrdist distance matrix for {len(tcrs)} clonotypes')
    D= np.array( [ tcrdist_calculator(x,y) for x in tcrs for y in tcrs ] ).reshape( (len(tcrs), len(tcrs)) )

    #D = np.matrix(D) # I dunno if this is important, I think it came from some example

    n_components = min( n_components_in, D.shape[0] )

    print( 'running KernelPCA', D.shape)

    pca = KernelPCA(kernel='precomputed', n_components=n_components)

    gram = 1 - ( D / D.max() )
    xy = pca.fit_transform(gram)

    #show the eigenvalues
    for ii in range(n_components):
        print( 'eigenvalue: {:3d} {:.3f}'.format( ii, pca.lambdas_[ii]))

    # this is the kpcs_file that conga.preprocess.read_dataset is expecting:
    #kpcs_file = clones_file[:-4]+'_AB.dist_50_kpcs'
    outfile = '{}_AB.dist_{}_kpcs'.format(clones_file[:-4], n_components_in)
    print( 'making:', outfile)
    out = open(outfile,'w')

    for ii in range(D.shape[0]):
        out.write('pc_comps: {} {}\n'\
                  .format( ids[ii], ' '.join( '{:.6f}'.format(xy[ii,j]) for j in range(n_components) ) ) )
    out.close()
    return



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--clones_file')
    parser.add_argument('--organism', choices=['mouse', 'human'], required=True)
    parser.add_argument('--n_components', type=int, default=50)

    args = parser.parse_args()
    assert exists(args.clones_file)

    compute_tcrdist_kpcs_from_clones_file(args.clones_file, args.organism, args.n_components)

