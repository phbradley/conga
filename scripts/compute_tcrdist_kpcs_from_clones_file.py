import argparse
import sys
import os
#from os import system
from os.path import exists
import numpy as np
from sklearn.decomposition import KernelPCA
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # so we can import conga
import conga


def compute_tcrdist_kpcs_from_clones_file(
        clones_file,
        organism,
        n_components_in=50
):
    assert organism in ['mouse', 'human']

    python2 = conga.util.PYTHON2_EXE
    script = conga.util.TCRDIST_REPO+'compute_distances.py'
    if not exists(script):
        print('cant find the tcr-dist script:', script)
        print('please install tcr-dist repo from github and set the location in conga/conga/util.py at the top')
        if not exists(python2):
            print('also the python2 executable, too')
        sys.exit(1)

    cmd = '{} {} --organism {} --dist_chains AB --clones_file {}'\
        .format( python2, script, organism, clones_file )
    conga.util.run_command(cmd) # I think this just calls system right now

    distfile = clones_file[:-4]+'_AB.dist'
    if not exists(distfile):
        print('computing the distances with tcr-dist failed!')
        print('the command was:', cmd)
        sys.exit(1)

    ## read distances
    counter=0
    D=[]
    ids = []
    for line in open( distfile,'r'):
        if counter%1000==0:
            print( 'reading', distfile, counter)
        counter+=1
        l = line.split()
        dists = [ float(x) for x in l[1:] ]
        D.append( dists )
        ids.append(l[0])

    D = np.matrix(D)

    n_components = min( n_components_in, D.shape[0] )

    print( 'running KernelPCA', D.shape)

    pca = KernelPCA(kernel='precomputed', n_components=n_components)

    gram = 1 - ( D / D.max() )
    xy = pca.fit_transform(gram)

    #show the eigenvalues
    for ii in range(n_components):
        print( 'eigenvalue: {:3d} {:.3f}'.format( ii, pca.lambdas_[ii]))

    outfile = '{}_{}_kpcs'.format(distfile, n_components_in)
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

