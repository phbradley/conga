import argparse
from os.path import exists

parser = argparse.ArgumentParser()

parser.add_argument('--clones_file')
parser.add_argument('--organism', choices=['mouse', 'human'], required=True)
parser.add_argument('--n_components', type=int, default=50)

args = parser.parse_args()
assert exists(args.clones_file)

## more imports-- slow, so after initial arg parse #######################33
from os import system
import sys
import numpy as np
from sklearn.decomposition import KernelPCA
import conga
# import conga.preprocess as pp
# import conga.correlations as cc
# import conga.plotting as pl
# import scanpy as sc
# import scanpy.neighbors
# from sklearn.metrics import pairwise_distances
# import numpy as np
# import pandas as pd
# from collections import Counter


# first compute the distance matrix
python2 = conga.util.PYTHON2_EXE
script = conga.util.TCRDIST_REPO+'compute_distances.py'
if not exists(script):
    print('cant find the tcr-dist script:', script)
    print('please install tcr-dist repo from github and set the location in conga')
    sys.exit(1)

cmd = '{} {} --organism {} --dist_chains AB --clones_file {}'\
    .format( python2, script, args.organism, args.clones_file )
print('run:', cmd)
sys.stdout.flush()
#system(cmd)

distfile = args.clones_file[:-4]+'_AB.dist'
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

n_components = min( args.n_components, D.shape[0] )

print( 'running KernelPCA', D.shape)

pca = KernelPCA(kernel='precomputed', n_components=n_components)

gram = 1 - ( D / D.max() )
xy = pca.fit_transform(gram)

#show the eigenvalues
for ii in range(n_components):
    print( 'eigenvalue: {:3d} {:.3f}'.format( ii, pca.lambdas_[ii]))

outfile = '{}_{}_kpcs'.format(distfile, args.n_components)
print( 'making:', outfile)
out = open(outfile,'w')

for ii in range(D.shape[0]):
    out.write('pc_comps: {} {}\n'\
              .format( ids[ii], ' '.join( '{:.6f}'.format(xy[ii,j]) for j in range(n_components) ) ) )
out.close()

