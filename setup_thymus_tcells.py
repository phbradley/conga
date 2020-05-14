import conga
import conga.preprocess as pp
import conga.correlations as cc
import conga.plotting as pl
import scanpy as sc
import scanpy.neighbors
from sklearn.metrics import pairwise_distances
import numpy as np
#from sys import exit
from collections import Counter
from os.path import exists
from anndata import AnnData

# this is silly but for some reason adata.uns can't hold a list of strings with alphas and betas,
# or at least we can't save such an adata to a file
#
def safestring(s):
    news = []
    for a in s:
        if ord(a)<128:
            news.append(a)
        elif ord(a)==945:
            news.append('a')
        elif ord(a)==946:
            news.append('b')
        elif ord(a)==947:
            news.append('g')
        elif ord(a)==948:
            news.append('d')
        # elif ord(a)==206:
        #     pass
        else:
            news.append('?')
    return ''.join(news)

# we are using the X_umap_3d array in lieu of GEX PCA
X_gex_tag = 'X_umap_3d'

gex_file = '/home/pbradley/csdat/tcr10x/other_vdj/thymus/all_tcell_clones.h5ad'
clones_file = '/home/pbradley/csdat/tcr10x/other_vdj/thymus/all_tcell_clones.tsv'

adata = pp.read_dataset(gex_file, 'h5ad', clones_file )
assert not adata.isview
assert 'X_pca_tcr' in adata.obsm_keys() # tcr-dist kPCA info
assert 'cdr3a' in adata.obs # tcr sequence (VDJ) info (plus other obs keys)
assert X_gex_tag in adata.obsm_keys()
print(adata)

# set organism
adata.uns['organism'] = 'human'

# the X matrix in adata is already scaled and log1p'ed. But it was scaled to 100,000 counts per cell, and
# we're used to looking at data scaled to 10,000 counts per cell. So we renormalize here:
X2 = np.expm1( adata.X )
celltotals = X2.sum(axis=1).A1
print('celltotals:', celltotals.shape, np.min(celltotals), np.max(celltotals))
assert np.min(celltotals) > 1e5-1 and np.max(celltotals) < 1e5+1
X2 *= 0.1

new_celltotals = X2.sum(axis=1).A1
assert min(new_celltotals) > 1e4-1 and max(new_celltotals) < 1e4+1

np.log1p( X2.data, out = X2.data )

adata_new = AnnData( X = X2, obs = adata.obs, var = adata.var )
adata.raw = adata_new
pp.set_raw_matrix_is_logged_to_true(adata)

# we are going to skip the reduce_to_single_cell_per_clone(...) step, so we need to do the same stuff:
# clone_sizes, X_igex, X_igex_genes
X_igex, good_genes = pp.setup_X_igex(adata)
adata.obsm['X_igex'] = X_igex
adata.uns['X_igex_genes'] = good_genes
num_clones = adata.shape[0] # not quite right
adata.obs['clone_sizes'] = [1]*num_clones # this is not quite right, actually. There are a few expanded clones.

# we are also skipping filter_and_scale
# and cluster_and_tsne_and_umap, so we need:
# X_pca_gex, X_gex_2d, X_tcr_2d, clusters_gex, and clusters_tcr
#
adata.obsm['X_pca_gex'] = adata.obsm[X_gex_tag].copy() # probably X_umap_3d
adata.obsm['X_gex_2d'] = adata.obsm['X_umap']

# make fake clusters_gex using the cell types array
celltypes = adata.obs['cell types']
celltypes_sorted = sorted(set(celltypes))
print('celltypes_sorted:', ' '.join( '{}: {}'.format(x,y) for x,y in enumerate(celltypes_sorted)))

celltype2cluster = { y:x for x,y in enumerate(celltypes_sorted) }

clusters_gex = np.array( [ celltype2cluster[x] for x in celltypes ] )
adata.obs['clusters_gex'] = clusters_gex
adata.uns['clusters_gex_names'] = [ safestring(x) for x in celltypes_sorted ] ## store these for annotation of logos

# compute X_tcr_2d and clusters_tcr
adata.obsm['X_pca'] = adata.obsm['X_pca_tcr']

print('tcr nbrs')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40,num_clones))

print('tcr umap')
sc.tl.umap(adata)
adata.obsm['X_tcr_2d'] = adata.obsm['X_umap']

print('tcr louvain')
sc.tl.louvain(adata, key_added='louvain', resolution=1.0)
adata.obs['clusters_tcr'] = np.copy( adata.obs['louvain'] ).astype(int)

adata.write_h5ad('thymus_tcells_checkpoint1.h5ad')

# we could go ahead and do the conga part now too:
print('Now you can run:\n')

cmd = '/home/pbradley/anaconda2/envs/scanpy_new/bin/python run_conga.py --organism human --from_checkpoint1 thymus_tcells_checkpoint1.h5ad --checkpoint --nbr_fracs 0.01 0.1 --outfile_prefix thymus_tcells_new'
print(cmd)

