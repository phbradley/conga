import argparse
import conga
#import conga.tcr_scores as tcr_scores

parser = argparse.ArgumentParser()

#type is str by default
parser.add_argument('--gex_data')
parser.add_argument('--gex_data_type', choices=['h5ad', '10x_mtx', '10x_h5'])
parser.add_argument('--clones_file')
parser.add_argument('--organism', choices=['mouse', 'human'])
parser.add_argument('--nbr_fracs', type=float, nargs='*', default=[0.01,0.1] )
parser.add_argument('--min_cluster_size', type=int, default=5)
parser.add_argument('--min_cluster_size_repsize', type=int, default=5000)
parser.add_argument('--outfile_prefix', required=True)
parser.add_argument('--checkpoint', action='store_true')
parser.add_argument('--from_checkpoint1')
parser.add_argument('--make_unfiltered_logos', action='store_true')
parser.add_argument('--make_avggood_logos', action='store_true')
parser.add_argument('--write_proj_info', action='store_true')
parser.add_argument('--filter_ribo_norm_low_cells', action='store_true')
parser.add_argument('--calc_clone_pmhc_pvals', action='store_true')
parser.add_argument('--find_nbrhood_overlaps', action='store_true')
parser.add_argument('--find_pmhc_nbrhood_overlaps', action='store_true') # only if pmhc info is present
parser.add_argument('--find_tcr_nbrhood_genes', action='store_true')
parser.add_argument('--find_tcr_cluster_genes', action='store_true')
parser.add_argument('--find_tcr_segment_genes', action='store_true')
parser.add_argument('--find_gex_nbrhood_scores', action='store_true')
parser.add_argument('--find_gex_cluster_scores', action='store_true')
parser.add_argument('--tenx_agbt', action='store_true')
parser.add_argument('--gex_nbrhood_tcr_score_names', type=str, nargs='*',
                    default=conga.tcr_scores.all_tcr_scorenames )

args = parser.parse_args()

if args.find_pmhc_nbrhood_overlaps or args.calc_clone_pmhc_pvals:
    # we need pmhc info for these analyses; right now that's restricted to the 10x AGBT dataset format
    assert args.tenx_agbt


import conga.preprocess as pp
import conga
import conga.correlations as cc
import conga.plotting as pl
import scanpy as sc
import scanpy.neighbors
from sklearn.metrics import pairwise_distances
import numpy as np
#from sys import exit
from collections import Counter
from os.path import exists

logfile = args.outfile_prefix+'_log.txt'
outlog = open(logfile, 'w')

if args.from_checkpoint1 is None:

    assert exists(args.gex_data)
    assert exists(args.clones_file)

    ## load the dataset
    adata = pp.read_dataset(args.gex_data, args.gex_data_type, args.clones_file )
    assert args.organism
    adata.uns['organism'] = args.organism
    assert 'organism' in adata.uns_keys()

    if args.tenx_agbt:
        adata.uns['pmhc_var_names'] = conga.pmhc_scoring.get_tenx_agbt_pmhc_var_names(adata)
        print('pmhc_var_names:', adata.uns['pmhc_var_names'])

    assert not adata.isview
    assert 'X_pca_tcr' in adata.obsm_keys() # tcr-dist kPCA info
    assert 'cdr3a' in adata.obs # tcr sequence (VDJ) info (plus other obs keys)

    print(adata)

    adata = pp.filter_and_scale( adata )

    if args.calc_clone_pmhc_pvals: # do this before condensing to a single clone per cell
        results_df = conga.pmhc_scoring.calc_clone_pmhc_pvals(adata)
        tsvfile = args.outfile_prefix+'_clone_pvals.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, sep='\t')

    if args.filter_ribo_norm_low_cells:
        adata = pp.filter_cells_by_ribo_norm( adata )

    adata = pp.reduce_to_single_cell_per_clone( adata )
    assert 'X_igex' in adata.obsm_keys()

    adata = pp.cluster_and_tsne_and_umap( adata )

    if args.checkpoint:
        adata.write_h5ad(args.outfile_prefix+'_checkpoint1.h5ad')
else:
    adata = sc.read_h5ad(args.from_checkpoint1)
    print('recover from checkpoint:', adata )

    if 'organism' not in adata.uns_keys():
        assert args.organism
        adata.uns['organism'] = args.organism


if args.write_proj_info:
    outfile = args.outfile_prefix+'_2d_proj_info.txt'
    pp.write_proj_info( adata, outfile )

clusters_gex = adata.obs['clusters_gex']
clusters_tcr = adata.obs['clusters_tcr']

tcrs = pp.retrieve_tcrs_from_adata(adata)
barcode2tcr = { x:y for x,y in zip( adata.obs_names, tcrs )}
num_clones = len(tcrs)

cc.compute_cluster_interactions( clusters_gex, clusters_tcr, adata.obs_names, barcode2tcr, outlog )

atcrs = sorted( set( x[0] for x in tcrs ) )
btcrs = sorted( set( x[1] for x in tcrs ) )
atcr2agroup = dict( (y,x) for x,y in enumerate(atcrs))
btcr2bgroup = dict( (y,x) for x,y in enumerate(btcrs))
agroups = np.array( [ atcr2agroup[x[0]] for x in tcrs] )
bgroups = np.array( [ btcr2bgroup[x[1]] for x in tcrs] )

print('compute D_gex')
D_gex = pairwise_distances( adata.obsm['X_pca_gex'], metric='euclidean' )

print('compute D_tcr')
D_tcr = pairwise_distances( adata.obsm['X_pca_tcr'], metric='euclidean' )

for ii,a in enumerate(agroups):
    D_gex[ii, (agroups==a) ] = 1e3
    D_tcr[ii, (agroups==a) ] = 1e3
for ii,b in enumerate(bgroups):
    D_gex[ii, (bgroups==b) ] = 1e3
    D_tcr[ii, (bgroups==b) ] = 1e3

bad_conga_score = -1*np.log10(num_clones)
conga_scores = np.full( (num_clones,3), bad_conga_score)

if args.find_nbrhood_overlaps or args.find_pmhc_nbrhood_overlaps:
    for nbrfrac in args.nbr_fracs:

        num_neighbors = max(1, int(nbrfrac*num_clones))
        nbrs_gex = np.argpartition( D_gex, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
        nbrs_tcr = np.argpartition( D_tcr, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
        assert nbrs_tcr.shape == (num_clones, num_neighbors) and nbrs_gex.shape == nbrs_tcr.shape

        if args.find_nbrhood_overlaps:

            pval_threshold = 1. # since they are being scaled by num_clones
            print('find_neighbor_neighbor_interactions:')
            results = cc.find_neighbor_neighbor_interactions( nbrs_gex, nbrs_tcr, agroups, bgroups, pval_threshold )

            for pval, ii, overlap in results:
                conga_scores[ii,0] = max(conga_scores[ii,0], -1*np.log10(pval) )

            print('find_neighbor_cluster_interactions:')
            results = cc.find_neighbor_cluster_interactions( nbrs_tcr, clusters_gex, agroups, bgroups, pval_threshold)
            for pval, ii, overlap in results:
                conga_scores[ii,1] = max(conga_scores[ii,1], -1*np.log10(pval) )

            print('find_neighbor_cluster_interactions:')
            results = cc.find_neighbor_cluster_interactions( nbrs_gex, clusters_tcr, agroups, bgroups, pval_threshold)
            for pval, ii, overlap in results:
                conga_scores[ii,2] = max(conga_scores[ii,2], -1*np.log10(pval) )

        if args.find_pmhc_nbrhood_overlaps:
            for tag, nbrs in [['gex', nbrs_gex], ['tcr', nbrs_tcr]]:
                conga.pmhc_scoring.compute_pmhc_versus_nbrs(
                    adata, nbrs, agroups, bgroups, '{}_{}'.format(tag, nbrfrac))

else:
    # compute nbrs for other analyses, using largest of nbr_fracs
    nbrfrac = max(args.nbr_fracs)
    num_neighbors = max(1, int(nbrfrac* num_clones))
    nbrs_gex = np.argpartition( D_gex, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
    nbrs_tcr = np.argpartition( D_tcr, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
    assert nbrs_tcr.shape == (num_clones, num_neighbors) and nbrs_gex.shape == nbrs_tcr.shape


adata.obsm['conga_scores'] = conga_scores

max_conga_score = np.max( conga_scores, axis=1 )
good_mask = (max_conga_score >= 0.)

outlog.write('num_good: {}\n'.format(np.sum(good_mask)))

adata.obs['good_score_mask'] = good_mask

if num_clones <= args.min_cluster_size_repsize:
    min_cluster_size = args.min_cluster_size
else:
    min_cluster_size = int( 0.5 + args.min_cluster_size * float(num_clones)/args.min_cluster_size_repsize )


if args.checkpoint:
    adata.write_h5ad(args.outfile_prefix+'_checkpoint2.h5ad')

if args.find_tcr_nbrhood_genes:
    pval_threshold = 0.5
    cc.tcr_nbrhood_rank_genes_fast( adata, nbrs_tcr, pval_threshold)

    ####cc.tcr_nbrhood_rank_genes( adata, nbrs_tcr, pval_threshold, rank_method=args.tcr_nbrhood_rank_genes_method)

if args.find_tcr_cluster_genes:
    # make some fake nbrs
    fake_nbrs_tcr = []
    seen = set()
    for ii, cl in enumerate(clusters_tcr):
        if cl in seen:
            fake_nbrs_tcr.append([])
        else:
            seen.add(cl)
            fake_nbrs_tcr.append(np.nonzero( clusters_tcr==cl )[0])
    pval_threshold = 0.5

    cc.tcr_nbrhood_rank_genes_fast( adata, fake_nbrs_tcr, pval_threshold, prefix_tag='clust')

if args.find_tcr_segment_genes:
    for iab,ab in enumerate('AB'):
        for iseg,seg in enumerate('VJ'):
            genes = [ x[iab][iseg] for x in tcrs ]
            genes = np.array([ x[:x.index('*')] for x in genes ])

            # make some fake nbrs
            fake_nbrs_tcr = []
            clone_display_names = []
            seen = set()
            for ii, g in enumerate(genes):
                if g in seen:
                    fake_nbrs_tcr.append([])
                    clone_display_names.append('')
                else:
                    seen.add(g)
                    fake_nbrs_tcr.append(np.nonzero( genes==g )[0] ) # this will include ii but don't think thats a prob
                    clone_display_names.append(g)

            pval_threshold = .5

            cc.tcr_nbrhood_rank_genes_fast( adata, fake_nbrs_tcr, pval_threshold, prefix_tag=seg+ab,
                                            clone_display_names=clone_display_names )



if args.find_gex_nbrhood_scores:
    pval_threshold = 0.5
    cc.gex_nbrhood_rank_tcr_scores( adata, nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold )


if args.find_gex_cluster_scores:
    # make some fake nbrs
    fake_nbrs_gex = []
    seen = set()
    for ii, cl in enumerate(clusters_gex):
        if cl in seen:
            fake_nbrs_gex.append([])
        else:
            seen.add(cl)
            fake_nbrs_gex.append(np.nonzero( clusters_gex==cl )[0])

    pval_threshold = 0.5
    cc.gex_nbrhood_rank_tcr_scores( adata, fake_nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold,
                                    prefix_tag = 'clust' )



clp_counts = Counter( (x,y) for x,y,m in zip(clusters_gex, clusters_tcr, good_mask) if m )
num_good_cluster_pairs = sum( 1 for x,y in clp_counts.items() if y>=min_cluster_size )

print('num_good_cluster_pairs:', num_good_cluster_pairs)

if num_good_cluster_pairs:
    # run rank_genes on most common clps
    rank_genes_uns_tag = 'rank_genes_good_cluster_pairs'
    cc.run_rank_genes_on_good_cluster_pairs( adata, good_mask, clusters_gex, clusters_tcr, min_count=min_cluster_size,
                                             key_added= rank_genes_uns_tag)


    rank_genes_filt_uns_tag = 'rank_genes_good_cluster_pairs_filtered'

    sc.tools.filter_rank_genes_groups(adata, key=rank_genes_uns_tag, groupby='test',
                                      key_added=rank_genes_filt_uns_tag,
                                      min_in_group_fraction=0.25, min_fold_change=2,
                                      max_out_group_fraction=1.1) # 1.1 means don't filter for this (dflt was 0.5)

    # these use nbrs from last nbrfrac
    gex_cluster_names = adata.uns['clusters_gex_names'] if 'clusters_gex_names' in adata.uns else None

    pl.make_logo_plots( adata, nbrs_gex, nbrs_tcr, min_cluster_size, args.outfile_prefix+'_good_logos_rgfilt.png',
                        gex_cluster_names = gex_cluster_names,
                        rank_genes_uns_tag = rank_genes_filt_uns_tag )

    if args.make_unfiltered_logos:
        pl.make_logo_plots( adata, nbrs_gex, nbrs_tcr, min_cluster_size, args.outfile_prefix+'_good_logos.png',
                            rank_genes_uns_tag = rank_genes_uns_tag )

# try clustering using an averaged distance matrix
if np.sum(good_mask) >= min_cluster_size:
    for subset_tag in ['good','full']:
        if subset_tag=='good':
            if not args.make_avggood_logos:
                continue
            if np.sum(good_mask) < min_cluster_size:
                continue
            adata_sub = adata[good_mask,:].copy()
            sub2full = dict( enumerate( np.nonzero(good_mask)[0] ) )

        else:
            adata_sub = adata.copy()
            sub2full = {i:i for i in range(adata.shape[0]) }


        print('masked:', adata_sub.shape, adata_sub.isview)

        num_neighbors = max(1, min(10, adata_sub.shape[0]//2) )

        D_gex = pairwise_distances( adata_sub.obsm['X_pca_gex'], metric='euclidean' )
        print('D_gex mean median max std', np.mean(D_gex), np.median(D_gex), np.max(D_gex), np.std(D_gex))

        D_tcr = pairwise_distances( adata_sub.obsm['X_pca_tcr'], metric='euclidean' )
        print('D_tcr mean median max std', np.mean(D_tcr), np.median(D_tcr), np.max(D_tcr), np.std(D_tcr))

        D_avg = (( D_gex/np.median(D_gex) )**2 + ( D_tcr/np.median(D_tcr) )**2)**0.5

        # this is borrowed from scanpy/neighbors/__init__.py
        #_distances = pairwise_distances(X, metric=metric, **metric_kwds)
        knn_indices, knn_distances = scanpy.neighbors.get_indices_distances_from_dense_matrix(
            D_avg, num_neighbors)

        distances, connectivities = scanpy.neighbors.compute_connectivities_umap(
            knn_indices, knn_distances, adata_sub.shape[0], num_neighbors)

        adata_sub.uns['neighbors'] = {}
        adata_sub.uns['neighbors']['params'] = {'n_neighbors': num_neighbors, 'method': 'umap'}
        adata_sub.uns['neighbors']['params']['metric'] = 'euclidean'
        adata_sub.uns['neighbors']['distances'] = distances
        adata_sub.uns['neighbors']['connectivities'] = connectivities

        print('run umap D_avg')
        sc.tl.umap(adata_sub)
        X_umap_avg_sub = adata_sub.obsm['X_umap']
        print('run louvain D_avg')
        sc.tl.louvain(adata_sub, key_added='louvain_avg', resolution=1.0)
        clusters_avg_sub = adata_sub.obs['louvain_avg'].astype(int)

        # now map back to the full array
        clusters_avg = [0]*adata.shape[0]
        X_umap_avg = [[0.,0.]]*adata.shape[0]
        for ii, jj in sub2full.items():
            clusters_avg[jj] = clusters_avg_sub[ii] + 1 ###### NOTE -- 0 is used for all the 'bad' clones
            X_umap_avg[jj] = list(X_umap_avg_sub[ii])
        X_umap_avg = np.array(X_umap_avg)

        adata.obsm['X_umap_avg_'+subset_tag] = X_umap_avg
        adata.obs['clusters_avg_'+subset_tag] = clusters_avg

        cluster_counts = Counter( x for x,m in zip(clusters_avg, good_mask) if m )
        print('cluster_counts:', cluster_counts.most_common())
        num_good_clusters = sum( 1 for x,y in cluster_counts.items() if y>=min_cluster_size )


        if num_good_clusters:
            # make logos using these clusters
            rank_genes_uns_tag = 'rank_genes_good_avg{}_clusters'.format(subset_tag)
            rank_genes_filt_uns_tag = rank_genes_uns_tag+'_filtered'

            clusters_tcr_fake = [1]*adata.shape[0]
            cc.run_rank_genes_on_good_cluster_pairs( adata, good_mask, clusters_avg, clusters_tcr_fake,
                                                     min_count= min_cluster_size,
                                                     key_added=rank_genes_uns_tag )

            sc.tools.filter_rank_genes_groups(adata, key=rank_genes_uns_tag, groupby='test',
                                              key_added=rank_genes_filt_uns_tag,
                                              min_in_group_fraction=0.25, min_fold_change=2,
                                              max_out_group_fraction=1.1) # 1.1 means don't filter for this

            pl.make_logo_plots( adata, nbrs_gex, nbrs_tcr, min_cluster_size,
                                '{}_good_avg{}_logos_rgfilt.png'.format(args.outfile_prefix, subset_tag),
                                clusters_gex = clusters_avg,
                                clusters_tcr = clusters_tcr_fake,
                                ignore_tcr_cluster_colors = True,
                                rank_genes_uns_tag = rank_genes_filt_uns_tag )


adata.write_h5ad(args.outfile_prefix+'_final.h5ad')
outlog.close()
