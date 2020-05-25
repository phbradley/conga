import sys
import os
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )
from collections import Counter
from os.path import exists
import argparse
import conga
import time
import conga.preprocess as pp
import conga.correlations as cc
import conga.plotting as pl
import scanpy as sc
import scanpy.neighbors
from sklearn.metrics import pairwise_distances
import numpy as np
import pandas as pd
import scipy.stats

start_time = time.time()
parser = argparse.ArgumentParser()

#type is str by default
parser.add_argument('--gex_data')
parser.add_argument('--gex_data_type', choices=['h5ad', '10x_mtx', '10x_h5'])
parser.add_argument('--clones_file')
parser.add_argument('--organism', choices=['mouse', 'human'])
parser.add_argument('--nbr_fracs', type=float, nargs='*', default=[0.01,0.1] )
parser.add_argument('--exclude_gex_clusters', type=int, nargs='*')
parser.add_argument('--min_cluster_size', type=int, default=5)
parser.add_argument('--min_cluster_size_fraction', type=float, default=0.001)
parser.add_argument('--outfile_prefix', required=True)
parser.add_argument('--bad_barcodes_file')
parser.add_argument('--checkpoint', action='store_true')
parser.add_argument('--from_checkpoint1')
parser.add_argument('--make_unfiltered_logos', action='store_true')
parser.add_argument('--make_avggood_logos', action='store_true')
parser.add_argument('--make_avgfull_logos', action='store_true')
parser.add_argument('--make_clone_plots', action='store_true')
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
parser.add_argument('--find_distance_correlations', action='store_true')

parser.add_argument('--skip_gex_header', action='store_true')
parser.add_argument('--skip_gex_header_raw', action='store_true')
parser.add_argument('--skip_gex_header_nbrZ', action='store_true')
parser.add_argument('--skip_tcr_scores_in_gex_header', action='store_true')
parser.add_argument('--tenx_agbt', action='store_true')
parser.add_argument('--include_alphadist_in_tcr_feature_logos', action='store_true')
parser.add_argument('--show_pmhc_info_in_logos', action='store_true')
parser.add_argument('--gex_header_tcr_score_names', type=str, nargs='*',
                    default= ['mhci2', 'cdr3len', 'cd8', 'nndists_tcr'])
parser.add_argument('--gex_nbrhood_tcr_score_names', type=str, nargs='*',
                    default=conga.tcr_scoring.all_tcr_scorenames )
parser.add_argument('--shuffle_tcr_kpcs', action='store_true') # shuffle the TCR kpcs to test for FDR

args = parser.parse_args()

## check consistency of args
if args.find_pmhc_nbrhood_overlaps or args.calc_clone_pmhc_pvals:
    # we need pmhc info for these analyses; right now that's restricted to the 10x AGBT dataset format
    assert args.tenx_agbt

if args.calc_clone_pmhc_pvals or args.bad_barcodes_file or args.filter_ribo_norm_low_cells:
    assert not args.from_checkpoint1


logfile = args.outfile_prefix+'_log.txt'
outlog = open(logfile, 'w')
outlog.write('sys.argv: {}\n'.format(' '.join(sys.argv)))

if args.from_checkpoint1 is None:

    assert exists(args.gex_data)
    assert exists(args.clones_file)

    ## load the dataset
    adata = pp.read_dataset(args.gex_data, args.gex_data_type, args.clones_file )
    assert args.organism
    adata.uns['organism'] = args.organism
    assert 'organism' in adata.uns_keys()

    if args.tenx_agbt:
        conga.pmhc_scoring.shorten_pmhc_var_names(adata)

        adata.uns['pmhc_var_names'] = conga.pmhc_scoring.get_tenx_agbt_pmhc_var_names(adata)
        print('pmhc_var_names:', adata.uns['pmhc_var_names'])

    if args.bad_barcodes_file:
        bad_barcodes = frozenset([x[:-1] for x in open(args.bad_barcodes_file,'rU')])
        bad_bc_mask = np.array( [x in bad_barcodes for x in adata.obs_names ] )
        num_bad = np.sum(bad_bc_mask)
        if num_bad:
            print('excluding {} bad barcodes found in {}'\
                  .format(num_bad, args.bad_barcodes_file))
            adata = adata[~bad_bc_mask,:].copy()
        else:
            print('WARNING:: no matched barcodes in bad_barcodes_file: {}'.format(args.bad_barcodes_file))


    assert not adata.isview
    assert 'X_pca_tcr' in adata.obsm_keys() # tcr-dist kPCA info
    assert 'cdr3a' in adata.obs # tcr sequence (VDJ) info (plus other obs keys)

    print(adata)

    adata = pp.filter_and_scale( adata )

    if args.filter_ribo_norm_low_cells:
        adata = pp.filter_cells_by_ribo_norm( adata )

    if args.calc_clone_pmhc_pvals: # do this before condensing to a single clone per cell
        # note that we are doing this after filtering out the ribo-low cells
        results_df = conga.pmhc_scoring.calc_clone_pmhc_pvals(adata)
        tsvfile = args.outfile_prefix+'_clone_pvals.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, sep='\t', index=False)

    if args.make_clone_plots:
        pngfile = args.outfile_prefix+'_clone_plots.png'
        pl.make_clone_plots(adata, 16, pngfile)

    adata = pp.reduce_to_single_cell_per_clone( adata )
    assert 'X_igex' in adata.obsm_keys()

    if args.shuffle_tcr_kpcs:
        X_pca_tcr = adata.obsm['X_pca_tcr']
        assert X_pca_tcr.shape[0] == adata.shape[0]
        reorder = np.random.permutation(X_pca_tcr.shape[0])
        adata.obsm['X_pca_tcr'] = X_pca_tcr[reorder,:]
        outlog.write('randomly permuting X_pca_tcr {}\n'.format(X_pca_tcr.shape))

    adata = pp.cluster_and_tsne_and_umap( adata )

    if args.checkpoint:
        adata.write_h5ad(args.outfile_prefix+'_checkpoint1.h5ad')
else:
    adata = sc.read_h5ad(args.from_checkpoint1)
    print('recover from checkpoint:', adata )

    if 'organism' not in adata.uns_keys():
        assert args.organism
        adata.uns['organism'] = args.organism


if args.exclude_gex_clusters:
    xl = args.exclude_gex_clusters
    clusters_gex = adata.obs['clusters_gex']
    mask = (clusters_gex==xl[0])
    for c in xl[1:]:
        mask |= (clusters_gex==c)
    print('exclude_gex_clusters: exclude {} cells in {} clusters: {}'.format(np.sum(mask), len(xl), xl))
    sys.stdout.flush()
    adata = adata[~mask,:].copy()

    adata = pp.cluster_and_tsne_and_umap( adata )

    if args.checkpoint:
        adata.write_h5ad(args.outfile_prefix+'_checkpoint1.h5ad')

################################################ DONE WITH INITIAL SETUP #########################################


if args.write_proj_info:
    outfile = args.outfile_prefix+'_2d_proj_info.txt'
    pp.write_proj_info( adata, outfile )

pp.setup_tcr_cluster_names(adata) #stores in adata.uns

clusters_gex = adata.obs['clusters_gex']
clusters_tcr = adata.obs['clusters_tcr']

num_clones = adata.shape[0]

#barcode2tcr = { x:y for x,y in zip( adata.obs_names, tcrs )}
#cc.compute_cluster_interactions( clusters_gex, clusters_tcr, adata.obs_names, barcode2tcr, outlog )

agroups, bgroups = pp.setup_tcr_groups(adata)

# all_nbrs is dict from nbr_frac to [nbrs_gex, nbrs_tcr]
# for nndist calculations, use a smallish nbr_frac, but not too small:
nbr_frac_for_nndists = min( x for x in args.nbr_fracs if x*num_clones>=10 or x==max(args.nbr_fracs) )
outlog.write(f'nbr_frac_for_nndists: {nbr_frac_for_nndists}\n')
all_nbrs, nndists_gex, nndists_tcr = pp.calc_nbrs( adata, args.nbr_fracs, also_calc_nndists=True,
                                                   nbr_frac_for_nndists=nbr_frac_for_nndists )
adata.obs['nndists_gex'] = nndists_gex
adata.obs['nndists_tcr'] = nndists_tcr

bad_conga_score = -1*np.log10(num_clones)
conga_scores = np.full( (num_clones,3), bad_conga_score)



if args.find_nbrhood_overlaps:
    all_results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]

        # look at the in-degree distribution
        expected_indegree = nbrs_gex.shape[1]
        gex_counts = Counter(nbrs_gex.flatten())
        gex_indegree_bias = np.array( [ gex_counts[x]/expected_indegree for x in range(num_clones) ] )
        print(f'gex_indegree_bias: nbr_frac= {nbr_frac:.4f}', scipy.stats.describe(gex_indegree_bias))
        tcr_counts = Counter(nbrs_tcr.flatten())
        tcr_indegree_bias = np.array( [ tcr_counts[x]/expected_indegree for x in range(num_clones) ] )
        print('tcr_indegree_bias: nbr_frac= {nbr_frac:.4f}', scipy.stats.describe(tcr_indegree_bias))
        # any correlation?
        print('indegree_bias_correlation:', scipy.stats.linregress(gex_indegree_bias, tcr_indegree_bias))


        pval_threshold = 1. # since they are being scaled by num_clones
        print('find_neighbor_neighbor_interactions:')
        results_df = cc.find_neighbor_neighbor_interactions(
            adata, nbrs_gex, nbrs_tcr, agroups, bgroups, pval_threshold)
        if not results_df.empty:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'nbr_nbr'
            all_results.append(results_df)
            for r in results_df.itertuples():
                conga_scores[r.clone_index,0] = max(conga_scores[r.clone_index,0], -1*np.log10(r.conga_score) )

        print('find_neighbor_cluster_interactions:')
        results_df = cc.find_neighbor_cluster_interactions(
            adata, nbrs_tcr, clusters_gex, agroups, bgroups, pval_threshold)
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'cluster_nbr'
            all_results.append(results_df)
            for r in results_df.itertuples():
                conga_scores[r.clone_index,1] = max(conga_scores[r.clone_index,1], -1*np.log10(r.conga_score) )

        print('find_neighbor_cluster_interactions:')
        results_df = cc.find_neighbor_cluster_interactions(
            adata, nbrs_gex, clusters_tcr, agroups, bgroups, pval_threshold)
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'nbr_cluster'
            all_results.append(results_df)
            for r in results_df.itertuples():
                conga_scores[r.clone_index,2] = max(conga_scores[r.clone_index,2], -1*np.log10(r.conga_score) )

    results_df = pd.concat(all_results, ignore_index=True)
    if results_df.shape[0]:
        indices = results_df['clone_index']
        results_df['gex_cluster'] = list(clusters_gex[indices])
        results_df['tcr_cluster'] = list(clusters_tcr[indices])
        for tag in 'va ja cdr3a vb jb cdr3b'.split():
            results_df[tag] = list(adata.obs[tag][indices])
        tsvfile = args.outfile_prefix+'_graph_graph_overlaps.tsv'
        results_df.to_csv(tsvfile, sep='\t', index=False)

if args.find_distance_correlations:
    pvalues, rvalues = cc.compute_distance_correlations(adata)
    results = []
    for ii, (pval, rval) in enumerate(zip(rvalues, pvalues)):
        if pval<1:
            results.append( dict( clone_index=ii, pvalue_adj=pval, rvalue=rval, gex_cluster=clusters_gex[ii],
                                  tcr_cluster=clusters_tcr[ii]))
    if results:
        results_df = pd.DataFrame(results)
        outfile = args.outfile_prefix+'_distance_correlations.tsv'
        results_df.to_csv(outfile, sep='\t', index=False)

if args.find_pmhc_nbrhood_overlaps:
    pmhc_nbrhood_overlap_results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        for tag, nbrs in [['gex', nbrs_gex], ['tcr', nbrs_tcr]]:
            results_df = conga.pmhc_scoring.compute_pmhc_versus_nbrs(adata, nbrs, agroups, bgroups )
            results_df['nbr_tag'] = tag
            results_df['nbr_frac'] = nbr_frac
            pmhc_nbrhood_overlap_results.append( results_df )

    tsvfile = args.outfile_prefix+'_pmhc_versus_nbrs.tsv'
    print('making:', tsvfile)
    pd.concat(pmhc_nbrhood_overlap_results).to_csv(tsvfile, index=False, sep='\t')

adata.obsm['conga_scores'] = conga_scores

max_conga_score = np.max( conga_scores, axis=1 )
good_mask = (max_conga_score >= 0.)

outlog.write('num_good: {}\n'.format(np.sum(good_mask)))

adata.obs['good_score_mask'] = good_mask

# take the LARGER of the two min_cluster_size thresholds
min_cluster_size = max( args.min_cluster_size, int( 0.5 + args.min_cluster_size_fraction * num_clones) )


# if args.checkpoint:
#     adata.write_h5ad(args.outfile_prefix+'_checkpoint2.h5ad')
if args.find_tcr_nbrhood_genes:
    pval_threshold = 1.
    results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results.append( cc.tcr_nbrhood_rank_genes_fast( adata, nbrs_tcr, pval_threshold))
        results[-1]['nbr_frac'] = nbr_frac

    tsvfile = args.outfile_prefix+'_tcr_nbrhood_genes.tsv'
    print('making:', tsvfile)
    results_df = pd.concat(results, ignore_index=True)
    results_df.to_csv(tsvfile, index=False, sep='\t')
    tcr_nbrhood_genes_results = results_df
    if results_df.shape[0]:
        pngfile = args.outfile_prefix+'_tcr_nbrhood_genes.png'
        print('making:', pngfile)
        exclude_strings = [pp.FUNNY_MOUSE_V_GENE] # bad mouse gene, actually a tcr v gene
        pl.plot_ranked_strings_on_cells(adata, results_df, 'X_tcr_2d', 'clone_index', 'mwu_pvalue_adj', 1.0, 'feature',
                                        pngfile, exclude_strings=exclude_strings)

        pngfile = args.outfile_prefix+'_tcr_nbrhood_gene_panels.png'
        print('making:', pngfile)
        #exclude_strings = ['5830405F06Rik'] # bad mouse gene, actually a tcr v gene
        pl.make_feature_panel_plots(adata, 'tcr', all_nbrs, results_df, pngfile)


if args.find_tcr_cluster_genes:
    # make some fake nbrs
    fake_nbrs_tcr = []
    seen = set()
    for cl in clusters_tcr:
        if cl in seen:
            fake_nbrs_tcr.append([])
        else:
            seen.add(cl)
            fake_nbrs_tcr.append(np.nonzero( clusters_tcr==cl )[0])
    pval_threshold = 1.

    results_df = cc.tcr_nbrhood_rank_genes_fast( adata, fake_nbrs_tcr, pval_threshold, prefix_tag='clust')
    if results_df.shape[0]:
        results_df['clone_index'] = -1
        tsvfile = args.outfile_prefix+'_tcr_cluster_genes.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, index=False, sep='\t')
        results_df['nbr_frac'] = 0.0
        tcr_cluster_genes_results = results_df
    else:
        tcr_cluster_genes_results = None

if args.find_tcr_segment_genes:
    tcrs = pp.retrieve_tcrs_from_adata(adata)

    for iab,ab in enumerate('AB'):
        for iseg,seg in enumerate('VJ'):
            genes = [ x[iab][iseg] for x in tcrs ]
            genes = np.array([ x[:x.index('*')] for x in genes ])

            # make some fake nbrs
            fake_nbrs_tcr = []
            clone_display_names = []
            seen = set()
            for g in genes:
                if g in seen:
                    fake_nbrs_tcr.append([])
                    clone_display_names.append('')
                else:
                    seen.add(g)
                    fake_nbrs_tcr.append(np.nonzero( genes==g )[0] ) # this will include self but dont think thats a prob
                    clone_display_names.append(g)

            pval_threshold = 1.

            cc.tcr_nbrhood_rank_genes_fast( adata, fake_nbrs_tcr, pval_threshold, prefix_tag=seg+ab,
                                            clone_display_names=clone_display_names )



if args.find_gex_nbrhood_scores:
    pval_threshold = 1.
    results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results.append( cc.gex_nbrhood_rank_tcr_scores(
            adata, nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold ))
        results[-1]['nbr_frac'] = nbr_frac
    results_df = pd.concat(results, ignore_index=True)

    tsvfile = args.outfile_prefix+'_gex_nbrhood_scores.tsv'
    print('making:', tsvfile)
    results_df.to_csv(tsvfile, index=False, sep='\t')
    gex_nbrhood_scores_results = results_df

    if results_df.shape[0]:
        pngfile = args.outfile_prefix+'_gex_nbrhood_scores.png'
        print('making:', pngfile)

        pl.plot_ranked_strings_on_cells(adata, results_df, 'X_gex_2d', 'clone_index', 'mwu_pvalue_adj', 1.0,
                                        'feature', pngfile, direction_column='ttest_stat')

        pngfile = args.outfile_prefix+'_gex_nbrhood_score_panels.png'
        print('making:', pngfile)
        pl.make_feature_panel_plots(adata, 'gex', all_nbrs, results_df, pngfile)


if args.find_gex_cluster_scores:
    # make some fake nbrs
    fake_nbrs_gex = []
    seen = set()
    for cl in clusters_gex:
        if cl in seen:
            fake_nbrs_gex.append([])
        else:
            seen.add(cl)
            fake_nbrs_gex.append(np.nonzero( clusters_gex==cl )[0])

    pval_threshold = 1.
    results_df = cc.gex_nbrhood_rank_tcr_scores( adata, fake_nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold,
                                                 prefix_tag = 'clust' )
    if results_df.shape[0]:
        results_df['clone_index'] = -1 # the clone_index values are not meaningful
        tsvfile = args.outfile_prefix+'_gex_cluster_scores.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, index=False, sep='\t')
        results_df['nbr_frac'] = 0.0

        gex_cluster_scores_results = results_df
    else:
        gex_cluster_scores_results = None


if args.find_nbrhood_overlaps and args.find_tcr_nbrhood_genes and args.find_gex_nbrhood_scores:
    pngfile = args.outfile_prefix+'_summary.png'
    print('making:', pngfile)

    if args.find_tcr_cluster_genes and tcr_cluster_genes_results is not None:
        tcr_genes_results = pd.concat( [tcr_nbrhood_genes_results, tcr_cluster_genes_results ], ignore_index=True )
    else:
        tcr_genes_results = tcr_nbrhood_genes_results

    if args.find_gex_cluster_scores and gex_cluster_scores_results is not None:
        gex_scores_results = pd.concat( [gex_nbrhood_scores_results, gex_cluster_scores_results], ignore_index=True )
    else:
        gex_scores_results = gex_nbrhood_scores_results

    # default pval thresholds are .05
    pl.make_summary_figure(adata, tcr_genes_results, gex_scores_results, pngfile )


clp_counts = Counter( (x,y) for x,y,m in zip(clusters_gex, clusters_tcr, good_mask) if m )
num_good_cluster_pairs = sum( 1 for x,y in clp_counts.items() if y>=min_cluster_size )

print('num_good_cluster_pairs:', num_good_cluster_pairs)

# for the logo plots, use the largest nbr_frac
nbrs_gex, nbrs_tcr = all_nbrs[ max(args.nbr_fracs) ]


if num_good_cluster_pairs:
    # calc tcr sequence features of good cluster pairs
    good_cluster_pair_tcr_scores = cc.calc_good_cluster_tcr_features( adata, good_mask, clusters_gex, clusters_tcr,
                                                                      args.gex_nbrhood_tcr_score_names,
                                                                      min_count=min_cluster_size )

    # run rank_genes on most common clps
    rank_genes_uns_tag = 'rank_genes_good_cluster_pairs'
    cc.run_rank_genes_on_good_cluster_pairs( adata, good_mask, clusters_gex, clusters_tcr, min_count=min_cluster_size,
                                             key_added= rank_genes_uns_tag)

    # not currently using these filtered results:
    rank_genes_filt_uns_tag = 'rank_genes_good_cluster_pairs_filtered'

    sc.tools.filter_rank_genes_groups(adata, key=rank_genes_uns_tag, groupby='test',
                                      key_added=rank_genes_filt_uns_tag,
                                      min_in_group_fraction=0.25, min_fold_change=2,
                                      max_out_group_fraction=1.1) # 1.1 means don't filter for this (dflt was 0.5)

    gex_header_tcr_score_names = [] if args.skip_tcr_scores_in_gex_header else args.gex_header_tcr_score_names

    pl.make_logo_plots( adata, nbrs_gex, nbrs_tcr, min_cluster_size, args.outfile_prefix+'_good_logos_rgfilt.png',
                        good_cluster_pair_tcr_scores=good_cluster_pair_tcr_scores,
                        make_gex_header = not args.skip_gex_header,
                        make_gex_header_raw = not args.skip_gex_header_raw,
                        make_gex_header_nbrZ = not args.skip_gex_header_nbrZ,
                        include_alphadist_in_tcr_feature_logos=args.include_alphadist_in_tcr_feature_logos,
                        rank_genes_uns_tag = rank_genes_uns_tag,
                        show_pmhc_info_in_logos = args.show_pmhc_info_in_logos,
                        gex_header_tcr_score_names = gex_header_tcr_score_names )
    #rank_genes_uns_tag = rank_genes_filt_uns_tag )

    # if args.make_unfiltered_logos:
    #     pl.make_logo_plots( adata, nbrs_gex, nbrs_tcr, min_cluster_size, args.outfile_prefix+'_good_logos.png',
    #                         make_gex_header = not args.skip_gex_header,
    #                         rank_genes_uns_tag = rank_genes_uns_tag )

# try clustering using an averaged distance matrix
if np.sum(good_mask) >= min_cluster_size and (args.make_avgfull_logos or args.make_avggood_logos):
    for subset_tag in ['good','full']:
        if subset_tag=='good':
            if not args.make_avggood_logos:
                continue
            if np.sum(good_mask) < min_cluster_size:
                continue
            adata_sub = adata[good_mask,:].copy()
            sub2full = dict( enumerate( np.nonzero(good_mask)[0] ) )

        else:
            if not args.make_avgfull_logos:
                continue
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
                                make_gex_header = not args.skip_gex_header,
                                ignore_tcr_cluster_colors = True,
                                rank_genes_uns_tag = rank_genes_filt_uns_tag )


adata.write_h5ad(args.outfile_prefix+'_final.h5ad')
adata.obs.to_csv(args.outfile_prefix+'_final_obs.tsv', sep='\t')

outlog.write('run_conga took {:.3f} minutes\n'.format((time.time()- start_time)/60))

outlog.close()
