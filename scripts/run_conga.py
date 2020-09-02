######################## MAX LINE LENGTH OF ABOUT 120 ##################################################################
import sys
import os
from collections import Counter
from os.path import exists
import argparse
import time
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # in order to import conga package
import matplotlib
matplotlib.use('Agg') # for remote calcs
import conga
import scanpy as sc
import scanpy.neighbors
from sklearn.metrics import pairwise_distances
import numpy as np
import pandas as pd

start_time = time.time()
parser = argparse.ArgumentParser()

#type is str by default
parser.add_argument('--gex_data')
parser.add_argument('--gex_data_type', choices=['h5ad', '10x_mtx', '10x_h5'])
parser.add_argument('--clones_file')
parser.add_argument('--organism', choices=['mouse', 'human', 'mouse_gd', 'human_gd', 'human_ig'])
parser.add_argument('--nbr_fracs', type=float, nargs='*', default=[0.01,0.1] )
parser.add_argument('--exclude_gex_clusters', type=int, nargs='*')
parser.add_argument('--min_cluster_size', type=int, default=5)
parser.add_argument('--min_cluster_size_fraction', type=float, default=0.001)
parser.add_argument('--outfile_prefix', required=True)
parser.add_argument('--clustering_method', choices=['louvain','leiden'])
parser.add_argument('--bad_barcodes_file')
parser.add_argument('--checkpoint', action='store_true')
parser.add_argument('--restart')
parser.add_argument('--make_unfiltered_logos', action='store_true')
parser.add_argument('--make_avggood_logos', action='store_true')
parser.add_argument('--make_avgfull_logos', action='store_true')
parser.add_argument('--make_clone_plots', action='store_true')
parser.add_argument('--write_proj_info', action='store_true')
parser.add_argument('--filter_ribo_norm_low_cells', action='store_true')
# the main modes of operation
parser.add_argument('--graph_vs_graph', action='store_true')
parser.add_argument('--graph_vs_tcr_features', action='store_true')
parser.add_argument('--graph_vs_gex_features', action='store_true')
# some extra analyses
parser.add_argument('--cluster_vs_cluster', action='store_true')
parser.add_argument('--calc_clone_pmhc_pvals', action='store_true')
parser.add_argument('--find_pmhc_nbrhood_overlaps', action='store_true') # only if pmhc info is present
parser.add_argument('--find_distance_correlations', action='store_true')
parser.add_argument('--find_gex_cluster_degs', action='store_true')
# configure things
parser.add_argument('--skip_gex_header', action='store_true')
parser.add_argument('--skip_gex_header_raw', action='store_true')
parser.add_argument('--skip_gex_header_nbrZ', action='store_true')
parser.add_argument('--skip_tcr_scores_in_gex_header', action='store_true')
parser.add_argument('--tenx_agbt', action='store_true')
parser.add_argument('--include_alphadist_in_tcr_feature_logos', action='store_true')
parser.add_argument('--show_pmhc_info_in_logos', action='store_true')
parser.add_argument('--gex_header_tcr_score_names', type=str, nargs='*',
                    default= ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'])
parser.add_argument('--gex_nbrhood_tcr_score_names', type=str, nargs='*',
                    default=conga.tcr_scoring.all_tcr_scorenames )
parser.add_argument('--shuffle_tcr_kpcs', action='store_true') # shuffle the TCR kpcs to test for FDR

args = parser.parse_args()

## check consistency of args
if args.find_pmhc_nbrhood_overlaps or args.calc_clone_pmhc_pvals:
    # we need pmhc info for these analyses; right now that's restricted to the 10x AGBT dataset format
    assert args.tenx_agbt

if args.calc_clone_pmhc_pvals or args.bad_barcodes_file or args.filter_ribo_norm_low_cells:
    assert not args.restart

logfile = args.outfile_prefix+'_log.txt'
outlog = open(logfile, 'w')
outlog.write('sys.argv: {}\n'.format(' '.join(sys.argv)))
sc.logging.print_versions() # goes to stderr

if args.restart is None:

    assert exists(args.gex_data)
    assert exists(args.clones_file)

    ## load the dataset
    adata = conga.preprocess.read_dataset(args.gex_data, args.gex_data_type, args.clones_file )
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

    adata = conga.preprocess.filter_and_scale( adata )

    if args.filter_ribo_norm_low_cells:
        adata = conga.preprocess.filter_cells_by_ribo_norm( adata )

    if args.calc_clone_pmhc_pvals: # do this before condensing to a single clone per cell
        # note that we are doing this after filtering out the ribo-low cells
        results_df = conga.pmhc_scoring.calc_clone_pmhc_pvals(adata)
        tsvfile = args.outfile_prefix+'_clone_pvals.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, sep='\t', index=False)

    if args.make_clone_plots:
        pngfile = args.outfile_prefix+'_clone_plots.png'
        conga.plotting.make_clone_plots(adata, 16, pngfile)

    adata = conga.preprocess.reduce_to_single_cell_per_clone( adata )
    assert 'X_igex' in adata.obsm_keys()

    if args.shuffle_tcr_kpcs:
        X_pca_tcr = adata.obsm['X_pca_tcr']
        assert X_pca_tcr.shape[0] == adata.shape[0]
        reorder = np.random.permutation(X_pca_tcr.shape[0])
        adata.obsm['X_pca_tcr'] = X_pca_tcr[reorder,:]
        outlog.write('randomly permuting X_pca_tcr {}\n'.format(X_pca_tcr.shape))

    adata = conga.preprocess.cluster_and_tsne_and_umap( adata, clustering_method=args.clustering_method )

    if args.checkpoint:
        adata.write_h5ad(args.outfile_prefix+'_checkpoint.h5ad')
else:
    assert exists(args.restart)
    adata = sc.read_h5ad(args.restart)
    print('recover from h5ad file:', args.restart, adata )

    if 'organism' not in adata.uns_keys():
        assert args.organism
        adata.uns['organism'] = args.organism


if args.exclude_gex_clusters:
    xl = args.exclude_gex_clusters
    clusters_gex = np.array(adata.obs['clusters_gex'])
    mask = (clusters_gex==xl[0])
    for c in xl[1:]:
        mask |= (clusters_gex==c)
    print('exclude_gex_clusters: exclude {} cells in {} clusters: {}'.format(np.sum(mask), len(xl), xl))
    sys.stdout.flush()
    adata = adata[~mask,:].copy()

    adata = conga.preprocess.cluster_and_tsne_and_umap( adata )

    if args.checkpoint:
        adata.write_h5ad(args.outfile_prefix+'_checkpoint.h5ad')

################################################ DONE WITH INITIAL SETUP #########################################



# all_nbrs is dict from nbr_frac to [nbrs_gex, nbrs_tcr]
# for nndist calculations, use a smallish nbr_frac, but not too small:
num_clones = adata.shape[0]
nbr_frac_for_nndists = min( x for x in args.nbr_fracs if x*num_clones>=10 or x==max(args.nbr_fracs) )
outlog.write(f'nbr_frac_for_nndists: {nbr_frac_for_nndists}\n')
all_nbrs, nndists_gex, nndists_tcr = conga.preprocess.calc_nbrs(
    adata, args.nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)

# stash these in obs array, they are used in a few places...
adata.obs['nndists_gex'] = nndists_gex
adata.obs['nndists_tcr'] = nndists_tcr
conga.preprocess.setup_tcr_cluster_names(adata) #stores in adata.uns



if args.graph_vs_graph: ############################################################################################
    # make these numpy arrays because there seems to be a problem with np.nonzero on pandas series...
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    # run the graph vs graph analysis
    results_df = conga.correlations.run_graph_vs_graph(adata, all_nbrs)

    if results_df.shape[0]:
        # add in some extra info that may be useful before writing to tsv file
        indices = results_df['clone_index']
        results_df['gex_cluster'] = list(clusters_gex[indices])
        results_df['tcr_cluster'] = list(clusters_tcr[indices])
        for tag in 'va ja cdr3a vb jb cdr3b'.split():
            results_df[tag] = list(adata.obs[tag][indices])
        tsvfile = args.outfile_prefix+'_graph_vs_graph_hits.tsv'
        results_df.to_csv(tsvfile, sep='\t', index=False)



    # the conga scores
    conga_scores = np.array(adata.obs['conga_scores'])
    good_mask = (conga_scores <= 1.0)


    adata.obs['good_score_mask'] = good_mask

    bic_counts = Counter( (x,y) for x,y,m in zip(clusters_gex, clusters_tcr, good_mask) if m )

    # take the LARGER of the two min_cluster_size thresholds
    min_cluster_size = max( args.min_cluster_size, int( 0.5 + args.min_cluster_size_fraction * num_clones) )

    num_good_biclusters = sum( 1 for x,y in bic_counts.items() if y>=min_cluster_size )

    outlog.write(f'num_gvg_hit_clonotypes: {np.sum(good_mask)} num_gvg_hit_biclusters: {num_good_biclusters}\n')
    print('num_good_biclusters:', num_good_biclusters)

    # for the logo plots, use the largest nbr_frac
    nbrs_gex, nbrs_tcr = all_nbrs[ max(args.nbr_fracs) ]


    if num_good_biclusters:
        # calc tcr sequence features of good cluster pairs
        good_bicluster_tcr_scores = conga.correlations.calc_good_cluster_tcr_features(
            adata, good_mask, clusters_gex, clusters_tcr, args.gex_nbrhood_tcr_score_names, min_count=min_cluster_size)

        # run rank_genes on most common bics
        rank_genes_uns_tag = 'rank_genes_good_biclusters'
        conga.correlations.run_rank_genes_on_good_biclusters(
            adata, good_mask, clusters_gex, clusters_tcr, min_count=min_cluster_size, key_added= rank_genes_uns_tag)

        gex_header_tcr_score_names = [] if args.skip_tcr_scores_in_gex_header else args.gex_header_tcr_score_names

        conga.plotting.make_logo_plots(
            adata, nbrs_gex, nbrs_tcr, min_cluster_size, args.outfile_prefix+'_bicluster_logos.png',
            good_bicluster_tcr_scores=good_bicluster_tcr_scores,
            make_gex_header = not args.skip_gex_header,
            make_gex_header_raw = not args.skip_gex_header_raw,
            make_gex_header_nbrZ = not args.skip_gex_header_nbrZ,
            include_alphadist_in_tcr_feature_logos=args.include_alphadist_in_tcr_feature_logos,
            rank_genes_uns_tag = rank_genes_uns_tag,
            show_pmhc_info_in_logos = args.show_pmhc_info_in_logos,
            gex_header_tcr_score_names = gex_header_tcr_score_names )


if args.graph_vs_gex_features: #######################################################################################

    ## first use the TCRdist kPCA nbr graph:
    pval_threshold = 1.
    results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results.append( conga.correlations.tcr_nbrhood_rank_genes_fast( adata, nbrs_tcr, pval_threshold))
        results[-1]['nbr_frac'] = nbr_frac

    tsvfile = args.outfile_prefix+'_tcr_nbr_graph_vs_gex_features.tsv'
    print('making:', tsvfile)
    results_df = pd.concat(results, ignore_index=True)
    results_df.to_csv(tsvfile, index=False, sep='\t')
    tcr_nbrhood_genes_results = results_df
    combo_results = []
    if results_df.shape[0]:
        combo_results.append( results_df)


    # now make a TCR cluster graph and use the nbrhoods in there
    # make some fake nbrs-- note that only one clone per cluster has a nonempty nbrhood
    fake_nbrs_tcr = conga.correlations.setup_fake_nbrs_from_clusters_for_graph_vs_features_analysis(clusters_tcr)
    pval_threshold = 1.
    results_df = conga.correlations.tcr_nbrhood_rank_genes_fast(
        adata, fake_nbrs_tcr, pval_threshold, prefix_tag='clust')
    if results_df.shape[0]:
        results_df['clone_index'] = -1
        tsvfile = args.outfile_prefix+'_tcr_cluster_graph_vs_gex_features.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, index=False, sep='\t')
        results_df['nbr_frac'] = 0.0
        tcr_cluster_genes_results = results_df
        combo_results.append(results_df)
    else:
        tcr_cluster_genes_results = None

    if combo_results:
        results_df = pd.concat(combo_results, ignore_index=True)
        pngfile = args.outfile_prefix+'_tcr_nbr_graph_vs_gex_features.png'
        print('making:', pngfile)
        conga.plotting.plot_ranked_strings_on_cells(
            adata, results_df, 'X_tcr_2d', 'clone_index', 'mwu_pvalue_adj', 1.0, 'feature', pngfile)

        pngfile = args.outfile_prefix+'_tcr_nbr_graph_vs_gex_features_panels.png'
        print('making:', pngfile)
        conga.plotting.make_feature_panel_plots(adata, 'tcr', all_nbrs, results_df, pngfile)



    ## now make another fake nbr graph defined by TCR gene segment usage
    tcrs = conga.preprocess.retrieve_tcrs_from_adata(adata)

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
                    # this will include self but dont think thats a problem
                    fake_nbrs_tcr.append(np.nonzero( genes==g )[0] )
                    clone_display_names.append(g)

            pval_threshold = 1.

            conga.correlations.tcr_nbrhood_rank_genes_fast( adata, fake_nbrs_tcr, pval_threshold, prefix_tag=seg+ab,
                                                            clone_display_names=clone_display_names )


if args.graph_vs_tcr_features: #######################################################################################
    pval_threshold = 1.
    results = []
    for nbr_frac in args.nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results.append( conga.correlations.gex_nbrhood_rank_tcr_scores(
            adata, nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold ))
        results[-1]['nbr_frac'] = nbr_frac
    results_df = pd.concat(results, ignore_index=True)

    tsvfile = args.outfile_prefix+'_gex_nbr_graph_vs_tcr_features.tsv'
    print('making:', tsvfile)
    results_df.to_csv(tsvfile, index=False, sep='\t')
    gex_nbrhood_scores_results = results_df

    combo_results = []
    if results_df.shape[0]:
        combo_results.append(results_df)

    # make some fake nbrs
    fake_nbrs_gex = conga.correlations.setup_fake_nbrs_from_clusters_for_graph_vs_features_analysis(clusters_gex)
    pval_threshold = 1.
    results_df = conga.correlations.gex_nbrhood_rank_tcr_scores(
        adata, fake_nbrs_gex, args.gex_nbrhood_tcr_score_names, pval_threshold, prefix_tag = 'clust' )
    if results_df.shape[0]:
        results_df['clone_index'] = -1 # the clone_index values are not meaningful
        tsvfile = args.outfile_prefix+'_gex_cluster_graph_vs_tcr_features.tsv'
        print('making:', tsvfile)
        results_df.to_csv(tsvfile, index=False, sep='\t')
        results_df['nbr_frac'] = 0.0

        gex_cluster_scores_results = results_df
        combo_results.append(results_df)
    else:
        gex_cluster_scores_results = None

    if combo_results:
        pngfile = args.outfile_prefix+'_gex_nbr_graph_vs_tcr_features.png'
        print('making:', pngfile)

        results_df = pd.concat(combo_results, ignore_index=True)

        conga.plotting.plot_ranked_strings_on_cells(
            adata, results_df, 'X_gex_2d', 'clone_index', 'mwu_pvalue_adj', 1.0, 'feature', pngfile,
            direction_column='ttest_stat')

        pngfile = args.outfile_prefix+'_gex_nbr_graph_vs_tcr_features_panels.png'
        print('making:', pngfile)
        conga.plotting.make_feature_panel_plots(adata, 'gex', all_nbrs, results_df, pngfile)

if args.graph_vs_graph and args.graph_vs_tcr_features and args.graph_vs_gex_features: ################################
    pngfile = args.outfile_prefix+'_summary.png'
    print('making:', pngfile)

    if tcr_cluster_genes_results is not None:
        tcr_genes_results = pd.concat( [tcr_nbrhood_genes_results, tcr_cluster_genes_results ], ignore_index=True )
    else:
        tcr_genes_results = tcr_nbrhood_genes_results

    if gex_cluster_scores_results is not None:
        gex_scores_results = pd.concat( [gex_nbrhood_scores_results, gex_cluster_scores_results], ignore_index=True )
    else:
        gex_scores_results = gex_nbrhood_scores_results

    # default pval thresholds are .05
    conga.plotting.make_summary_figure(adata, tcr_genes_results, gex_scores_results, pngfile )


## some extra analyses
if args.cluster_vs_cluster:
    tcrs = conga.preprocess.retrieve_tcrs_from_adata(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    barcodes = list(adata.obs_names)
    barcode2tcr = dict(zip(barcodes,tcrs))
    conga.correlations.compute_cluster_interactions( clusters_gex, clusters_tcr, barcodes, barcode2tcr, outlog )

if args.find_gex_cluster_degs: # look at differentially expressed genes in gex clusters
    import matplotlib.pyplot as plt
    obs_tag = 'genex_clusters'
    adata.obs[obs_tag] = [ str(x) for x in adata.obs['clusters_gex']]#.astype('category')
    key_added = 'degs_for_gex_clusters'
    rank_method = 'wilcoxon'
    all_clusters = sorted(set(adata.obs[obs_tag]))
    sc.tl.rank_genes_groups(adata, groupby=obs_tag, method=rank_method, groups=all_clusters, reference='rest',
                            key_added=key_added)
    n_genes = 25
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, show=False, key=key_added)
    pngfile = args.outfile_prefix+'_gex_cluster_degs.png'
    plt.savefig(pngfile, bbox_inches="tight")
    print('made:', pngfile)


    new_rank_genes_genes, var_group_positions, var_group_labels = [],[],[]
    allow_gene_repeats = False
    min_rank_genes_log2fold_change = 1.0
    max_rank_genes_pval_adj=0.05
    n_genes_for_plotting = 5

    for group in all_clusters:
        my_genes = []
        for igene,gene in enumerate( adata.uns[key_added]['names'][group] ):
            log2fold = adata.uns[key_added]['logfoldchanges'][group][igene]
            pval_adj = adata.uns[key_added]['pvals_adj'][group][igene]
            #print('rank_gene:',group, igene, gene, log2fold, pval_adj)
            if len(my_genes) >= n_genes_for_plotting:
                continue
            if gene in new_rank_genes_genes and not allow_gene_repeats:
                continue # no repeats
            elif gene.startswith('MT-'):
                continue
            elif gene[:3] in ['RPL','RPS'] and gene[3].isdigit():
                continue
            elif abs(log2fold) < min_rank_genes_log2fold_change:
                continue
            elif pval_adj > max_rank_genes_pval_adj:
                continue
            print('log2fold: {:.2f} pval_adj: {:9.1e} score: {:.1f} {} {}'\
                  .format( log2fold, pval_adj, adata.uns[key_added]['scores'][group][igene],
                           gene, group ) )
            my_genes.append( gene )
        if my_genes:
            var_group_positions.append( ( len(new_rank_genes_genes),
                                          len(new_rank_genes_genes)+len(my_genes)-1 ) )
            var_group_labels.append( group )
            new_rank_genes_genes.extend( my_genes )

    if new_rank_genes_genes:
        sc.pl.stacked_violin( adata, var_names = new_rank_genes_genes, groupby=obs_tag,
                              figsize=(10,n_genes_for_plotting*10),
                              use_raw = True,
                              stripplot=True, show=False, swap_axes=True,
                              var_group_positions = var_group_positions,
                              var_group_labels = var_group_labels,
                              var_group_rotation = 1.0 )
        pngfile = args.outfile_prefix+'_gex_cluster_degs_violin.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:',pngfile)

        sc.pl.dotplot(adata, var_names=new_rank_genes_genes, groupby=obs_tag, show=False,
                      var_group_labels=var_group_labels,
                      var_group_positions=var_group_positions)
        pngfile = args.outfile_prefix+'_gex_cluster_degs_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

        sc.pl._tools.plot_scatter( adata, 'gex_2d', ncols = 6, color = new_rank_genes_genes, show=False,
                                   use_raw = True, s=40)
        pngfile = args.outfile_prefix+'_gex_cluster_degs_tsne.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)


    if adata.uns['organism'] == 'human_ig':
        genes_lines = """GC-Bs BCL6, RGS13, MEF2B, STMN1, ELL3, SERPINA9
        PBs XBP1, IRF4, SEC11C, FKBP11, JCHAIN, PRDM1
        naive TCL1A, IL4R, CCR7, IGHM, IGHD
        act-Bs TBX21, FCRL5, ITGAX, NKG7, ZEB2, CR2
        rest TNFRSF13B, CD27, CD24
        misc IGHA1 IGHA2 IGHG1 IGHG2 IGHG3 IGHG4 IGHE""".replace(',',' ').split('\n')
        genes, var_group_positions, var_group_labels = [], [], []
        for line in genes_lines:
            my_genes = [ x for x in line.split()[1:] if x in adata.raw.var_names]
            print(len(my_genes), line.split())
            if my_genes:
                var_group_positions.append( (len(genes), len(genes)+len(my_genes)-1) )
                var_group_labels.append( line.split()[0])
                genes.extend(my_genes)
        sc.pl.dotplot(adata, var_names=genes, groupby=obs_tag, show=False, var_group_labels=var_group_labels,
                      var_group_positions=var_group_positions)
        pngfile = args.outfile_prefix+'_gex_cluster_bcell_genes_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

# just out of curiosity:
conga.correlations.check_nbr_graphs_indegree_bias(all_nbrs)

if args.find_distance_correlations:
    pvalues, rvalues = conga.correlations.compute_distance_correlations(adata)
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
    agroups, bgroups = conga.preprocess.setup_tcr_groups(adata)

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


if args.write_proj_info:
    outfile = args.outfile_prefix+'_2d_proj_info.txt'
    conga.preprocess.write_proj_info( adata, outfile )

adata.write_h5ad(args.outfile_prefix+'_final.h5ad')
adata.obs.to_csv(args.outfile_prefix+'_final_obs.tsv', sep='\t')

outlog.write('run_conga took {:.3f} minutes\n'.format((time.time()- start_time)/60))

outlog.close()
print('DONE')
