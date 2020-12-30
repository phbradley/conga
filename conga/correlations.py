import numpy as np
from scipy import stats
from sklearn.metrics import pairwise_distances
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm
#from scipy.sparse import issparse, csr_matrix
import scipy.sparse as sps
from collections import Counter, OrderedDict
import scanpy as sc
from . import preprocess as pp
from . import tcr_scoring
from . import util
from .tcrdist.all_genes import all_genes
import sys
import pandas as pd
from sys import exit
import time #debugging

def find_neighbor_neighbor_interactions(
        adata,
        nbrs_gex,
        nbrs_tcr,
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True,
        verbose=False
):
    ''' Returns a pandas dataframe with results
    of all clones with nbr-nbr overlaps having pvalues <= pval_threshold

    AND a numpy array of the adjusted_pvals

    '''
    pp.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_mait']

    num_clones = len(nbrs_gex)

    pval_rescale = num_clones if scale_pvals_by_num_clones else 1.0

    results = []
    adjusted_pvalues = [pval_rescale]*num_clones # initialize to big value

    for ii in range(num_clones):
        is_ii_group = (agroups==agroups[ii]) | (bgroups==bgroups[ii])
        assert is_ii_group[ii]
        actual_num_clones = num_clones - np.sum( is_ii_group ) - counts_correction

        ii_nbrs_gex = nbrs_gex[ii] # might be slice of array, or list entry if use_sym_nbrs
        ii_nbrs_tcr = nbrs_tcr[ii]
        ii_nbrs_gex_set = frozenset( ii_nbrs_gex)
        assert ii not in ii_nbrs_gex_set

        num_neighbors_gex = len(ii_nbrs_gex)
        num_neighbors_tcr = len(ii_nbrs_tcr)

        overlap = sum( 1 for x in ii_nbrs_tcr if x in ii_nbrs_gex_set )
        expected_overlap = float( num_neighbors_gex*num_neighbors_tcr )/actual_num_clones
        if overlap and overlap > expected_overlap:
            nbr_pval = pval_rescale * hypergeom.sf( overlap-1, actual_num_clones, num_neighbors_gex, num_neighbors_tcr)
        else:
            nbr_pval = pval_rescale

        if verbose:
            print('nbr_overlap:', ii, overlap, nbr_pval)

        adjusted_pvalues[ii] = nbr_pval

        if nbr_pval > pval_threshold:
            continue ## NOTE

        double_nbrs = [ x for x in ii_nbrs_tcr if x in ii_nbrs_gex_set ]
        assert overlap == len(double_nbrs)

        if verbose:
            print('nbr_overlap_nbrs:', ii, overlap, nbr_pval, double_nbrs)

        overlap_corrected = overlap
        if correct_overlaps_for_groups:
            overlap_corrected = min( len(set(agroups[double_nbrs])), len(set(bgroups[double_nbrs])) )

            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                nbr_pval = pval_rescale * hypergeom.sf( overlap_corrected-1, actual_num_clones, num_neighbors_gex-delta,
                                                        num_neighbors_tcr-delta)
                adjusted_pvalues[ii] = nbr_pval # update
                if nbr_pval > pval_threshold:
                    continue ## NOTE

        results.append( dict( conga_score=nbr_pval,
                              num_neighbors_gex=num_neighbors_gex,
                              num_neighbors_tcr=num_neighbors_tcr,
                              overlap=overlap,
                              overlap_corrected=overlap_corrected,
                              mait_fraction=np.sum(is_mait[double_nbrs])/overlap,
                              clone_index=ii ))#double_nbrs ] )

    return pd.DataFrame(results), np.array(adjusted_pvalues)


def find_neighbor_cluster_interactions(
        adata,
        nbrs,
        clusters, # for example, or could be swapped: tcr/gex
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number of them
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True
):
    ''' Returns a pandas dataframe with results
    of all clones with nbr-cluster overlaps having pvalues <= pval_threshold

    AND a numpy array of the adjusted_pvals

    '''
    pp.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_mait']

    num_clones = len(nbrs)

    pval_rescale = num_clones if scale_pvals_by_num_clones else 1.0

    results = []
    adjusted_pvalues = [pval_rescale]*num_clones # initialize to big pval

    for ii in range(num_clones):
        is_ii_group = (agroups==agroups[ii]) | (bgroups==bgroups[ii])
        assert is_ii_group[ii]
        actual_num_clones = num_clones - np.sum( is_ii_group ) - counts_correction

        ii_nbrs = nbrs[ii]
        ii_cluster = clusters[ii]
        ii_cluster_clustersize = np.sum(clusters==ii_cluster) - np.sum( (clusters==ii_cluster)&(is_ii_group) )

        num_neighbors = ii_nbrs.shape[0]
        overlap = np.sum( clusters[ii_nbrs] == ii_cluster )
        expected = float(ii_cluster_clustersize*num_neighbors)/actual_num_clones
        if overlap and overlap>expected:
            nbr_pval = pval_rescale * hypergeom.sf( overlap-1, actual_num_clones, num_neighbors, ii_cluster_clustersize )
        else:
            nbr_pval = pval_rescale

        adjusted_pvalues[ii] = nbr_pval
        if nbr_pval > pval_threshold:
            continue ## NOTE

        same_cluster_nbrs = ii_nbrs[ clusters[ii_nbrs] == ii_cluster ]
        assert len(same_cluster_nbrs)== overlap

        overlap_corrected = overlap
        if correct_overlaps_for_groups:
            overlap_corrected = min( len(set(agroups[same_cluster_nbrs])), len(set(bgroups[same_cluster_nbrs])) )
            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                nbr_pval = pval_rescale * hypergeom.sf(overlap_corrected-1, actual_num_clones, num_neighbors-delta,
                                                       ii_cluster_clustersize-delta )
                adjusted_pvalues[ii] = nbr_pval
                if nbr_pval > pval_threshold:
                    continue

        results.append( dict( conga_score=nbr_pval,
                              num_neighbors=num_neighbors,
                              cluster_size=ii_cluster_clustersize,
                              overlap=overlap,
                              overlap_corrected=overlap_corrected,
                              mait_fraction=np.sum(is_mait[same_cluster_nbrs])/overlap,
                              clone_index=ii ))#double_nbrs ] )

    return pd.DataFrame(results), np.array(adjusted_pvalues)



def check_nbr_graphs_indegree_bias(all_nbrs):
    ''' all_nbrs is a dict mapping from nbr_frac to (gex_nbrs, tcr_nbrs) setup by preprocess.calc_nbrs
    '''
    for nbr_frac in all_nbrs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        num_clones = nbrs_gex.shape[0]

        # look at the in-degree distribution
        expected_indegree = nbrs_gex.shape[1]
        gex_counts = Counter(nbrs_gex.flatten())
        gex_indegree_bias = np.array( [ gex_counts[x]/expected_indegree for x in range(num_clones) ] )
        print(f'gex_indegree_bias: nbr_frac= {nbr_frac:.4f}', stats.describe(gex_indegree_bias))

        tcr_counts = Counter(nbrs_tcr.flatten())
        tcr_indegree_bias = np.array( [ tcr_counts[x]/expected_indegree for x in range(num_clones) ] )
        print(f'tcr_indegree_bias: nbr_frac= {nbr_frac:.4f}', stats.describe(tcr_indegree_bias))

        # any correlation?
        # if there is strong correlation, this could skew the graph-vs-graph analysis
        print(f'indegree_bias_correlation: nbr_frac= {nbr_frac:.4f}',
              stats.linregress(gex_indegree_bias, tcr_indegree_bias))


def run_graph_vs_graph(
        adata,
        all_nbrs,
        pval_threshold=1.0, #pvals are crude bonferroni corrected by multiplying by num_clones, ie they can be >> 1
        verbose=False
):
    ''' Runs graph-vs-graph analysis for each nbr_frac in the all_nbrs dictionary

    Returns a (possibly empty) pandas dataframe with the results

    Also sets up the

    'conga_scores'  array in adata.obs

    '''

    num_clones = adata.shape[0]
    agroups, bgroups = pp.setup_tcr_groups(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    bad_conga_score = num_clones # since conga_score= raw_pval*num_clones, max is num_clones if raw_pval=1

    conga_scores = np.full( (num_clones,), bad_conga_score, dtype=float)

    all_results = []

    for nbr_frac in all_nbrs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]


        print('find_neighbor_neighbor_interactions:')
        results_df, adjusted_pvalues = find_neighbor_neighbor_interactions(
            adata, nbrs_gex, nbrs_tcr, agroups, bgroups, pval_threshold, verbose=verbose)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)

        if not results_df.empty:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'nbr_nbr'
            all_results.append(results_df)

        print('find_neighbor_cluster_interactions:')
        results_df, adjusted_pvalues = find_neighbor_cluster_interactions(
            adata, nbrs_tcr, clusters_gex, agroups, bgroups, pval_threshold)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'cluster_nbr'
            all_results.append(results_df)

        print('find_neighbor_cluster_interactions:')
        results_df, adjusted_pvalues = find_neighbor_cluster_interactions(
            adata, nbrs_gex, clusters_tcr, agroups, bgroups, pval_threshold)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['overlap_type'] = 'nbr_cluster'
            all_results.append(results_df)

    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)
    else:
        results_df = pd.DataFrame([])

    adata.obs['conga_scores'] = conga_scores

    return results_df


def run_rank_genes_on_good_biclusters(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        rank_method='wilcoxon',
        rg_tag = 'test',
        neg_tag='none',
        min_count=5,
        key_added = 'rank_genes_good_biclusters'
):
    num_clones = adata.shape[0]

    clp_counts = Counter( (x,y) for x,y,z in zip( clusters_gex, clusters_tcr, good_mask ) if z )
    print( clp_counts.most_common())

    vals = [ neg_tag ]*num_clones

    for clp, count in clp_counts.most_common():
        if count< min_count:
            break
        tag = 'clp_{}_{}'.format(clp[0], clp[1])
        inds = np.nonzero( [ x==clp[0] and y==clp[1] and z for x,y,z in zip(clusters_gex, clusters_tcr, good_mask)])[0]
        for ii in inds:
            vals[ii] = tag

    pos_tags = set(vals)
    if neg_tag in pos_tags:
        pos_tags.remove(neg_tag)

    if not pos_tags:
        print('run_rank_genes_on_good_biclusters: no good cluster pairs')
        return

    adata.obs[rg_tag] = vals

    print('run rank_genes_groups', Counter(vals).most_common() )

    sc.tl.rank_genes_groups( adata, groupby=rg_tag, method=rank_method, groups=pos_tags, reference='rest',
                             key_added = key_added )


def calc_good_cluster_tcr_features(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        tcr_score_names,
        min_count=5,
        verbose=True,
):
    num_clones = adata.shape[0]

    # there seems to be a problem with np.nonzero on a pandas series, which is what these might be if
    # taken from adata.obs
    clusters_gex = np.array(clusters_gex)
    clusters_tcr = np.array(clusters_tcr)
    good_mask = np.array(good_mask)

    clp_counts = Counter( (x,y) for x,y,z in zip( clusters_gex, clusters_tcr, good_mask ) if z )
    print( clp_counts.most_common())

    good_clps = [ x for x,y in clp_counts.items() if y>=min_count]

    fake_nbrs_gex = []
    seen = set()
    for cl_gex, cl_tcr, m in zip( clusters_gex, clusters_tcr, good_mask ):
        clp=(cl_gex, cl_tcr)
        if m and clp in good_clps and clp not in seen:
            seen.add(clp)
            fake_nbrs_gex.append( np.nonzero( (clusters_gex==cl_gex) & (clusters_tcr==cl_tcr) &good_mask )[0] )
        else:
            fake_nbrs_gex.append([])

    pval_threshold = 1.
    results_df = gex_nbrhood_rank_tcr_scores( adata, fake_nbrs_gex, tcr_score_names, pval_threshold,
                                              verbose=verbose, prefix_tag = 'good_clp' )


    all_tcr_features = {}
    for clp in good_clps:
        all_tcr_features[clp] = [] # in case results_df doesn't have any rows for this clp

    for row in results_df.itertuples():
        #clp = (row.gex_cluster, row.tcr_cluster)
        clp = (clusters_gex[row.clone_index], clusters_tcr[row.clone_index])
        assert clp in good_clps
        all_tcr_features[clp].append( ( row.mwu_pvalue_adj, row.ttest_stat, row.feature))


    for clp in all_tcr_features:
        # plotting code expects name then stat then pvalue
        # but we want them sorted by significance
        all_tcr_features[clp] = [ (x[2], x[1], x[0]) for x in sorted(all_tcr_features[clp]) ]

    return all_tcr_features




def get_split_mean_var( X, X_sq, mask, mean, mean_sq ):
    # use: var = (mean_sq - mean**2)
    #
    N = mask.shape[0]
    assert X.shape[0] == N
    wt_fg = np.sum(mask) / N
    wt_bg = 1. - wt_fg
    mean_fg = X[mask].mean(axis=0)
    mean_sq_fg = X_sq[mask].mean(axis=0)
    if sps.issparse(X):
        mean_fg = mean_fg.A1
        mean_sq_fg = mean_sq_fg.A1
    mean_bg = (mean - wt_fg*mean_fg)/wt_bg
    mean_sq_bg = (mean_sq - wt_fg*mean_sq_fg)/wt_bg
    var_fg = (mean_sq_fg - mean_fg**2)
    var_bg = (mean_sq_bg - mean_bg**2)
    return mean_fg, var_fg, mean_bg, var_bg



def gex_nbrhood_rank_tcr_scores(
        adata,
        nbrs_gex,
        tcr_score_names,
        pval_threshold,
        prefix_tag='nbr',
        min_num_fg=3,
        verbose=True,
        ttest_pval_threshold_for_mwu_calc=None
):
    ''' pvalues are bonferroni corrected (actually just multiplied by numtests)
    '''
    num_clones = adata.shape[0]

    tcrs = pp.retrieve_tcrs_from_adata(adata)
    pp.add_mait_info_to_adata_obs(adata)
    is_mait = adata.obs['is_mait']
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold*10

    print('making tcr score table:', tcr_score_names)
    score_table = tcr_scoring.make_tcr_score_table(adata, tcr_score_names)
    score_table_sq = np.multiply(score_table, score_table)
    mean = score_table.mean(axis=0)
    mean_sq = score_table_sq.mean(axis=0)

    num_nonempty_nbrhoods = sum(1 for x in nbrs_gex if len(x)>0)
    pval_rescale = num_nonempty_nbrhoods * len(tcr_score_names)

    nbrhood_mask = np.full( (num_clones,), False)

    results = []

    for ii in range(num_clones):
        if len(nbrs_gex[ii])==0:
            continue
        nbrhood_mask.fill(False)
        nbrhood_mask[ nbrs_gex[ii] ] = True
        nbrhood_mask[ ii ] = True
        mean_fg, var_fg, mean_bg, var_bg = get_split_mean_var(score_table, score_table_sq, nbrhood_mask, mean, mean_sq)
        num_fg = np.sum(nbrhood_mask)
        if num_fg < min_num_fg:
            continue
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_fg, std1=np.sqrt(np.maximum(var_fg, 1e-12)), nobs1=num_fg,
            mean2=mean_bg, std2=np.sqrt(np.maximum(var_bg, 1e-12)), nobs2=num_clones-num_fg,
            equal_var=False  # Welch's
        )

        scores[np.isnan(scores)] = 0
        pvals[np.isnan(pvals)] = 1

        # crude bonferroni
        pvals *= pval_rescale

        nbrhood_clusters_gex, nbrhood_clusters_tcr = None,None # lazy

        for ind in np.argsort(pvals):
            pval = pvals[ind]

            if pval>ttest_pval_threshold_for_mwu_calc:
                continue

            _,mwu_pval = mannwhitneyu( score_table[:,ind][nbrhood_mask], score_table[:,ind][~nbrhood_mask],
                                       alternative='two-sided')
            mwu_pval_adj = mwu_pval * pval_rescale

            if min(pval, mwu_pval_adj) <= pval_threshold:
                if nbrhood_clusters_gex is None: # lazy
                    nbrhood_clusters_gex = clusters_gex[nbrhood_mask]
                    nbrhood_clusters_tcr = clusters_tcr[nbrhood_mask]
                    nbrhood_is_mait = is_mait[nbrhood_mask]

                # get info about the clones most contributing to this skewed score
                score_name = tcr_score_names[ind]
                score = scores[ind] # ie the t-statistic

                num_top = max(1,num_fg//4)
                if score>0: # score is high
                    top_indices = np.argpartition( score_table[:,ind][nbrhood_mask], -num_top)[-num_top:]
                else: # score is low
                    top_indices = np.argpartition( score_table[:,ind][nbrhood_mask], num_top-1)[:num_top]


                gex_cluster = Counter( nbrhood_clusters_gex[ top_indices ]).most_common(1)[0][0]
                tcr_cluster = Counter( nbrhood_clusters_tcr[ top_indices ]).most_common(1)[0][0]
                mait_fraction = np.sum(nbrhood_is_mait[ top_indices ] )/len(top_indices)

                if verbose and mwu_pval_adj <= pval_threshold:
                    print('gex_{}_score: {:9.2e} {:9.2e} {:7.2f} clp {:2d} {:2d} {:7.3f} {:7.3f} {:15s} {} {} mf: {:.3f} {}'\
                          .format( prefix_tag, pval, mwu_pval_adj, score, gex_cluster, tcr_cluster, mean_fg[ind],
                                   mean_bg[ind], score_name,
                                   ' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3]), mait_fraction, ii ))

                results.append( dict(ttest_pvalue_adj=pval,
                                     ttest_stat=score,
                                     mwu_pvalue_adj=mwu_pval_adj,
                                     gex_cluster=gex_cluster,
                                     tcr_cluster=tcr_cluster,
                                     num_fg=num_fg,
                                     mean_fg=mean_fg[ind],
                                     mean_bg=mean_bg[ind],
                                     feature=score_name,
                                     mait_fraction=mait_fraction,
                                     clone_index=ii) )

        sys.stdout.flush()
    return pd.DataFrame(results)


def tcr_nbrhood_rank_genes_fast(
        adata,
        nbrs_tcr,
        pval_threshold,
        top_n=50,
        verbose=True,
        prefix_tag='nbr',
        min_num_fg=3,
        clone_display_names=None,
        ttest_pval_threshold_for_mwu_calc=None
):
    ''' Modeled on scanpy rank_genes_groups
    All pvals are crude bonferroni corrected for:
    * number of non-empty nbrhoods in nbrs_tcr and number of genes in adata.raw.X with at least 3 nonzero cells
      (see pval_rescale below)
    '''

    ## unpack from adata
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    tcrs = pp.retrieve_tcrs_from_adata(adata)
    pp.add_mait_info_to_adata_obs(adata)
    is_mait = adata.obs['is_mait']
    organism = adata.uns['organism']
    ## done unpacking ###############################

    if clone_display_names is None:
        clone_display_names = [ '{} {}'.format(' '.join(x[0][:3]), ' '.join(x[1][:3])) for x in tcrs ]

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold * 10

    assert pp.check_if_raw_matrix_is_logged(adata)

    rankby_abs = False # following scanpy: this means that we only look at enriched/upregulated/higher score values

    num_clones = adata.shape[0]

    genes = list(adata.raw.var_names)
    num_real_genes = len(genes)
    X = adata.raw.X
    assert sps.issparse(X)
    assert X.shape[1] == len(genes)

    X_csc = adata.raw.X.tocsc()

    X_sq = X.multiply(X)

    mean = X.mean(axis=0)
    mean_sq = X_sq.mean(axis=0)

    ## add some extra fake genes
    if 'nndists_gex' in adata.obs_keys():
        nndists_gex = np.array(adata.obs['nndists_gex'])
    else:
        print('WARNING nndists_gex not in adata.obs!')
        nndists_gex = np.zeros(num_clones)

    if len(set(adata.obs['clone_sizes']))==1: # only one clone_size
        genes2 = ['nndists_gex_rank']
        X2 = np.log1p(np.argsort(-1*nndists_gex))[:,np.newaxis]
    else:
        genes2 = ['clone_sizes', 'nndists_gex_rank']
        X2 = np.vstack( [np.log1p(np.array(adata.obs['clone_sizes'])),
                         np.log1p(np.argsort(-1*nndists_gex))] ).transpose()
    assert X2.shape == (num_clones,len(genes2))
    X2_sq = X2*X2
    mean2 = X2.mean(axis=0)
    mean2_sq = X2_sq.mean(axis=0)

    # need this?
    mean = mean.A1
    mean_sq = mean_sq.A1

    num_nonempty_nbrhoods = sum(1 for x in nbrs_tcr if len(x)>0)

    # len(genes) is probably too hard since lots of the genes are all zeros
    min_nonzero_cells = 3
    gene_nonzero_counts = Counter( X.nonzero()[1] )
    bad_gene_mask = np.array([ gene_nonzero_counts[x] < min_nonzero_cells for x in range(len(genes)) ]+
                             [False]*len(genes2))
    n_genes_eff = np.sum(~bad_gene_mask)
    pval_rescale = num_nonempty_nbrhoods * n_genes_eff

    results = []

    genes.extend(genes2) # since we are hstacking the vars, etc
    reference_indices = np.arange(len(genes), dtype=int)

    for ii in range(num_clones):
        if len(nbrs_tcr[ii])==0:
            continue
        nbrhood_mask = np.full( (num_clones,), False)
        nbrhood_mask[ nbrs_tcr[ii] ] = True
        nbrhood_mask[ ii ] = True

        mean_fg, var_fg, mean_bg, var_bg = get_split_mean_var(X, X_sq, nbrhood_mask, mean, mean_sq)
        mean2_fg, var2_fg, mean2_bg, var2_bg = get_split_mean_var(X2, X2_sq, nbrhood_mask, mean2, mean2_sq)

        num_fg = np.sum(nbrhood_mask)
        if num_fg < min_num_fg:
            continue
        mean_fg = np.hstack([mean_fg, mean2_fg])
        mean_bg = np.hstack([mean_bg, mean2_bg]) # note that we dont do the variances...
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_fg, std1=np.sqrt(np.maximum(np.hstack([var_fg, var2_fg]), 1e-12)), nobs1=num_fg,
            mean2=mean_bg, std2=np.sqrt(np.maximum(np.hstack([var_bg, var2_fg]), 1e-12)), nobs2=num_clones-num_fg,
            equal_var=False  # Welch's
        )

        # scanpy code:
        scores[np.isnan(scores)] = 0.
        pvals [np.isnan(pvals)] = 1.
        pvals [bad_gene_mask] = 1.
        logfoldchanges = np.log2((np.expm1(mean_fg) + 1e-9) / (np.expm1(mean_bg) + 1e-9))

        pvals_adj = pvals * pval_rescale

        scores_sort = np.abs(scores) if rankby_abs else scores
        partition = np.argpartition(scores_sort, -top_n)[-top_n:]
        partial_indices = np.argsort(scores_sort[partition])[::-1]
        global_indices = reference_indices[partition][partial_indices]

        nbrhood_clusters_gex, nbrhood_clusters_tcr = None,None

        for igene, ind in enumerate(global_indices):
            gene = genes[ind]
            if util.is_vdj_gene(gene, organism, include_constant_regions=True):
                continue
            pval_adj = pvals_adj[ind]
            log2fold= logfoldchanges[ind]

            if pval_adj > ttest_pval_threshold_for_mwu_calc:
                continue

            is_real_gene = ind < num_real_genes
            # here we are looking for genes (or clone_sizes/inverted nndists) that are LARGER in the forground (fg)
            if is_real_gene:
                col = X_csc[:,ind][nbrhood_mask]
                noncol = X_csc[:,ind][~nbrhood_mask]
                _, mwu_pval = mannwhitneyu( col.todense(), noncol.todense(), alternative='greater' )
            else:
                col = X2[:,ind-num_real_genes][nbrhood_mask]
                noncol = X2[:,ind-num_real_genes][~nbrhood_mask]
                _, mwu_pval = mannwhitneyu(col, noncol, alternative='greater')
            mwu_pval_adj = mwu_pval * pval_rescale

            if min(mwu_pval_adj, pval_adj) < pval_threshold:
                if nbrhood_clusters_gex is None: # lazy
                    nbrhood_clusters_gex = clusters_gex[nbrhood_mask]
                    nbrhood_clusters_tcr = clusters_tcr[nbrhood_mask]
                    nbrhood_is_mait = is_mait[nbrhood_mask]


                # better annotation of the enriched tcrs...
                num_top = num_fg//4
                if is_real_gene: # col is sparse...
                    if len(col.data)>num_top: # more than a quarter non-zero
                        top_indices = col.indices[ np.argpartition(col.data, -num_top)[-num_top:] ]
                        #bot_indices = col.indices[ np.argpartition(col.data, num_top-1)[:num_top] ]
                        #assert np.mean(col[top_indices]) > np.mean(col[bot_indices])
                    else:
                        top_indices = col.indices
                else:
                    top_indices = np.argpartition(col, -num_top)[-num_top:]

                #top_indices = np.nonzero(nbrhood_mask)[0][col_top_inds]

                gex_cluster = Counter( nbrhood_clusters_gex[ top_indices ]).most_common(1)[0][0]
                tcr_cluster = Counter( nbrhood_clusters_tcr[ top_indices ]).most_common(1)[0][0]
                mait_fraction = np.sum(nbrhood_is_mait[ top_indices ] )/len(top_indices)


                if verbose and mwu_pval_adj<=pval_threshold:
                    print('tcr_{}_gene: {:9.2e} {:9.2e} {:7.3f} clp {:2d} {:2d} {:8s} {:.4f} {:.4f} {:4d} {} mf: {:.3f} {} {}'\
                          .format( prefix_tag, pval_adj, mwu_pval_adj, log2fold, gex_cluster, tcr_cluster, gene,
                                   mean_fg[ind], mean_bg[ind], num_fg, clone_display_names[ii], mait_fraction,
                                   ii, igene ))
                results.append( dict(ttest_pvalue_adj=pval_adj,
                                     mwu_pvalue_adj=mwu_pval_adj,
                                     log2enr=log2fold,
                                     gex_cluster=gex_cluster,
                                     tcr_cluster=tcr_cluster,
                                     feature=gene,
                                     mean_fg=mean_fg[ind],
                                     mean_bg=mean_bg[ind],
                                     num_fg=num_fg,
                                     clone_index=ii,
                                     mait_fraction=mait_fraction ) )

        sys.stdout.flush()

    return pd.DataFrame(results)

def compute_distance_correlations( adata, verbose=False ):
    ''' return pvalues, rvalues  (each 1 1d numpy array of shape (num_clones,))
    '''
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    agroups, bgroups = pp.setup_tcr_groups(adata)

    print('compute D_gex', adata.shape[0])
    D_gex = pairwise_distances( adata.obsm['X_pca_gex'], metric='euclidean' )

    print('compute D_tcr', adata.shape[0])
    D_tcr = pairwise_distances( adata.obsm['X_pca_tcr'], metric='euclidean' )

    pval_rescale = adata.shape[0]
    print('compute distance correlations' )
    results = []
    for ii, (agroup, bgroup)  in enumerate( zip( agroups, bgroups ) ):
        if verbose and ii%1000==0:
            print(ii)
        mask = (agroups != agroup)&(bgroups != bgroup)
        reg = linregress( D_gex[ii, mask], D_tcr[ii, mask])
        pval, rval = reg.pvalue, reg.rvalue
        pval *= pval_rescale
        results.append( [ pval, rval ] )
        if verbose and pval < 1:
            print(f'distcorr: {pval:9.2e} {rval:7.3f} {clusters_gex[ii]:2d} {clusters_tcr[ii]:2d} {ii:4d}')

    results = np.array(results)
    return results[:,0], results[:,1]


def compute_cluster_interactions( aclusters_in, bclusters_in, barcodes_in, barcode2tcr, outlog, max_pval = 1.0 ):

    ''' Compute the cluster-cluster intxn (positive) with the lowest hypergeometric pval, and report.
    Then eliminate one of the two interacting clusters: first try taking the one with the fewest cells
    outside the pairwise interaction... Then iterate.

    pval includes a rescaling factor = num_a_clusters * num_b_clusters, max_pval applies to this rescaled pval

    This code is not currently being used in the conga pipeline but might be useful someday

    '''


    # condense a/b tcr groups within cluster-cluster intxns
    # confusing: the a/b here is alpha vs beta not aclusters vs bclusters (ie gex vs tcr)
    tcrs_list = [ barcode2tcr[x] for x in barcodes_in ]
    atcrs = sorted( set( x[0] for x in tcrs_list ) )
    btcrs = sorted( set( x[1] for x in tcrs_list ) )
    atcr2agroup = dict( (y,x) for x,y in enumerate(atcrs))
    btcr2bgroup = dict( (y,x) for x,y in enumerate(btcrs))
    groups1 = np.array( [ atcr2agroup[x[0]] for x in tcrs_list ] ) # alpha tcr groups
    groups2 = np.array( [ btcr2bgroup[x[1]] for x in tcrs_list ] ) # beta tcr groups

    aclusters = np.copy( aclusters_in )
    bclusters = np.copy( bclusters_in )
    barcodes = np.copy( barcodes_in )

    pval_rescale = ( np.max(aclusters)+1 ) * ( np.max(bclusters)+1 )

    while True:
        if aclusters.shape[0] == 0:
            break
        num_a = np.max(aclusters)+1
        num_b = np.max(bclusters)+1

        overlaps1 = {}
        overlaps2 = {}

        for ii in range(len(barcodes)):
            cl = (aclusters[ii],bclusters[ii])
            if cl not in overlaps1:
                overlaps1[cl] = set()
                overlaps2[cl] = set()
            overlaps1[cl].add( groups1[ii] )
            overlaps2[cl].add( groups2[ii] )

        abcounts = np.zeros( (num_a,num_b), dtype=int )
        for i in range(num_a):
            for j in range(num_b):
                if (i,j) in overlaps1:
                    abcounts[i,j] = min( len(overlaps1[(i,j)]), len(overlaps2[(i,j)] ) )

        acounts = np.sum( abcounts, axis=1 )
        bcounts = np.sum( abcounts, axis=0 )
        assert acounts.shape[0] == num_a
        assert bcounts.shape[0] == num_b


        pvals = np.ones( ( num_a, num_b ) )

        total = np.sum( acounts ) # or bcounts

        for a in range(num_a):
            for b in range(num_b):
                overlap = abcounts[a,b]
                expected = acounts[a] * bcounts[b] / float(total)

                if overlap>expected and overlap>1: # otherwise just stays at 1.0
                    pvals[a,b] = hypergeom.sf( overlap-1, total, acounts[a], bcounts[b] )


        # most significant interaction:
        a,b = np.unravel_index( np.argmin( pvals ), pvals.shape )
        overlap = abcounts[a,b]
        if not overlap:
            break
        overlap_barcodes = barcodes[ ( aclusters==a ) & ( bclusters==b ) ]
        expected = acounts[a] * bcounts[b] / float(total)
        pval = pval_rescale * pvals[a,b]

        if pval > max_pval:
            break

        outlog.write('clusclus2_intxn: {:2d} {:2d} {:8.1e} {:3d} {:6.1f} {:4d} {:4d} {:2d} {:2d} {:5d} ex {}\n'\
                     .format( a,b,pval,overlap,expected,acounts[a],bcounts[b],
                              len(set(aclusters)),len(set(bclusters)),len(barcodes),
                              len(overlap_barcodes)-overlap ))

        # now iterate
        if acounts[a] <= bcounts[b]:
            mask = aclusters!=a
        else:
            mask = bclusters!=b
        aclusters = aclusters[ mask ]
        bclusters = bclusters[ mask ]
        barcodes = barcodes[ mask ]
        groups1 = groups1[ mask ]
        groups2 = groups2[ mask ]
    return



def setup_fake_nbrs_from_clusters_for_graph_vs_features_analysis( clusters ):
    ''' Make a fake nbr graph in which one clone in each cluster has a set of nbrs to the other
    cluster members. Everybody else has empty nbr lists. For graph-vs-feature correlation analyses.
    '''
    clusters = np.array(clusters) # just in case; there's a problem w/ np.nonzero on some pandas series

    fake_nbrs = []
    seen = set()
    for cl in clusters:
        if cl in seen:
            fake_nbrs.append([])
        else:
            seen.add(cl)
            fake_nbrs.append(np.nonzero( clusters==cl )[0])
    return fake_nbrs

def find_hotspot_features(
        X,
        nbrs,
        features,
        pval_threshold,
        verbose=False
):
    """ My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """
    print('START computing H matrix', type(X)) # we want this to be a sps.csr_matrix since we are working with rows...
    assert type(X) is sps.csr_matrix # right now anyhow; not a strict requirement

    num_clones, num_features = X.shape
    assert len(features) == num_features

    # compute mean and sdev of raw gene expression matrix
    X_sq = X.multiply(X)

    X_mean = X.mean(axis=0).A1
    X_sq_mean = X_sq.mean(axis=0).A1
    X_mean_sq = X_mean**2
    X_var = X_sq_mean - X_mean_sq

    H = sps.csr_matrix( np.zeros((num_features,)) )
    indegrees = np.zeros((num_clones,))

    last_time = time.time()
    for ii in range(num_clones):
        if ii%500==0:
            elapsed_time = time.time() - last_time
            if ii:
                rate = 10*elapsed_time
            else:
                rate = 0
            print(f'computing H matrix {ii:6d} {num_clones:6d} {rate:12.6f}')
            sys.stdout.flush()
            last_time = time.time()
        X_ii = X[ii,:]
        ii_nbrs = nbrs[ii]
        if len(ii_nbrs)==0:
            continue
        indegrees[ii_nbrs] += 1
        H += X_ii.multiply( X[ii_nbrs,:].sum(axis=0) )

    # multiply the indegree matrix by the raw X matrix by the means
    indegrees_mat = sps.csr_matrix( indegrees[:, np.newaxis] )
    assert indegrees_mat.shape == (num_clones,1)

    X_mean_mat = sps.csr_matrix( X_mean )
    assert X_mean_mat.shape == (1, num_features)

    Y = X.multiply( indegrees_mat ).multiply( X_mean_mat ).sum(axis=0)

    assert H.shape == Y.shape

    # Y is a numpy.matrix
    #
    H = (H-Y).A1
    mask1 = (H==0)
    mask2 = (X_var==0)
    mask3 = mask2 & (~mask1) # H nonzero but stddev 0
    print('zeros:', np.sum(mask1), np.sum(mask2), np.sum(mask3))

    H /= np.maximum(1e-9, X_var)

    H[X_var==0] = 0 # set to zero if stddev is 0

    inds = np.argsort(H)[::-1] # decreasing

    # the simple estimate for the variance of H is the total number of neighbors
    # compute the variance
    nbrs_sets = []
    for ii in range(num_clones):
        nbrs_sets.append( frozenset(nbrs[ii]) )

    H_var = 0
    print('compute H_var')
    for ii in range(num_clones):
        for jj in nbrs[ii]:
            if ii in nbrs_sets[jj]:
                H_var += 2
            else:
                H_var += 1

    H /= np.sqrt(H_var)

    results = []
    for ind in inds:
        feature = features[ind]
        true_std = np.std(X[:,ind].toarray()[:,0])
        if true_std<1e-6:
            print('WHOAH var prob? {} {} =?= {}'.format(feature, X_var[ind], true_std**2))
            continue
        Z = H[ind]
        pvalue_adj = num_features * norm.sf(Z)
        if pvalue_adj > pval_threshold:
            break
        if verbose:
            print('top_var: {} =?= {} {} {} {} {}'.format(X_var[ind], true_std**2, true_std, Z, pvalue_adj, feature))
        results.append(OrderedDict(Z=Z, pvalue_adj=pvalue_adj, feature=feature))

    return pd.DataFrame(results)




def find_hotspot_genes(
        adata,
        nbrs_tcr,
        pval_threshold
):
    """ My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    clusters_gex = np.array(adata.obs['clusters_gex'])
    X = adata.raw.X
    genes = adata.raw.var_names

    num_clones = adata.shape[0]
    num_clusters = np.max(clusters_gex)+1
    Y = np.zeros((num_clones, num_clusters))
    for ii in range(num_clusters):
        Y[:,ii] = (clusters_gex==ii).astype(float)
    Y = sps.csr_matrix(Y)

    print('stacking extra columns!', X.shape, Y.shape)
    X = sps.hstack([X,Y]).tocsr()
    print('DONE stacking extra columns!')


    df = find_hotspot_features(X, nbrs_tcr, list(genes)+['gex_cluster{}'.format(x) for x in range(num_clusters)],
                               pval_threshold)

    if df.shape[0]==0:
        return df

    # filter out VDJ genes
    mask = [ not util.is_vdj_gene(x.feature, organism, include_constant_regions=True) for x in df.itertuples()]

    df = df[mask]
    df['feature_type'] = 'gex'

    # show the top 100 hits
    for ii,l in enumerate(df[:100].itertuples()):
        print('hotspot_gene: {:4d} {:9.3f} {:8.1e} {:10s}'.format(ii, l.Z, l.pvalue_adj, l.feature))

    return df


def find_hotspot_tcr_features(
        adata,
        nbrs_gex,
        pval_threshold,
        min_gene_count=5
):
    """ My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    tcrs = pp.retrieve_tcrs_from_adata(adata)
    num_clusters = np.max(adata.obs['clusters_tcr'])+1

    organism_genes = all_genes[organism]
    counts = Counter( [ organism_genes[x[i_ab][j_vj]].count_rep
                        for x in tcrs for i_ab in range(2) for j_vj in range(2)] )
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    features = tcr_scoring.all_tcr_scorenames + count_reps + ['tcr_cluster{}'.format(x) for x in range(num_clusters)]
    score_table = tcr_scoring.make_tcr_score_table(adata, features)
    X = sps.csr_matrix(score_table)
    #print('find_hotspot_tcr_features: X=', X)

    df = find_hotspot_features(X, nbrs_gex, features, pval_threshold)#, verbose=True)
    df['feature_type'] = 'tcr'

    # show the top 100 hits
    for ii,l in enumerate(df[:100].itertuples()):
        print('hotspot_tcr_feature: {:4d} {:9.3f} {:8.1e} {:10s}'.format(ii, l.Z, l.pvalue_adj, l.feature))

    return df


def find_hotspot_nbrhoods(
        adata,
        nbrs_gex,
        nbrs_tcr,
        pval_threshold,
        also_use_cluster_graphs=False
):
    """ My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    tcrs = pp.retrieve_tcrs_from_adata(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    nbrs_gex_clusters = setup_fake_nbrs_from_clusters_for_graph_vs_features_analysis( clusters_gex )
    nbrs_tcr_clusters = setup_fake_nbrs_from_clusters_for_graph_vs_features_analysis( clusters_tcr )
    # create csr_matrix
    # >>> row = np.array([0, 0, 1, 2, 2, 2])
    # >>> col = np.array([0, 2, 2, 0, 1, 2])
    # >>> data = np.array([1, 2, 3, 4, 5, 6])
    # >>> csr_matrix((data, (row, col)), shape=(3, 3)).toarray()
    # array([[1, 0, 2],
    #        [0, 0, 3],
    #        [4, 5, 6]])

    num_clones = adata.shape[0]

    dfl = []
    for feature_nbr_tag, feature_nbrs, graph_tag, graph_nbrs in [ ['gex', nbrs_gex, 'graph', nbrs_tcr],
                                                                  ['gex', nbrs_gex, 'clust', nbrs_tcr_clusters],
                                                                  ['tcr', nbrs_tcr, 'graph', nbrs_gex],
                                                                  ['tcr', nbrs_tcr, 'clust', nbrs_gex_clusters]]:
        if graph_tag == 'clust' and not also_use_cluster_graphs:
            continue

        rows = [] # list of (1,num_clones) csr_matrix'es
        for ii in range(num_clones):
            ii_nbrs = feature_nbrs[ii]
            num_nbrs = len(ii_nbrs)
            data = np.full((num_nbrs,),1.0)
            row_ind = np.full((num_nbrs,),0).astype(int)
            new_row = sps.csr_matrix((data, (row_ind, ii_nbrs)), shape=(1,num_clones))
            new_row[0,ii] = 1.
            rows.append(new_row)
        X = sps.vstack(rows)
        assert X.shape == (num_clones, num_clones)
        # find_hotspot_features expects X to look like the GEX matrix, ie
        #  shape=(num_clones,num_features) ie the features are the columns
        # but we set it up with the features (the nbrhoods) as the rows
        # so we need to transpose
        X = X.transpose().tocsr()

        feature_names = [str(x) for x in range(num_clones)]
        #
        df = find_hotspot_features(X, graph_nbrs, feature_names, pval_threshold)#, verbose=True)

        if df.shape[0] == 0:
            continue

        feature_type = '{}_nbrs_vs_{}'.format(feature_nbr_tag, graph_tag)
        df['feature_type'] = feature_type
        df['clone_index'] = df.feature.astype(int)

        # show the top 100 hits
        for ii,l in enumerate(df[:100].itertuples()):
            atcr, btcr = tcrs[l.clone_index]
            print('hotspot_{}: {:4d} {:9.3f} {:8.1e} {:2d} {:2d} {:4d} {} {} {} {} {} {}'\
                  .format(feature_type, ii, l.Z, l.pvalue_adj,
                          clusters_gex[l.clone_index], clusters_tcr[l.clone_index],
                          l.clone_index,
                          atcr[0], atcr[1], atcr[2], btcr[0], btcr[1], btcr[2]))
        dfl.append(df)

    if dfl:
        df = pd.concat(dfl)
        return df.drop(columns=['feature']) # feature is just str(clone_index), ie superfluous
    else:
        return pd.DataFrame()

def find_hotspots( adata, nbrs, pval_threshold = None ):
    """
    wrapper function combining Hotspot calculations for TCR and GEX. Returns combined df.
    nbrs: tuple of (nbrs_gex,nbrs_tcr) from preprocess.calc_nbrs used to look for correlations
    pval_threshold: pval for Bonferroni test. default is 0.05
    """

    if pval_threshold is None:
        pval_threshold = 0.05
    print(f'conga.correlations.find_hotspots:: Using pval_threshold = {pval_threshold}')

    nbrs_gex, nbrs_tcr = nbrs

    gene_df = find_hotspot_genes(adata , nbrs_tcr, pval_threshold )
    gene_df['feature_type'] = 'gex'

    tcr_df = find_hotspot_tcr_features(adata , nbrs_gex , pval_threshold  )
    tcr_df['feature_type'] = 'tcr'

    hotspot_df = pd.concat([gene_df, tcr_df])

    return hotspot_df

def find_batch_biases(
        adata,
        all_nbrs,
        pval_threshold=0.05,
        exclude_batch_keys=None
):
    ''' Look for graph neighborhoods that are biased in their batch distribution

    Returns a tuple of two dataframes: nbrhood_results, hotspot_results

    in nbrhood_results, the hits are grouped by clusters_tcr assignment
    and shared biased batches into
    groups using single-linkage clustering, stored in the 'cluster_group' column

    '''
    if 'batch_keys' not in adata.uns_keys():
        print('find_batch_biases:: no batch_keys in adata.uns!!!')
        return

    tcrs = pp.retrieve_tcrs_from_adata(adata)

    # for grouping the hit clones
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    num_clones = adata.shape[0]
    batch_keys = adata.uns['batch_keys']

    nbrhood_results = []

    hotspot_results = []

    is_significant = np.full((num_clones,), False)

    for bkey in batch_keys:
        if exclude_batch_keys and bkey in exclude_batch_keys:
            continue
        bcounts = np.array(adata.obsm[bkey])
        num_choices = bcounts.shape[1]
        assert bcounts.shape == (num_clones, num_choices)
        clone_sizes = np.sum(bcounts,axis=1)
        bfreqs = bcounts.astype(float)/clone_sizes[:,np.newaxis]

        ## hotspot analysis
        X = sps.csr_matrix(bfreqs)

        for nbr_frac in all_nbrs:
            nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
            for nbrs_tag, nbrs in [['gex', nbrs_gex], ['tcr', nbrs_tcr]]:
                if nbrs_tag=='gex':
                    continue
                ## look for neighborhoods with skewed distribution of scores
                num_nbrs = nbrs.shape[1]
                bfreqs_nbr_avged = bfreqs[nbrs[:, :, np.newaxis],
                                          np.arange(num_choices)[np.newaxis, np.newaxis, :]]
                assert bfreqs_nbr_avged.shape == (num_clones, num_nbrs, num_choices)
                bfreqs_nbr_avged = (bfreqs + bfreqs_nbr_avged.sum(axis=1))/(num_nbrs+1)
                assert bfreqs_nbr_avged.shape == (num_clones, num_choices)

                bfreqs_mean = np.mean(bfreqs, axis=0)
                bfreqs_std = np.std(bfreqs, axis=0)

                zscores = (bfreqs_nbr_avged - bfreqs_mean[np.newaxis,:])/bfreqs_std[np.newaxis,:]

                for ib in range(num_choices):
                    if bfreqs_std[ib] < 1e-6:
                        continue # no variation at this choice
                    for ii in np.argsort(-1*zscores[:,ib]):
                        zscore = zscores[ii,ib]
                        if zscore<1e-2: # .1 could be significant but not really .01
                            break
                        nbr_scores = list(bfreqs[:,ib][nbrs[ii]])+[bfreqs[ii,ib]]
                        nbrs_mask = np.full((num_clones,), False)
                        nbrs_mask[nbrs[ii]] = True
                        nbrs_mask[ii] = True
                        non_nbr_scores = bfreqs[:,ib][~nbrs_mask]
                        _,mwu_pval1 = mannwhitneyu( nbr_scores, non_nbr_scores, alternative='greater')
                        #_,mwu_pval2 = mannwhitneyu( nbr_scores, non_nbr_scores, alternative='less')
                        mwu_pval1 *= num_clones
                        if mwu_pval1 < pval_threshold:
                            is_significant[ii] = True
                            nbrhood_results.append( OrderedDict(
                                batch_key=bkey,
                                batch_choice=ib,
                                nbrs_tag=nbrs_tag, # 'gex' or 'tcr'
                                nbr_frac=nbr_frac,
                                clone_index=ii,
                                pvalue_adj=mwu_pval1,
                                zscore=zscore,
                                nbrs_mean=np.mean(nbr_scores),
                                non_nbrs_mean=np.mean(non_nbr_scores)
                                ))

                            print('nbr_batch_bias: {:9.1e} {:7.3f} {:7.4f} {:7.4f} {} {} {} {} {} {}'\
                                  .format(mwu_pval1, zscore, np.mean(nbr_scores), np.mean(non_nbr_scores),
                                          bkey, ib, nbr_frac, nbrs_tag,
                                          ' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3])))

                ## hotspot
                features = ['{}_{}'.format(bkey, x) for x in range(num_choices)]
                print(f'find hotspot features for {bkey} vs {nbrs_tag} at {nbr_frac:.4f}')
                results = find_hotspot_features(X, nbrs, features, pval_threshold)

                for ii, l in enumerate(results.itertuples()):
                    if l.pvalue_adj <= pval_threshold:
                        hotspot_results.append( OrderedDict(
                            batch_key='_'.join(l.feature.split('_')[:-1]),
                            batch_choice=int(l.feature.split('_')[-1]),
                            nbrs_tag=nbrs_tag,
                            nbr_frac=nbr_frac,
                            pvalue_adj=l.pvalue_adj,
                            zscore=l.Z,
                            ))

                    print('hotspot_batch: {:2d} {:9.3f} {:8.1e} {} {} {:.4f}'\
                          .format(ii, l.Z, l.pvalue_adj, l.feature, nbrs_tag, nbr_frac))
    sys.stdout.flush()

    # now try identifying groups of related nbrhoods
    nbrhood_results = pd.DataFrame(nbrhood_results)
    if nbrhood_results.shape[0]:
        #
        significant_inds = np.nonzero(is_significant)[0]

        all_sig_nbrs = {}
        all_batches = {}
        for ii in significant_inds:
            all_sig_nbrs[ii] = set([ii])
            # the batch-key combos for which ii has a significant bias
            # itertuples with index=False and name=None returns simple row tuples
            all_batches[ii] = set( nbrhood_results[nbrhood_results.clone_index==ii][['batch_key','batch_choice']]\
                                   .itertuples(index=False, name=None))

        for ii in significant_inds:
            # for nbrs_tag, nbr_frac in set(nbrhood_results[nbrhood_results.clone_index==ii][['nbrs_tag','nbr_frac']]\
            #                               .itertuples(index=False, name=None)):
            #     nbrs = all_nbrs[nbr_frac][0] if nbrs_tag=='gex' else all_nbrs[nbr_frac][1]
            for nbrs_tag in set(nbrhood_results[nbrhood_results.clone_index==ii].nbrs_tag): # should just be 'tcr' now
                assert nbrs_tag=='tcr' # tmp hack
                clusters = clusters_gex if nbrs_tag=='gex' else clusters_tcr
                for jj in np.nonzero( (clusters==clusters[ii]) & is_significant )[0]:
                    if len(all_batches[ii]&all_batches[jj]):
                        all_sig_nbrs[ii].add(jj)
                        all_sig_nbrs[jj].add(ii)


        all_smallest_nbr = {}
        for ii in significant_inds:
            all_smallest_nbr[ii] = min(all_sig_nbrs[ii])

        while True:
            updated = False
            for ii in significant_inds:
                nbr = all_smallest_nbr[ii]
                new_nbr = min(nbr, np.min([all_smallest_nbr[x] for x in all_sig_nbrs[ii]]))
                if nbr != new_nbr:
                    all_smallest_nbr[ii] = new_nbr
                    updated = True
            if not updated:
                break
        # define clusters, choose cluster centers
        clusters = np.array([0]*num_clones) # 0 if not clumped

        cluster_number=0
        for ii in significant_inds:
            nbr = all_smallest_nbr[ii]
            if ii==nbr:
                cluster_number += 1
                members = [ x for x,y in all_smallest_nbr.items() if y==nbr]
                clusters[members] = cluster_number

        for ii, nbrs in all_sig_nbrs.items():
            for nbr in nbrs:
                assert clusters[ii] == clusters[nbr] # confirm single-linkage clusters

        assert not np.any(clusters[is_significant]==0)
        assert np.all(clusters[~is_significant]==0)

        nbrhood_results['cluster_group'] = [ clusters[x.clone_index] for x in nbrhood_results.itertuples()]


    return nbrhood_results, pd.DataFrame(hotspot_results) # first is already converted to DF
