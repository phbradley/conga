################################################################################
import numpy as np
from scipy import stats
from sklearn.metrics import pairwise_distances
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind
#from scipy.sparse import issparse, csr_matrix
import scipy.sparse as sps
from statsmodels.stats.multitest import multipletests
from collections import Counter, OrderedDict
import scanpy as sc
from . import preprocess
from . import tcr_scoring
from . import util
from .tcrdist.all_genes import all_genes
from .tags import *
import sys
import pandas as pd
from sys import exit
import time #debugging

MIN_CONGA_SCORE = 1e-100

def _find_neighbor_neighbor_interactions(
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
    ''' This is a helper function used in graph-vs-graph analysis

    it runs the GEX KNN graph vs TCR KNN graph comparison

    Returns a pandas dataframe with results
    of all clones with nbr-nbr overlaps having pvalues <= pval_threshold

    AND a numpy array of the adjusted_pvals

    '''
    preprocess.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_invariant']

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
            # this is heuristic
            overlap_corrected = min(len(set(agroups[double_nbrs])),
                                    len(set(bgroups[double_nbrs])) )

            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                nbr_pval = pval_rescale * hypergeom.sf(
                    overlap_corrected-1, actual_num_clones,
                    num_neighbors_gex-delta, num_neighbors_tcr-delta)
                adjusted_pvalues[ii] = nbr_pval # update
                if nbr_pval > pval_threshold:
                    continue ## NOTE

        nbr_pval = max(nbr_pval, MIN_CONGA_SCORE) # no 0s
        results.append(dict(conga_score=nbr_pval,
                            num_neighbors_gex=num_neighbors_gex,
                            num_neighbors_tcr=num_neighbors_tcr,
                            overlap=overlap,
                            overlap_corrected=overlap_corrected,
                            mait_fraction=np.sum(is_mait[double_nbrs])/overlap,
                            clone_index=ii ))#double_nbrs ] )

    adjusted_pvalues = np.maximum(np.array(adjusted_pvalues), MIN_CONGA_SCORE)

    return pd.DataFrame(results), adjusted_pvalues

def _find_neighbor_cluster_interactions(
        adata,
        nbrs,
        clusters, # for example, or could be swapped: tcr/gex
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number of them
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True,
):
    ''' This is a helper function used in graph-vs-graph analysis

    It computes KNN graph vs cluster graph overlaps

    Returns a pandas dataframe with results
    of all clones with nbr-cluster overlaps having pvalues <= pval_threshold

    AND a numpy array of the adjusted_pvals

    '''
    preprocess.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_invariant']

    num_clones = len(nbrs)

    pval_rescale = num_clones if scale_pvals_by_num_clones else 1.0

    results = []
    adjusted_pvalues = [pval_rescale]*num_clones # initialize to big pval

    for ii in range(num_clones):
        is_ii_group = (agroups==agroups[ii]) | (bgroups==bgroups[ii])
        assert is_ii_group[ii]
        actual_num_clones = num_clones - np.sum(is_ii_group) - counts_correction

        ii_nbrs = nbrs[ii]
        ii_cluster = clusters[ii]
        ii_cluster_clustersize = (np.sum(clusters==ii_cluster) -
                                  np.sum((clusters==ii_cluster)&(is_ii_group)))

        num_neighbors = ii_nbrs.shape[0]
        overlap = np.sum( clusters[ii_nbrs] == ii_cluster )
        expected = float(ii_cluster_clustersize*num_neighbors)/actual_num_clones
        if overlap and overlap>expected:
            nbr_pval = pval_rescale * hypergeom.sf(
                overlap-1, actual_num_clones, num_neighbors,
                ii_cluster_clustersize )
        else:
            nbr_pval = pval_rescale

        adjusted_pvalues[ii] = nbr_pval
        if nbr_pval > pval_threshold:
            continue ## NOTE

        same_cluster_nbrs = ii_nbrs[ clusters[ii_nbrs] == ii_cluster ]
        assert len(same_cluster_nbrs)== overlap

        overlap_corrected = overlap
        if correct_overlaps_for_groups:
            overlap_corrected = min(len(set(agroups[same_cluster_nbrs])),
                                    len(set(bgroups[same_cluster_nbrs])))
            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                nbr_pval = pval_rescale * hypergeom.sf(
                    overlap_corrected-1, actual_num_clones, num_neighbors-delta,
                    ii_cluster_clustersize-delta )
                adjusted_pvalues[ii] = nbr_pval
                if nbr_pval > pval_threshold:
                    continue

        mait_fraction=np.sum(is_mait[same_cluster_nbrs])/overlap
        nbr_pval = max(nbr_pval, MIN_CONGA_SCORE) # no 0s
        results.append(dict(conga_score=nbr_pval,
                            num_neighbors=num_neighbors,
                            cluster_size=ii_cluster_clustersize,
                            overlap=overlap,
                            overlap_corrected=overlap_corrected,
                            mait_fraction=mait_fraction,
                            clone_index=ii ))#double_nbrs ] )

    adjusted_pvalues = np.maximum(np.array(adjusted_pvalues), MIN_CONGA_SCORE)

    return pd.DataFrame(results), adjusted_pvalues



def check_nbr_graphs_indegree_bias(all_nbrs):
    ''' this routine looks at bias in the number of neighbor edges
    going *into* each vertex (clonotype). By definition, each vertex
    will have the same number going out, but not necessarily the same
    number going in. This is especially true of the GEX graph.

    Generally there is less bias in the TCR graph. If both graphs were biased
    in the same direction (toward the same nodes), this could create
    spurious graph-vs-graph signal.

    So we look here at the correlation between the two biases.

    all_nbrs is a dict mapping from nbr_frac to (gex_nbrs, tcr_nbrs)
       setup by preprocess.calc_nbrs
    '''
    for nbr_frac in all_nbrs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        num_clones = nbrs_gex.shape[0]

        # look at the in-degree distribution
        expected_indegree = nbrs_gex.shape[1]
        gex_counts = Counter(nbrs_gex.flatten())
        gex_indegree_bias = np.array( [ gex_counts[x]/expected_indegree
                                        for x in range(num_clones) ] )
        print(f'gex_indegree_bias: nbr_frac= {nbr_frac:.4f}',
              stats.describe(gex_indegree_bias))

        tcr_counts = Counter(nbrs_tcr.flatten())
        tcr_indegree_bias = np.array( [ tcr_counts[x]/expected_indegree
                                        for x in range(num_clones) ] )
        print(f'tcr_indegree_bias: nbr_frac= {nbr_frac:.4f}',
              stats.describe(tcr_indegree_bias))

        # any correlation?
        # if there is strong correlation, this could skew the
        #  graph-vs-graph analysis
        print(f'indegree_bias_correlation: nbr_frac= {nbr_frac:.4f}',
              stats.linregress(gex_indegree_bias, tcr_indegree_bias))



def _make_csr_nbrs(nbrs):
    row = []
    for i, inbrs in enumerate(nbrs):
        row.extend([i]*len(inbrs))
    col = np.hstack(nbrs)
    assert len(row) == len(col)
    data = np.full((len(col),), 1)
    return sps.csr_matrix((data, (row, col)), shape=(len(nbrs), len(nbrs)))

def _compute_nbr_overlap_slow(gex_nbrs, tcr_nbrs):
    overlap = 0
    for g_nbrs, t_nbrs in zip(gex_nbrs, tcr_nbrs):
        g_nbrs_set = frozenset(g_nbrs)
        overlap += sum(x in g_nbrs_set for x in t_nbrs)
    return overlap

def _compute_graph_overlap_stats(
        gex_nbrs,
        tcr_nbrs,
        num_random_repeats,
        verbose=False,
        swaptags=False,
        max_calculation_time=2000,# in seconds
):
    ''' Helper function for graph-graph overlap summary analysis

    see compute_graph_vs_graph_stats(...) function below
    '''
    starttime = time.time()
    gtag,ttag = 'gex','tcr'
    if swaptags:
        gtag, ttag = ttag, gtag
    N = len(gex_nbrs)
    assert N == len(tcr_nbrs)

    ## if this will take too long, we may just estimate things
    gex_edges = sum(len(x) for x in gex_nbrs)
    tcr_edges = sum(len(x) for x in tcr_nbrs)
    expected_overlap = (sum(len(x)*len(y) for x,y in zip(gex_nbrs, tcr_nbrs))/
                        (N-1))

    # compute the bias in the number of incoming edges in the gex graph
    expected_indegree = gex_edges/N
    gex_counts = Counter()
    for nbrs in gex_nbrs:
        gex_counts.update(nbrs)
    gex_indegree_bias = np.array([gex_counts[x]/expected_indegree
                                  for x in range(N)])
    gex_indegree_bias_stats = stats.describe(gex_indegree_bias)
    if verbose:
        print('gex_indegree_bias:', gex_indegree_bias_stats)

    # compute the bias in the number of incoming edges in the tcr graph
    expected_indegree = tcr_edges/N
    tcr_counts = Counter()
    for nbrs in tcr_nbrs:
        tcr_counts.update(nbrs)
    tcr_indegree_bias = np.array([tcr_counts[x]/expected_indegree
                                  for x in range(N)])
    tcr_indegree_bias_stats = stats.describe(tcr_indegree_bias)
    if verbose:
        print('tcr_indegree_bias:', tcr_indegree_bias_stats)
    indegree_correlation = linregress(gex_indegree_bias, tcr_indegree_bias)

    # this is a little silly: it's basically C*gex_edges * tcr_edges/nodes**2
    # from smf.ols(f'log10_calculation_time ~ log10_nodes + log10_gex_edges + log10_tcr_edges'
    # Intercept         -4.453273
    # log10_nodes       -1.914652
    # log10_gex_edges    0.958948
    # log10_tcr_edges    1.044386
    # this was fitted with num_random_repeats = 100
    estimated_log10_calculation_time = (-1.9147 * np.log10(N)
                                        +0.9589 * np.log10(gex_edges)
                                        +1.0444 * np.log10(tcr_edges)
                                        -4.4533)
    estimated_calculation_time = (10**estimated_log10_calculation_time*
                                  num_random_repeats/100.)
    if estimated_calculation_time <= max_calculation_time:
        M0 = _make_csr_nbrs(gex_nbrs)
        M1 = _make_csr_nbrs(tcr_nbrs).tocoo()
        M2 = M1.copy()
        overlaps = []
        for r in range(num_random_repeats+1):
            if r:
                p = np.random.permutation(N)
            else:
                p = np.arange(N)
            M1.row = p[M2.row]
            M1.col = p[M2.col]
            overlap = M0.multiply(M1.tocsr()).sum()
            if verbose and r%10==0:
                print(f'{r:2d} {overlap:6d}')
            overlaps.append(overlap)
        o0 = overlaps[0]
        m,s = np.mean(overlaps[1:]), np.std(overlaps[1:])
        zscore_source = 'shuffling'
    else:
        zscore_source = 'fitting'
        o0 = _compute_nbr_overlap_slow(gex_nbrs, tcr_nbrs)

    ## params for log10_s determined with
    ## statsmodels.formula.api.ols(f'log10_overlap_sdev ~
    ##     log10_expected_overlap + total_log10_indegree_variance,...)
    # Intercept                       -0.340085
    # log10_expected_overlap           0.691433
    # total_log10_indegree_variance    0.253497
    total_log10_indegree_variance = (
        np.log10(gex_indegree_bias_stats.variance)+
        np.log10(tcr_indegree_bias_stats.variance))
    log10_s_fitted = (0.691433 * np.log10(expected_overlap)
                      +0.253497 * total_log10_indegree_variance
                      -0.340085)
    s_fitted = 10**log10_s_fitted

    if zscore_source=='fitting':
        s = s_fitted
        m = expected_overlap
        z_fitted = (o0-m)/s
        z = z_fitted
    else:
        z = (o0-m)/s
        z_fitted = (o0-expected_overlap)/s_fitted

    total_seconds = time.time() - starttime

    return {
        'overlap':o0,
        'expected_overlap':expected_overlap,
        'overlap_mean':m,
        'overlap_sdev':s,
        'overlap_zscore':z,
        'overlap_zscore_fitted':z_fitted,
        'overlap_zscore_source':zscore_source,
        'nodes':N,
        'calculation_time':total_seconds,
        'calculation_time_fitted':10**estimated_log10_calculation_time,
        f'{gtag}_edges':gex_edges,
        f'{ttag}_edges':tcr_edges,
        f'{gtag}_indegree_variance':gex_indegree_bias_stats.variance,
        f'{gtag}_indegree_skewness':gex_indegree_bias_stats.skewness,
        f'{gtag}_indegree_kurtosis':gex_indegree_bias_stats.kurtosis,
        f'{ttag}_indegree_variance':tcr_indegree_bias_stats.variance,
        f'{ttag}_indegree_skewness':tcr_indegree_bias_stats.skewness,
        f'{ttag}_indegree_kurtosis':tcr_indegree_bias_stats.kurtosis,
        'indegree_correlation_R':indegree_correlation.rvalue,
        'indegree_correlation_P':indegree_correlation.pvalue,
    }

# def _compute_graph_overlap_stats_old(
#         gex_nbrs_array,
#         tcr_nbrs,
#         num_random_repeats,
#         verbose = True
# ):
#     '''
#     I'm not sure about memory/compute efficiency here, for big datasets
#     '''
#     N = len(gex_nbrs_array)
#     gex_nbrs = [frozenset(x) for x in gex_nbrs_array]
#     overlaps = []
#     for r in range(num_random_repeats+1):
#         if r:
#             tcr_shuffle = np.random.permutation(N)
#         else:
#             tcr_shuffle = np.arange(N)
#         overlap = _compute_graph_overlap(gex_nbrs, tcr_nbrs, tcr_shuffle)
#         if verbose and r%10==0:
#             print(f'compute_graph_overlap_stats: rep {r:2d} {overlap:6d}')
#         overlaps.append(overlap)
#     o0 = overlaps[0]
#     m,s = np.mean(overlaps[1:]), np.std(overlaps[1:])
#     z = (o0-m)/s
#     gex_edges = sum(len(x) for x in gex_nbrs)
#     tcr_edges = sum(len(x) for x in tcr_nbrs)
#     return {
#         'overlap':o0,
#         'mean':m,
#         'sdev':s,
#         'zscore':z,
#         'nodes':N,
#         'gex_edges':gex_edges,
#         'tcr_edges':tcr_edges,
#     }

def compute_graph_vs_graph_stats(
        adata,
        all_nbrs,
        num_random_repeats = 100,
        outfile_prefix = None,
):
    '''Here we are assessing overall graph-vs-graph correlation by looking at
    the shared edges between TCR and GEX neighbor graphs and comparing
    that observed number to the number we would expect if the graphs were
    completely uncorrelated. Our null model for uncorrelated graphs is to
    take the vertices of one graph and randomly renumber them (permute their
    labels). We compare the observed overlap to that expected under this null
    model by computing a Z-score, either by permuting one of the graph's
    vertices many times to get a mean and standard deviation of the overlap
    distribution, or, for large graphs where this is time consuming,
    by using a regression model for the
    standard deviation.

    This is different from graph-vs-graph analysis which looks at graph
    overlap on a node-by-node basis

    returns a pandas dataframe with the results, and also stores them in
    adata.uns['conga_results'][conga.tags.GRAPH_VS_GRAPH_STATS]
    '''


    num_clones = adata.shape[0]
    agroups, bgroups = preprocess.setup_tcr_groups(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    gex_cluster_nbrs = []
    tcr_cluster_nbrs = []

    for i in range(num_clones):
        ag, bg = agroups[i], bgroups[i]
        gc, tc = clusters_gex[i], clusters_tcr[i]
        gex_cluster_nbrs.append(
            np.nonzero((clusters_gex==gc)&(agroups!=ag)&(bgroups!=bg))[0])
        tcr_cluster_nbrs.append(
            np.nonzero((clusters_tcr==tc)&(agroups!=ag)&(bgroups!=bg))[0])


    dfl = []
    for nbr_frac in all_nbrs:
        gex_nbrs, tcr_nbrs = all_nbrs[nbr_frac]
        stats = _compute_graph_overlap_stats(
            gex_nbrs, tcr_nbrs, num_random_repeats)
        stats['nbr_frac'] = nbr_frac
        stats['graph_overlap_type'] = 'gex_nbr_vs_tcr_nbr'
        dfl.append(stats)

        stats = _compute_graph_overlap_stats(
            gex_nbrs, tcr_cluster_nbrs, num_random_repeats)
        stats['nbr_frac'] = nbr_frac
        stats['graph_overlap_type'] = 'gex_nbr_vs_tcr_cluster'
        dfl.append(stats)

        stats = _compute_graph_overlap_stats(
            tcr_nbrs, gex_cluster_nbrs, num_random_repeats, swaptags=True)
        stats['nbr_frac'] = nbr_frac
        stats['graph_overlap_type'] = 'gex_cluster_vs_tcr_nbr'
        dfl.append(stats)

    results = pd.DataFrame(dfl)


    help_message = """
Here we are assessing overall graph-vs-graph correlation by looking at
the shared edges between TCR and GEX neighbor graphs and comparing
that observed number to the number we would expect if the graphs were
completely uncorrelated. Our null model for uncorrelated graphs is to
take the vertices of one graph and randomly renumber them (permute their
labels). We compare the observed overlap to that expected under this null
model by computing a Z-score, either by permuting one of the graph's
vertices many times to get a mean and standard deviation of the overlap
distribution, or, for large graphs where this is time consuming,
by using a regression model for the
standard deviation. The different rows of this table correspond to the
different graph-graph comparisons that we make in the conga graph-vs-graph
analysis: we compare K-nearest-neighbor graphs for GEX and TCR at different
K values ("nbr_frac" aka neighbor-fraction, which reports K as a fraction
of the total number of clonotypes) to each other and to GEX and TCR "cluster"
graphs in which each clonotype is connected to all the other clonotypes with
the same (GEX or TCR) cluster assignment. For two K values (the default),
this gives 2*3=6 comparisons: GEX KNN graph vs TCR KNN graph, GEX cluster
graph vs TCR KNN graph, and GEX KNN graph vs TCR cluster graph, for each of the
two K values (aka nbr_fracs).

The column to look at is *overlap_zscore*. Higher values indicate more
significant GEX/TCR covariation, with "interesting" levels starting around
zscores of 3-5.

Columns in more detail:

graph_overlap_type: KNN ("nbr") or cluster versus KNN ("nbr") or cluster

nbr_frac: the K value for the KNN graph, as a fraction of total clonotypes

overlap: the observed overlap (number of shared edges) between GEX and TCR
graphs

expected_overlap: the expected overlap under a shuffled null model.

overlap_zscore: a Z-score for the observed overlap computed by subtracting
the expected overlap and dividing by the standard deviation estimated from
shuffling.
"""

    ## store results in adata.uns
    table_tag = GRAPH_VS_GRAPH_STATS
    adata.uns.setdefault('conga_results', {})[table_tag] = results

    # store the help message
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = help_message

    if outfile_prefix is not None:
        util.save_table_and_helpfile(table_tag, adata, outfile_prefix)

    return results


def run_graph_vs_graph(
        adata,
        all_nbrs,
        pval_threshold=1.0, #pvals are multiplied by num_clones, so can be >> 1
        verbose=False,
        outfile_prefix=None,
):
    ''' Runs graph-vs-graph analysis for each nbr_frac in the all_nbrs dictionary

    Returns a (possibly empty) pandas dataframe with the results

    Also sets up the

    'conga_scores'  and 'conga_fdr_values' arrays in adata.obs

    Note that the "pvalues" aka conga_scores have been crudely Bonferroni
    corrected by multiplying by the number of clones. That's why the
    default pval_threshold of 1.0 makes any kind of sense.

    also stores the results dataframe in
        adata.uns['conga_results'][GRAPH_VS_GRAPH]

    and saves them to a tsvfile if outfile_prefix is not None

    '''

    num_clones = adata.shape[0]
    agroups, bgroups = preprocess.setup_tcr_groups(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    nbr_fracs = sorted(all_nbrs.keys())

    # since conga_score= raw_pval*num_clones, max is num_clones if raw_pval=1
    bad_conga_score = num_clones

    conga_scores = np.full((num_clones,), bad_conga_score, dtype=float)

    all_results = []
    all_raw_pvalues = np.full((num_clones, 3*len(nbr_fracs)), 1.0) # for FDR

    for inbr_frac, nbr_frac in enumerate(nbr_fracs):
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]

        print('find_neighbor_neighbor_interactions:')
        results_df, adjusted_pvalues = _find_neighbor_neighbor_interactions(
            adata, nbrs_gex, nbrs_tcr, agroups, bgroups, pval_threshold,
            verbose=verbose)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)
        all_raw_pvalues[:, 3*inbr_frac] = adjusted_pvalues/num_clones

        if not results_df.empty:
            results_df['nbr_frac'] = nbr_frac
            results_df['graph_overlap_type'] = 'gex_nbr_vs_tcr_nbr'
            all_results.append(results_df)

        print('find_neighbor_cluster_interactions:')
        results_df, adjusted_pvalues = _find_neighbor_cluster_interactions(
            adata, nbrs_tcr, clusters_gex, agroups, bgroups, pval_threshold)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)
        all_raw_pvalues[:, 3*inbr_frac+1] = adjusted_pvalues/num_clones
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['graph_overlap_type'] = 'gex_cluster_vs_tcr_nbr'
            results_df.rename(
                columns={'num_neighbors':'num_neighbors_tcr'}, inplace=True)
            all_results.append(results_df)

        print('find_neighbor_cluster_interactions:')
        results_df, adjusted_pvalues = _find_neighbor_cluster_interactions(
            adata, nbrs_gex, clusters_tcr, agroups, bgroups, pval_threshold)
        conga_scores = np.minimum(conga_scores, adjusted_pvalues)
        all_raw_pvalues[:, 3*inbr_frac+2] = adjusted_pvalues/num_clones
        if results_df.shape[0]:
            results_df['nbr_frac'] = nbr_frac
            results_df['graph_overlap_type'] = 'gex_nbr_vs_tcr_cluster'
            results_df.rename(
                columns={'num_neighbors':'num_neighbors_gex'}, inplace=True)
            all_results.append(results_df)

    # fdr calculation
    _, fdr_values, _, _ = multipletests(
        all_raw_pvalues.reshape((num_clones*3*len(nbr_fracs),)),
        alpha=0.05, method='fdr_bh')

    fdr_values = fdr_values.reshape((num_clones, 3*len(nbr_fracs)))
    fdr_values = fdr_values.min(axis=1)

    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)
        ## add some more info to the results table, and sort by conga score
        indices = results_df['clone_index']
        results_df['gex_cluster'] = list(clusters_gex[indices])
        results_df['tcr_cluster'] = list(clusters_tcr[indices])
        for tag in 'va ja cdr3a vb jb cdr3b'.split():
            results_df[tag] = list(adata.obs[tag][indices])
        results_df.sort_values('conga_score', inplace=True)

    else:
        results_df = pd.DataFrame([])

    adata.obs['conga_scores'] = conga_scores
    adata.obs['conga_fdr_values'] = fdr_values

    help_message = """Graph vs graph analysis looks for correlation between GEX and TCR space
by finding statistically significant overlap between two similarity graphs,
one defined by GEX similarity and one by TCR sequence similarity.

Overlap is defined one node (clonotype) at a time by looking for overlap
between that node's neighbors in the GEX graph and its neighbors in the
TCR graph. The null model is that the two neighbor sets are chosen
independently at random.

CoNGA looks at two kinds of graphs: K nearest neighbor (KNN) graphs, where
K = neighborhood size is specified as a fraction of the number of
clonotypes (defaults for K are 0.01 and 0.1), and cluster graphs, where
each clonotype is connected to all the other clonotypes in the same
(GEX or TCR) cluster. Overlaps are computed 3 ways (GEX KNN vs TCR KNN,
GEX KNN vs TCR cluster, and GEX cluster vs TCR KNN), for each of the
K values (called nbr_fracs short for neighbor fractions).

Columns (depend slightly on whether hit is KNN v KNN or KNN v cluster):
conga_score = P value for GEX/TCR overlap * number of clonotypes
mait_fraction = fraction of the overlap made up of 'invariant' T cells
num_neighbors* = size of neighborhood (K)
cluster_size = size of cluster (for KNN v cluster graph overlaps)
clone_index = 0-index of clonotype in adata object

    """

    ## store results in adata.uns
    table_tag = GRAPH_VS_GRAPH
    adata.uns.setdefault('conga_results', {})[table_tag] = results_df

    # store the help message
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = help_message

    if outfile_prefix is not None:
        util.save_table_and_helpfile(table_tag, adata, outfile_prefix)


    return results_df


def run_rank_genes_on_good_biclusters(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        rank_method='wilcoxon',
        rg_tag = 'test', # temporary, removed at the end
        neg_tag='none',
        min_count=5,
        key_added = 'rank_genes_good_biclusters'
):
    ''' Find differentially expressed genes for conga clusters
    each cluster is defined by:
    * good_mask must be True for all clonotypes in cluster
    * all clonotypes in cluster have same GEX (clusters_gex) and
      TCR (clusters_tcr) cluster assignment
    * size of cluster must be at least min_count

    for a while these were called 'biclusters' or 'cluster pairs'
    'cluster pairs' gets shortened to 'clp' in places in the code
    currently we are just calling them 'CoNGA clusters'

    '''
    num_clones = adata.shape[0]

    clp_counts = Counter(
        (x,y) for x,y,z in zip(clusters_gex, clusters_tcr, good_mask) if z)

    #print( clp_counts.most_common())

    vals = [ neg_tag ]*num_clones

    for clp, count in clp_counts.most_common():
        if count< min_count:
            break
        tag = 'clp_{}_{}'.format(clp[0], clp[1])
        inds = np.nonzero(
            [x==clp[0] and y==clp[1] and z
             for x,y,z in zip(clusters_gex, clusters_tcr, good_mask)])[0]
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

    sc.tl.rank_genes_groups(
        adata, groupby=rg_tag, method=rank_method, groups=pos_tags,
        reference='rest', key_added = key_added)

    # remove the temporary obs column
    adata.obs.drop(columns=[rg_tag], inplace=True)


def calc_good_cluster_tcr_features(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        tcr_score_names,
        min_count=5,
        verbose=False, # was True,
):
    ''' Find biased (significantly elevated or lowered) TCR sequence scores
    for conga clusters. Each cluster is defined by:
    * good_mask must be True for all clonotypes in cluster
    * all clonotypes in cluster have same GEX (clusters_gex) and
      TCR (clusters_tcr) cluster assignment
    * size of cluster must be at least min_count

    for a while these were called 'biclusters' or 'cluster pairs'
    'cluster pairs' gets shortened to 'clp' in places in the code
    currently we are just calling them 'CoNGA clusters'

    '''
    num_clones = adata.shape[0]

    # there seems to be a problem with np.nonzero on a pandas series,
    #   which is what these might be if taken from adata.obs:
    clusters_gex = np.array(clusters_gex)
    clusters_tcr = np.array(clusters_tcr)
    good_mask = np.array(good_mask)

    clp_counts = Counter((x,y) for x,y,z in zip(clusters_gex,
                                                clusters_tcr,
                                                good_mask) if z)

    good_clps = [ x for x,y in clp_counts.items() if y>=min_count]

    fake_nbrs_gex = []
    seen = set()
    for cl_gex, cl_tcr, m in zip( clusters_gex, clusters_tcr, good_mask ):
        clp=(cl_gex, cl_tcr)
        if m and clp in good_clps and clp not in seen:
            seen.add(clp)
            fake_nbrs_gex.append(np.nonzero((clusters_gex==cl_gex) &
                                            (clusters_tcr==cl_tcr) &
                                            good_mask)[0])
        else:
            fake_nbrs_gex.append([])

    pval_threshold = 1.
    results_df = gex_nbrhood_rank_tcr_scores(
        adata, fake_nbrs_gex, tcr_score_names, pval_threshold,
        verbose=verbose, prefix_tag = 'good_clp')


    all_tcr_features = {}
    for clp in good_clps:
        # in case results_df doesn't have any rows for this clp
        all_tcr_features[clp] = []

    for row in results_df.itertuples():
        #clp = (row.gex_cluster, row.tcr_cluster)
        clp = (clusters_gex[row.clone_index], clusters_tcr[row.clone_index])
        assert clp in good_clps
        all_tcr_features[clp].append((row.mwu_pvalue_adj, row.ttest_stat,
                                      row.feature))


    for clp in all_tcr_features:
        # plotting code expects name then stat then pvalue
        # but we want them sorted by significance
        all_tcr_features[clp] = [(x[2], x[1], x[0])
                                 for x in sorted(all_tcr_features[clp])]

    return all_tcr_features




def _get_split_mean_var( X, X_sq, mask, mean, mean_sq ):
    ''' Helper function to quickly compute mean and variance
    '''
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
        verbose=False,
        ttest_pval_threshold_for_mwu_calc=None,
        also_return_help_message=False,
):
    ''' Run graph-vs-features analysis comparing the GEX neighbor graph
    to TCR features.

    We also use this to find differential TCR features for conga clusters,
    by passing in a fake GEX neighbor graph.

    pvalues are bonferroni corrected (actually just multiplied by numtests)

    returns a pandas dataframe with the results,

    OR if also_return_help_message,

    returns results, help_message
    '''
    num_clones = adata.shape[0]

    assert len(tcr_score_names) == len(set(tcr_score_names)) # no dups

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    preprocess.add_mait_info_to_adata_obs(adata)
    is_mait = adata.obs['is_invariant']
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold*10

    print('making tcr score table, #features=', len(tcr_score_names))
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
        mean_fg, var_fg, mean_bg, var_bg = _get_split_mean_var(
            score_table, score_table_sq, nbrhood_mask, mean, mean_sq)
        num_fg = np.sum(nbrhood_mask)
        if num_fg < min_num_fg:
            continue
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_fg,
            std1=np.sqrt(np.maximum(var_fg, 1e-12)),
            nobs1=num_fg,
            mean2=mean_bg,
            std2=np.sqrt(np.maximum(var_bg, 1e-12)),
            nobs2=num_clones-num_fg,
            equal_var=False,  # Welch's
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

            _,mwu_pval = mannwhitneyu(score_table[:,ind][nbrhood_mask],
                                      score_table[:,ind][~nbrhood_mask],
                                      alternative='two-sided')
            mwu_pval_adj = mwu_pval * pval_rescale

            # make more stringent
            #if min(pval, mwu_pval_adj) <= pval_threshold:

            if ((mwu_pval_adj <= pval_threshold) or
                (pval <= pval_threshold and
                 mwu_pval_adj <= 10*pval_threshold)):
                if nbrhood_clusters_gex is None: # lazy
                    nbrhood_clusters_gex = clusters_gex[nbrhood_mask]
                    nbrhood_clusters_tcr = clusters_tcr[nbrhood_mask]
                    nbrhood_is_mait = is_mait[nbrhood_mask]

                # get info about the clones most contributing to this skewed
                #  score
                score_name = tcr_score_names[ind]
                score = scores[ind] # ie the t-statistic

                num_top = max(1,num_fg//4)
                if score>0: # score is high
                    top_indices = np.argpartition(
                        score_table[:,ind][nbrhood_mask], -num_top)[-num_top:]
                else: # score is low
                    top_indices = np.argpartition(
                        score_table[:,ind][nbrhood_mask], num_top-1)[:num_top]


                gex_cluster = Counter( nbrhood_clusters_gex[ top_indices ])\
                              .most_common(1)[0][0]
                tcr_cluster = Counter( nbrhood_clusters_tcr[ top_indices ])\
                              .most_common(1)[0][0]
                mait_fraction = np.sum(nbrhood_is_mait[ top_indices ] )\
                                /len(top_indices)

                if verbose and mwu_pval_adj <= pval_threshold:
                    print('gex_{}_score: {:9.2e} {:9.2e} {:7.2f} clp {:2d} {:2d} {:7.3f} {:7.3f} {:15s} {} {} mf: {:.3f} {}'\
                          .format(prefix_tag, pval, mwu_pval_adj, score,
                                  gex_cluster, tcr_cluster, mean_fg[ind],
                                  mean_bg[ind], score_name,
                                  ' '.join(tcrs[ii][0][:3]),
                                  ' '.join(tcrs[ii][1][:3]),
                                  mait_fraction, ii ))

                results.append(dict(ttest_pvalue_adj=pval,
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
    results_df = pd.DataFrame(results)

    help_message = """
    This table has results from a graph-vs-features analysis in which we
    look at the distribution of a set of TCR-defined features over the GEX
    neighbor graph. We look for neighborhoods in the graph that have biased
    score distributions, as assessed by a ttest first, for speed, and then
    by a mannwhitneyu test for nbrhood/score combinations whose ttest P-value
    passes an initial threshold (default is 10* the pvalue threshold).

    Each row of the table represents a single significant association, in other
    words a neighborhood (defined by the central clonotype index) and a
    tcr feature.

    The columns are as follows:

    ttest_pvalue_adj= ttest_pvalue * number of comparisons
    ttest_stat= ttest statistic (sign indicates where feature is up or down)
    mwu_pvalue_adj= mannwhitney-U P-value * number of comparisons
    gex_cluster= the consensus GEX cluster of the clonotypes w/ biased scores
    tcr_cluster= the consensus TCR cluster of the clonotypes w/ biased scores
    num_fg= the number of clonotypes in the neighborhood (including center)
    mean_fg= the mean value of the feature in the neighborhood
    mean_bg= the mean value of the feature outside the neighborhood
    feature= the name of the TCR score
    mait_fraction= the fraction of the skewed clonotypes that have an invariant
        TCR
    clone_index= the index in the anndata dataset of the clonotype that is the
        center of the neighborhood.

    """

    if also_return_help_message: # a little hacky...
        return results_df, help_message
    else:
        return results_df


def tcr_nbrhood_rank_genes_fast(
        adata,
        nbrs_tcr,
        pval_threshold,
        top_n=50,
        verbose=False,
        prefix_tag='nbr',
        min_num_fg=3,
        clone_display_names=None,
        ttest_pval_threshold_for_mwu_calc=None,
        also_return_help_message=False,
):
    ''' Run graph-vs-features analysis comparing the TCR neighbor graph
    to GEX features (ie, expression levels of the different genes)


    Modeled on scanpy rank_genes_groups
    All pvals are crude bonferroni corrected for:
    * number of non-empty nbrhoods in nbrs_tcr and number of genes in adata.raw.X with at least 3 nonzero cells
      (see pval_rescale below)

    returns a pandas dataframe with the results

    -OR- if also_return_help_message is True

    return results_df, help_message

    '''

    ## unpack from adata
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    preprocess.add_mait_info_to_adata_obs(adata)
    is_mait = adata.obs['is_invariant']
    organism = adata.uns['organism']
    ## done unpacking ###############################

    if clone_display_names is None:
        clone_display_names = [
            '{} {}'.format(' '.join(x[0][:3]),' '.join(x[1][:3])) for x in tcrs]

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold * 10

    assert preprocess.check_if_raw_matrix_is_logged(adata)

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
        X2 = np.vstack([np.log1p(np.array(adata.obs['clone_sizes'])),
                        np.log1p(np.argsort(-1*nndists_gex))]).transpose()
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
    bad_gene_mask = np.array([gene_nonzero_counts[x] < min_nonzero_cells
                              for x in range(len(genes))]+
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

        mean_fg, var_fg, mean_bg, var_bg = _get_split_mean_var(
            X, X_sq, nbrhood_mask, mean, mean_sq)
        mean2_fg, var2_fg, mean2_bg, var2_bg = _get_split_mean_var(
            X2, X2_sq, nbrhood_mask, mean2, mean2_sq)

        num_fg = np.sum(nbrhood_mask)
        if num_fg < min_num_fg:
            continue
        mean_fg = np.hstack([mean_fg, mean2_fg])
        mean_bg = np.hstack([mean_bg, mean2_bg]) # note that we dont do the variances...
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_fg,
            std1=np.sqrt(np.maximum(np.hstack([var_fg, var2_fg]),1e-12)),
            nobs1=num_fg,
            mean2=mean_bg,
            std2=np.sqrt(np.maximum(np.hstack([var_bg, var2_bg]), 1e-12)),
            nobs2=num_clones-num_fg,
            equal_var=False  # Welch's
        )

        # scanpy code:
        scores[np.isnan(scores)] = 0.
        pvals [np.isnan(pvals)] = 1.
        pvals [bad_gene_mask] = 1.
        logfoldchanges = np.log2((np.expm1(mean_fg) + 1e-9) /
                                 (np.expm1(mean_bg) + 1e-9))

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
            # here we are looking for genes (or clone_sizes/inverted nndists)
            #  that are LARGER in the forground (fg)
            if is_real_gene:
                col = X_csc[:,ind][nbrhood_mask]
                noncol = X_csc[:,ind][~nbrhood_mask]
                _, mwu_pval = mannwhitneyu( col.todense(), noncol.todense(),
                                            alternative='greater' )
            else:
                col = X2[:,ind-num_real_genes][nbrhood_mask]
                noncol = X2[:,ind-num_real_genes][~nbrhood_mask]
                _, mwu_pval = mannwhitneyu(col, noncol, alternative='greater')
            mwu_pval_adj = mwu_pval * pval_rescale

            # 2021-06-28 make this more stringent: it used to be either/or
            #if min(mwu_pval_adj, pval_adj) < pval_threshold:
            #
            # sometimes MWU seems a little wonky, so allow good ttests also
            # if MWU is not terrible
            if ((mwu_pval_adj < pval_threshold) or
                (pval_adj < pval_threshold and
                 mwu_pval_adj < 10*pval_threshold)):

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

                gex_cluster = Counter(
                    nbrhood_clusters_gex[ top_indices ]).most_common(1)[0][0]
                tcr_cluster = Counter(
                    nbrhood_clusters_tcr[ top_indices ]).most_common(1)[0][0]
                mait_fraction = (np.sum(nbrhood_is_mait[top_indices]) /
                                 len(top_indices))

                if verbose and mwu_pval_adj<=pval_threshold:
                    print(f'tcr_{prefix_tag}_gene: {pval_adj:9.2e} '
                          f'{mwu_pval_adj:9.2e} {log2fold:7.3f} clp '
                          f'{gex_cluster:2d} {tcr_cluster:2d} {gene:8s} '
                          f'{mean_fg[ind]:.4f} {mean_bg[ind]:.4f} {num_fg:4d} '
                          f'{clone_display_names[ii]} mf: {mait_fraction:.3f} '
                          f'{ii} {igene}')

                results.append(dict(ttest_pvalue_adj=pval_adj,
                                    mwu_pvalue_adj=mwu_pval_adj,
                                    log2enr=log2fold,
                                    gex_cluster=gex_cluster,
                                    tcr_cluster=tcr_cluster,
                                    feature=gene,
                                    mean_fg=mean_fg[ind],
                                    mean_bg=mean_bg[ind],
                                    num_fg=num_fg,
                                    clone_index=ii,
                                    mait_fraction=mait_fraction))

        sys.stdout.flush()

    results_df = pd.DataFrame(results)



    help_message = """
    This table has results from a graph-vs-features analysis in which we
    look for genes that are differentially expressed (elevated) in specific
    neighborhoods of the TCR neighbor graph. Differential expression is
    assessed by a ttest first, for speed, and then
    by a mannwhitneyu test for nbrhood/score combinations whose ttest P-value
    passes an initial threshold (default is 10* the pvalue threshold).

    Each row of the table represents a single significant association, in other
    words a neighborhood (defined by the central clonotype index) and a
    gene.

    The columns are as follows:

    ttest_pvalue_adj= ttest_pvalue * number of comparisons
    mwu_pvalue_adj= mannwhitney-U P-value * number of comparisons
    log2enr = log2 fold change of gene in neighborhood (will be positive)
    gex_cluster= the consensus GEX cluster of the clonotypes w/ biased scores
    tcr_cluster= the consensus TCR cluster of the clonotypes w/ biased scores
    num_fg= the number of clonotypes in the neighborhood (including center)
    mean_fg= the mean value of the feature in the neighborhood
    mean_bg= the mean value of the feature outside the neighborhood
    feature= the name of the gene
    mait_fraction= the fraction of the skewed clonotypes that have an invariant
        TCR
    clone_index= the index in the anndata dataset of the clonotype that is the
        center of the neighborhood.

    """

    if also_return_help_message: # a little hacky...
        return results_df, help_message
    else:
        return results_df


    return pd.DataFrame(results)

def setup_fake_nbrs_from_clusters( clusters ):
    ''' Make a fake nbr graph in which one clone in each cluster has a set
    of nbrs to the other cluster members. Everybody else has empty nbr lists.
    For graph-vs-feature correlation analyses.
    '''
    clusters = np.array(clusters) # just in case: pandas Series problem

    fake_nbrs = []
    seen = set()
    for cl in clusters:
        if cl in seen:
            fake_nbrs.append([])
        else:
            seen.add(cl)
            fake_nbrs.append(np.nonzero( clusters==cl )[0])
    return fake_nbrs

## wrapper for graph vs features analysis
def run_graph_vs_features(
        adata,
        all_nbrs,
        pval_threshold= 1., # 'pvalues' are raw_pvalue * num_tests
        outfile_prefix= None,
):
    ''' This runs graph-vs-features analysis comparing the TCR graph to
    GEX features and the GEX graph to TCR features. It also looks for genes
    that are associated with particular TCR V or J segments

    results are stored in adata.uns['conga_results'] under the tags

    GEX_GRAPH_VS_TCR_FEATURES
    TCR_GRAPH_VS_GEX_FEATURES
    TCR_GENES_VS_GEX_FEATURES

    and written to tsvfiles if outfile_prefix is not None

    '''

    util.setup_uns_dicts(adata) # make life easier

    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    nbr_fracs = sorted(all_nbrs.keys())

    #### TCR GRAPHS VS GEX FEATURES ###################
    tcr_graph_results = []

    ## first use the TCRdist KNN nbr graph:

    for nbr_frac in nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results_df, tcr_graph_help_message = tcr_nbrhood_rank_genes_fast(
            adata, nbrs_tcr, pval_threshold, also_return_help_message=True)

        results_df['nbr_frac'] = nbr_frac
        results_df['graph_type'] = 'tcr_nbr'

        tcr_graph_results.append(results_df)



    # now make a TCR cluster graph and use the nbrhoods in there
    # make some fake nbrs-- note that only one clone per cluster has
    #  a nonempty nbrhood
    fake_nbrs_tcr = setup_fake_nbrs_from_clusters(clusters_tcr)
    results_df = tcr_nbrhood_rank_genes_fast(
        adata, fake_nbrs_tcr, pval_threshold, prefix_tag='clust')

    results_df['clone_index'] = -1
    results_df['nbr_frac'] = 0.0
    results_df['graph_type'] = 'tcr_cluster'

    tcr_graph_results.append(results_df)

    results_df = pd.concat(tcr_graph_results, ignore_index=True)
    if not results_df.empty:
        results_df.sort_values('mwu_pvalue_adj', inplace=True)
    table_tag = TCR_GRAPH_VS_GEX_FEATURES

    adata.uns['conga_results'][table_tag] = results_df
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = tcr_graph_help_message

    ##
    ## now make another fake nbr graph defined by TCR gene segment usage
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    tcr_genes_results = []
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
                    fake_nbrs_tcr.append(np.nonzero( genes==g )[0])
                    clone_display_names.append(g)

            results_df = tcr_nbrhood_rank_genes_fast(
                adata, fake_nbrs_tcr, pval_threshold, prefix_tag=seg+ab,
                clone_display_names=clone_display_names)

            if results_df.shape[0]:
                inds = np.array(results_df['clone_index'])
                results_df['gene_segment'] = genes[inds]
                results_df['clone_index'] = -1
                results_df['graph_type'] = 'tcr_genes'
                tcr_genes_results.append(results_df)

    if tcr_genes_results:
        results_df = pd.concat(tcr_genes_results, ignore_index=True)
        results_df.sort_values('mwu_pvalue_adj', inplace=True)
        table_tag = TCR_GENES_VS_GEX_FEATURES
        adata.uns['conga_results'][table_tag] = results_df
        tcr_genes_help_message = tcr_graph_help_message+"""
        In this analysis the TCR graph is defined by
        connecting all clonotypes that have the same VA/JA/VB/JB-gene segment
        (it's run four times, once with each gene segment type)
        """
        adata.uns['conga_results'][table_tag+HELP_SUFFIX]=tcr_genes_help_message


    #### GEX GRAPHS VS TCR FEATURES ###################
    tcr_score_names = tcr_scoring.all_tcr_scorenames[:]
    # also add on the genes that occur in at least 5 clonotypes
    min_gene_count = 5
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    organism = adata.uns['organism']
    organism_genes = all_genes[organism]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs
                      for i_ab in range(2)
                      for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y>=min_gene_count ]
    tcr_score_names += count_reps

    gex_graph_results = []
    for nbr_frac in nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        results_df, gex_graph_help_message = gex_nbrhood_rank_tcr_scores(
            adata, nbrs_gex, tcr_score_names, pval_threshold,
            also_return_help_message = True)
        results_df['nbr_frac'] = nbr_frac
        results_df['graph_type'] = 'gex_nbr'
        gex_graph_results.append(results_df)


    # make some fake nbrs
    fake_nbrs_gex = setup_fake_nbrs_from_clusters(
        clusters_gex)
    results_df = gex_nbrhood_rank_tcr_scores(
        adata, fake_nbrs_gex, tcr_score_names, pval_threshold,
        prefix_tag = 'clust' )
    results_df['clone_index'] = -1
    results_df['nbr_frac'] = 0.0
    results_df['graph_type'] = 'gex_cluster'
    gex_graph_results.append(results_df)

    results_df = pd.concat(gex_graph_results, ignore_index=True)
    if not results_df.empty:
        results_df.sort_values('mwu_pvalue_adj', inplace=True)
    table_tag = GEX_GRAPH_VS_TCR_FEATURES

    adata.uns['conga_results'][table_tag] = results_df
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = gex_graph_help_message


    # write the tables to files:
    if outfile_prefix is not None:
        for table_tag in [TCR_GRAPH_VS_GEX_FEATURES,
                          TCR_GENES_VS_GEX_FEATURES,
                          GEX_GRAPH_VS_TCR_FEATURES]:
            util.save_table_and_helpfile(
                table_tag, adata, outfile_prefix)


    return # all done with graph-vs-features analysis




def find_hotspot_features(
        X,
        nbrs,
        features,
        pval_threshold,
        verbose=False
):
    """ My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of
       Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """
    # right now anyhow; not a strict requirement:
    assert type(X) is sps.csr_matrix

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
        if ii%1000==0:
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
    #print('zeros:', np.sum(mask1), np.sum(mask2), np.sum(mask3))

    H /= np.maximum(1e-9, X_var)

    H[X_var==0] = 0 # set to zero if stddev is 0

    inds = np.argsort(H)[::-1] # decreasing

    # the simple estimate for the variance of H is the total number of neighbors
    # compute the variance
    nbrs_sets = []
    for ii in range(num_clones):
        nbrs_sets.append( frozenset(nbrs[ii]) )

    H_var = 0
    #print('compute H_var')
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
            print('WHOAH var prob? {} {} =?= {}'\
                  .format(feature, X_var[ind], true_std**2))
            continue
        Z = H[ind]
        pvalue_adj = num_features * norm.sf(Z)
        if pvalue_adj > pval_threshold:
            break
        if verbose:
            print('top_var: {} =?= {} {} {} {} {}'\
                  .format(X_var[ind], true_std**2, true_std, Z,
                          pvalue_adj, feature))
        results.append(OrderedDict(Z=Z, pvalue_adj=pvalue_adj, feature=feature))

    return pd.DataFrame(results)




def find_hotspot_genes(
        adata,
        nbrs_tcr,
        pval_threshold,
        verbose=False,
):
    """ Find genes that show a biased distribution across the TCR neighbor graph
    using a simplified version of the Hotspot method from the Yosef lab.

    This is my hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of
       Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    returns a pandas dataframe

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

    X = sps.hstack([X,Y]).tocsr()


    df = find_hotspot_features(
        X, nbrs_tcr,
        list(genes)+['gex_cluster{}'.format(x) for x in range(num_clusters)],
        pval_threshold)

    if df.shape[0]==0:
        return df

    # filter out VDJ genes
    mask = [not util.is_vdj_gene(x, organism, include_constant_regions=True)
            for x in df.feature]

    df = df[mask]
    df['feature_type'] = 'gex'

    if verbose: # show the top 100 hits
        for ii,l in enumerate(df[:100].itertuples()):
            print('hotspot_gene: {:4d} {:9.3f} {:8.1e} {:10s}'\
                  .format(ii, l.Z, l.pvalue_adj, l.feature))

    return df


def find_hotspot_tcr_features(
        adata,
        nbrs_gex,
        pval_threshold,
        min_gene_count=5,
        verbose=False,
):
    """ Find TCR features that show a biased distribution across the
    GEX neighbor graph, using a simplified version of the Hotspot method
    from the Yosef lab.

    This is my hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of
       Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    num_clusters = np.max(adata.obs['clusters_tcr'])+1

    organism_genes = all_genes[organism]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs for i_ab in range(2) for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    features = (tcr_scoring.all_tcr_scorenames +
                count_reps +
                ['tcr_cluster{}'.format(x) for x in range(num_clusters)])

    score_table = tcr_scoring.make_tcr_score_table(adata, features)
    X = sps.csr_matrix(score_table)

    df = find_hotspot_features(X, nbrs_gex, features, pval_threshold)
    df['feature_type'] = 'tcr'

    if verbose: # show the top 100 hits
        for ii,l in enumerate(df[:100].itertuples()):
            print('hotspot_tcr_feature: {:4d} {:9.3f} {:8.1e} {:10s}'\
                  .format(ii, l.Z, l.pvalue_adj, l.feature))

    return df


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
    #gene_df['feature_type'] = 'gex' ## now done in find_hotspot_genes

    tcr_df = find_hotspot_tcr_features(adata , nbrs_gex , pval_threshold  )
    #tcr_df['feature_type'] = 'tcr' ## now done in find_hotspot_tcr_features

    hotspot_df = pd.concat([gene_df, tcr_df], ignore_index=True)

    return hotspot_df

def find_hotspots_wrapper(
        adata,
        all_nbrs,
        pval_threshold=None,
        outfile_prefix=None,
        nbr_fracs = None,
):
    if nbr_fracs is None:
        nbr_fracs = sorted(all_nbrs.keys())
    all_results = []
    for nbr_frac in nbr_fracs:
        results = find_hotspots(adata, all_nbrs[nbr_frac], pval_threshold)
        results['nbr_frac'] = nbr_frac
        all_results.append(results)
    results = pd.concat(all_results, ignore_index=True)
    if not results.empty:
        results.sort_values('pvalue_adj', inplace=True)
    table_tag = HOTSPOT_FEATURES
    adata.uns.setdefault('conga_results',{})[table_tag] = results
    help_message = 'Need a better help message here'
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = help_message

    if outfile_prefix is not None:
        util.save_table_and_helpfile(table_tag, adata, outfile_prefix)

    return results
