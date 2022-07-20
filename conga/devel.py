
################################################################################
##
## this file is a place to put code that is experimental
##  or otherwise under development / not ready for prime time
##
##

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind, poisson
import scipy.sparse as sps
from statsmodels.stats.multitest import multipletests
from collections import Counter, OrderedDict
import scanpy as sc
from . import preprocess
from . import tcr_scoring
from . import util
from . import correlations
from . import plotting
from . import tcr_clumping
from .tcrdist.all_genes import all_genes
import sys
import pandas as pd
from sys import exit
import time #debugging
import random
from os.path import exists
import os
from pathlib import Path

## read a list of human tfs downloaded from
## http://humantfs.ccbr.utoronto.ca/download.php
## which are from Lambert et al (PMID:29425488)
## thanks to Matt Weirauch and collaborators
##
try:
    _tfs_listfile = util.path_to_data / 'Lambert_et_al_PMID_29425488_TF_names_v_1.01.txt'
    human_tf_gene_names = [x.split()[0] for x in open(_tfs_listfile,'r')]
except:
    print('conga.devel:: failed to read human TFs list; prob not a big deal')
    human_tf_gene_names = []

T_class_module_genes = Path.joinpath( Path(util.path_to_data),'T_cell_prediction_markers.tsv')
assert exists(T_class_module_genes)

models_path = Path.joinpath( Path(util.path_to_data), 'prediction_models')
assert exists(models_path)


def compute_distance_correlations( adata, verbose=False ):
    ''' return pvalues, rvalues  (each 1 1d numpy array of shape (num_clones,))
    '''
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    agroups, bgroups = preprocess.setup_tcr_groups(adata)

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




def compute_cluster_interactions(
        aclusters_in,
        bclusters_in,
        barcodes_in,
        barcode2tcr,
        outlog,
        max_pval = 1.0
):

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

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

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
                        _,mwu_pval1 = mannwhitneyu(
                            nbr_scores, non_nbr_scores, alternative='greater',
                            **util.mannwhitneyu_kwargs)

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
                results = correlations.find_hotspot_features(
                    X, nbrs, features, pval_threshold)

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


def calc_nbr_zscore_threshold(nbrs, num_repeats = 50):
    num_clones, num_nbrs = nbrs.shape
    mxl = []
    for i in range(num_repeats):
        scores = np.random.randn(num_clones)
        scores = (scores-np.mean(scores))/np.std(scores) # silly
        scores = (scores + scores[nbrs].sum(axis=1))/(num_nbrs+1)
        mxl.append(np.max(scores))
    threshold = np.median(mxl)
    return threshold


def get_subset_from_scores(
        scores,
        nbrs,
        nbr_zscore_threshold=None,
        require_discrete_features_to_be_present=False,
        verbose=False,
):
    ''' returns a boolean mask

    for a binary (present/absent) feature, the nbr_zscore_threshold is not used

    '''
    num_clones, num_nbrs = nbrs.shape

    if nbr_zscore_threshold is None:
        nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)

    feature_values = set(scores)
    if len(feature_values)==1:
        return np.full((num_clones,), False)

    elif (len(feature_values)==2 and
          abs(min(feature_values)-0)<.01 and
          abs(max(feature_values)-1)<.01):
        # discrete, 0 or 1 feature
        # figure out how many of this gene we would expect in a nbrhood
        # by chance... use that to set a stringent threshold
        gene_count = int(sum(scores))
        threshold = 1
        while threshold < gene_count and threshold <= num_nbrs:
            # what are the odds of seeing a nbrhood > threshold
            pval = hypergeom.sf(threshold, num_clones, num_nbrs+1,
                                gene_count)
            if verbose:
                print('set_threshold_for_gene:', threshold, pval,
                      num_clones, num_nbrs+1, gene_count, feature)
            if pval*num_clones<1:
                break
            threshold += 1
        nbrhood_scores = scores + scores[nbrs].sum(axis=1)
        if require_discrete_features_to_be_present:
            mask = (nbrhood_scores >= threshold-0.5) & (scores>0.5)
        else:
            mask = (nbrhood_scores >= threshold-0.5)
    else:
        zscores = (scores-np.mean(scores))/np.std(scores)
        nbrhood_zscores = (zscores + zscores[nbrs].sum(axis=1))/(num_nbrs+1)

        mask = ((nbrhood_zscores >= nbr_zscore_threshold) &
                (zscores >= nbr_zscore_threshold))
    return mask



def find_distinctive_tcr_features_for_subset(
        adata,
        subset_mask, # boolean mask
        min_gene_count = 5,
):
    num_clones = adata.shape[0]
    subset_mask = np.array(subset_mask)
    assert subset_mask.shape == (adata.shape[0],)

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    organism_genes = all_genes[organism]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs for i_ab in range(2) for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    tcr_scorenames = tcr_scoring.all_tcr_scorenames[:]
    features = tcr_scorenames + count_reps
    score_table = tcr_scoring.make_tcr_score_table(adata, features)

    # from tcr_scoring.py: the feature scores:
    #all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait',
    #  'inkt', 'nndists_tcr'] + list(aa_props_df.columns)

    dfl = []

    ## use wilcoxon for the non-gene features
    for idx, feature in enumerate(features):
        if len(set(score_table[:,idx]))==1:
            continue # only one value
        t_stat, t_pval = ttest_ind(score_table[:,idx][subset_mask], #2-sided
                                   score_table[:,idx][~subset_mask])
        if feature in tcr_scorenames:
            # MWW
            _, mwu_pval = mannwhitneyu(score_table[:,idx][subset_mask],
                                       score_table[:,idx][~subset_mask],
                                       alternative='two-sided',
                                       **util.mannwhitneyu_kwargs)

            # correct for the number of tcr scores
            pvalue_adj = mwu_pval * len(tcr_scorenames)
            t_pval *= len(tcr_scorenames)
            pvalue_type='wilcoxon'
        else:
            # use hypergeometric for gene features
            gene_present = (score_table[:,idx] > 0.5)
            overlap = np.sum(gene_present & subset_mask)
            hgeom_pval = hypergeom.sf(overlap-1, num_clones,
                                      np.sum(subset_mask), np.sum(gene_present))
            pvalue_adj = hgeom_pval * len(count_reps)
            t_pval *= len(count_reps)
            pvalue_type='hypergeometric'

        if pvalue_adj < pvalue_adj_threshold:
            dfl.append(dict(
                feature=feature,
                pvalue_adj=pvalue_adj,
                ttest_stat=t_stat,
                ttest_pvalue_adj=t_pval,
                pvalue_type=pvalue_type,
            ))
    return pd.DataFrame(dfl)




def analyze_interesting_genes(
        adata,
        genes,
        nbrs, # probably TCR nbrs
        min_extreme_clonotypes=5,
        min_gene_count=5,
        pvalue_adj_threshold = 0.05, # threshold on the rank-genes pvalue_adj
        verbose = False,
):
    print('here1')
    nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)
    print('here2')

    var_names = list(adata.raw.var_names)

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    organism_genes = all_genes[adata.uns['organism']]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs for i_ab in range(2) for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    tcr_scorenames = tcr_scoring.all_tcr_scorenames[:]
    tcr_features = tcr_scorenames + count_reps
    score_table = tcr_scoring.make_tcr_score_table(adata, tcr_features)
    print('here3')

    # from tcr_scoring.py: the feature scores:
    #all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait',
    #  'inkt', 'nndists_tcr'] + list(aa_props_df.columns)

    dfl = []

    for gene in genes:
        if gene not in var_names:
            print(f'ERROR: {gene} not present in adata.raw.var_names')
            continue
        gene_index = var_names.index(gene)
        gene_scores = adata.raw.X[:,gene_index].toarray()[:,0]
        mask = get_subset_from_scores(gene_scores, nbrs, nbr_zscore_threshold)

        if verbose:
            print('num_extreme:', gene, np.sum(mask))

        if np.sum(mask) < min_extreme_clonotypes:
            continue

        ## use wilcoxon for the non-gene tcr_features
        for idx, feature in enumerate(tcr_features):
            if len(set(score_table[:,idx]))==1:
                continue # only one value
            t_stat, t_pval = ttest_ind(score_table[:,idx][mask], # 2-sided
                                       score_table[:,idx][~mask])
            if feature in tcr_scorenames:
                # MWW
                _, mwu_pval = mannwhitneyu(score_table[:,idx][mask],
                                           score_table[:,idx][~mask],
                                           alternative='two-sided',
                                           **util.mannwhitneyu_kwargs)

                # correct for the number of tcr scores
                pvalue_adj = mwu_pval * len(tcr_scorenames)
                t_pval *= len(tcr_scorenames)
                pvalue_type='wilcoxon'
            else:
                # use hypergeometric for gene tcr_features
                gene_present = (score_table[:,idx] > 0.5)
                overlap = np.sum(gene_present & mask)
                hgeom_pval = hypergeom.sf(overlap-1, adata.shape[0],
                                          np.sum(mask),
                                          np.sum(gene_present))
                pvalue_adj = hgeom_pval * len(count_reps)
                t_pval *= len(count_reps)
                pvalue_type='hypergeometric'

            if pvalue_adj < pvalue_adj_threshold:
                dfl.append(dict(
                    gene=gene,
                    feature=feature,
                    pvalue_adj=pvalue_adj,
                    pvalue_type=pvalue_type,
                    ttest_stat=t_stat,
                    ttest_pvalue_adj=t_pval,
                    subset_size=np.sum(mask),
                ))
    return pd.DataFrame(dfl)




def analyze_interesting_tcr_features(
        adata,
        features,
        nbrs, # probably GEX nbrs
        min_extreme_clonotypes=5,
        only_high_scores=True,
        pvalue_adj_threshold = 0.05, # threshold on the rank-genes pvalue_adj
        verbose = False,
        require_gene_features_to_be_present = False,
):
    # calculate the feature scores
    scoretable = tcr_scoring.make_tcr_score_table(adata, features)

    # first set a heuristic threshold on the nbr-avged Z-scores
    nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)
    if verbose:
        print('nbr_zscore_threshold:', nbr_zscore_threshold)


    num_clones, num_nbrs = nbrs.shape

    pos_tag, neg_tag = 'pos','neg' # keep same length, we put into np array
    rg_tag = 'test'
    rank_method = 'wilcoxon'
    rg_key = 'rg_feature_'

    dfl = []

    for idx,feature in enumerate(features):
        for sign in [1,-1]:
            if only_high_scores and sign==-1:
                continue

            scores = scoretable[:,idx]
            if sign==-1 and len(set(scores))<=2:
                continue
            mask = get_subset_from_scores(
                scores, nbrs, nbr_zscore_threshold,
                require_discrete_features_to_be_present=
                require_gene_features_to_be_present,
                )

            if verbose:
                print('num_extreme:', feature, sign, np.sum(mask))
            if np.sum(mask) < min_extreme_clonotypes:
                continue

            vals = np.full((num_clones,), neg_tag)
            vals[mask] = pos_tag

            adata.obs[rg_tag] = vals

            # run rank_genes on the cells with extreme values for this feature

            sc.tl.rank_genes_groups(
                adata, groupby=rg_tag, method=rank_method, groups=[pos_tag],
                reference='rest', key_added = rg_key, n_genes = 5000)

            for ii, (name, pvalue_adj, logfoldchange) in enumerate(zip(
                    adata.uns[rg_key]['names'],
                    adata.uns[rg_key]['pvals_adj'],
                    adata.uns[rg_key]['logfoldchanges'])):
                if pvalue_adj[0] > pvalue_adj_threshold:
                    break
                dfl.append(dict(
                    feature=feature,
                    sign=sign,
                    gene=name[0],
                    rank=ii,
                    pvalue_adj=pvalue_adj[0],
                    logfoldchange=logfoldchange[0],
                    subset_size=np.sum(mask),
                    ))

    results = pd.DataFrame(dfl)
    return results



def find_hotspot_nbrhoods(
        adata,
        nbrs_gex,
        nbrs_tcr,
        pval_threshold,
        also_use_cluster_graphs=False
):
    """ This is some experimental code looking to use the HotSpot method
    from the Yosef lab to find nbr-nbr interactions. Still under development.


    My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    nbrs_gex_clusters = setup_fake_nbrs_from_clusters( clusters_gex )
    nbrs_tcr_clusters = setup_fake_nbrs_from_clusters( clusters_tcr )
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
    comparisons = [['gex', nbrs_gex, 'graph', nbrs_tcr],
                   ['gex', nbrs_gex, 'clust', nbrs_tcr_clusters],
                   ['tcr', nbrs_tcr, 'graph', nbrs_gex],
                   ['tcr', nbrs_tcr, 'clust', nbrs_gex_clusters]]

    for feature_nbr_tag, feature_nbrs, graph_tag, graph_nbrs in comparisons:
        if graph_tag == 'clust' and not also_use_cluster_graphs:
            continue

        rows = [] # list of (1,num_clones) csr_matrix'es
        for ii in range(num_clones):
            ii_nbrs = feature_nbrs[ii]
            num_nbrs = len(ii_nbrs)
            data = np.full((num_nbrs,),1.0)
            row_ind = np.full((num_nbrs,),0).astype(int)
            new_row = sps.csr_matrix((data, (row_ind, ii_nbrs)),
                                     shape=(1,num_clones))
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
        df = correlations.find_hotspot_features(X, graph_nbrs, feature_names, pval_threshold)#, verbose=True)

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
        df = pd.concat(dfl, ignore_index=True)
        return df.drop(columns=['feature']) # feature is just str(clone_index), ie superfluous
    else:
        return pd.DataFrame()


def find_gex_cluster_degs(
        adata,
        outfile_prefix,
        min_cluster_size = 3,
):
    if len(set(adata.obs['clusters_gex'])) <= 1:
        print('conga.devel.find_gex_cluster_degs:: too few clusters')
        return

    # look at differentially expressed genes in gex clusters
    obs_tag = 'clusters_gex_for_degs'
    adata.obs[obs_tag] = [str(x) for x in adata.obs['clusters_gex']]
    key_added = 'degs_for_gex_clusters'
    rank_method = 'wilcoxon'
    good_clusters = sorted(
        [x for x,y in Counter(adata.obs[obs_tag]).most_common()
         if y>=min_cluster_size],
        key=int)
    sc.tl.rank_genes_groups(
        adata,
        groupby=obs_tag,
        method=rank_method,
        groups=good_clusters,
        reference='rest',
        key_added=key_added,
    )
    n_genes = 25
    sc.pl.rank_genes_groups(
        adata,
        n_genes=n_genes,
        sharey=False,
        show=False,
        key=key_added,
    )
    pngfile = outfile_prefix+'_gex_cluster_degs.png'
    plt.savefig(pngfile, bbox_inches="tight")
    print('made:', pngfile)


    new_rank_genes_genes, var_group_positions, var_group_labels = [],[],[]
    allow_gene_repeats = False
    min_rank_genes_log2fold_change = 1.0
    max_rank_genes_pval_adj=0.05
    n_genes_for_plotting = 5

    for group in good_clusters:
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
                  .format( log2fold, pval_adj,
                           adata.uns[key_added]['scores'][group][igene],
                           gene, group ) )
            my_genes.append( gene )
        if my_genes:
            var_group_positions.append((len(new_rank_genes_genes),
                                        len(new_rank_genes_genes)+
                                        len(my_genes)-1 ) )
            var_group_labels.append( group )
            new_rank_genes_genes.extend( my_genes )

    if new_rank_genes_genes:
        sc.pl.stacked_violin(
            adata, var_names = new_rank_genes_genes, groupby=obs_tag,
            figsize=(10,n_genes_for_plotting*10),
            use_raw = True, stripplot=True, show=False, swap_axes=True,
            var_group_positions = var_group_positions,
            var_group_labels = var_group_labels,
            var_group_rotation = 1.0 )
        pngfile = outfile_prefix+'_gex_cluster_degs_violin.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:',pngfile)

        sc.pl.dotplot(
            adata, var_names=new_rank_genes_genes, groupby=obs_tag, show=False,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions)
        pngfile = outfile_prefix+'_gex_cluster_degs_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

        # this plot_scatter seems to have moved in scanpy; need to update
        #sc.pl._tools.plot_scatter(adata, 'gex_2d', ncols = 6,
        #                          color = new_rank_genes_genes, show=False,
        #                          use_raw = True, s=40)
        #pngfile = args.outfile_prefix+'_gex_cluster_degs_tsne.png'
        #plt.savefig(pngfile, bbox_inches="tight")
        #print('made:', pngfile)


    if adata.uns['organism'] == 'human_ig':
        # list of B cell marker genes from "Human germinal centres engage
        # memory and naive B cells after influenza vaccination"
        # Turner...Ellebedy, Nature 2020:
        # https://doi.org/10.1038/s41586-020-2711-0
        # note that they say activated B cells are distinguished by lack of CR2
        genes_lines = """GC-Bs BCL6, RGS13, MEF2B, STMN1, ELL3, SERPINA9
        PBs XBP1, IRF4, SEC11C, FKBP11, JCHAIN, PRDM1
        naive TCL1A, IL4R, CCR7, IGHM, IGHD
        act-Bs TBX21, FCRL5, ITGAX, NKG7, ZEB2, CR2
        rest TNFRSF13B, CD27, CD24
        misc IGHA1 IGHA2 IGHG1 IGHG2 IGHG3 IGHG4 IGHE"""\
            .replace(',',' ').split('\n')
        genes, var_group_positions, var_group_labels = [], [], []
        for line in genes_lines:
            my_genes = [x for x in line.split()[1:] if x in adata.raw.var_names]
            print(len(my_genes), line.split())
            if my_genes:
                var_group_positions.append((len(genes),
                                            len(genes)+len(my_genes)-1))
                var_group_labels.append( line.split()[0])
                genes.extend(my_genes)
        sc.pl.dotplot(
            adata, var_names=genes, groupby=obs_tag, show=False,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions)
        pngfile = outfile_prefix+'_gex_cluster_bcell_genes_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

    # show some of our marker genes
    organism = adata.uns['organism']
    genes = (plotting.default_logo_genes[organism] +
             plotting.default_gex_header_genes[organism])
    genes = sorted(set(x for x in genes if x in adata.raw.var_names))
    sc.pl.dotplot(adata, var_names=genes, groupby=obs_tag, show=False)
    pngfile = outfile_prefix+'_gex_cluster_marker_genes_dotplot.png'
    plt.savefig(pngfile, bbox_inches="tight")
    print('made:', pngfile)




def find_hotspot_nbrhoods(
        adata,
        all_nbrs,
        outfile_prefix,
        make_hotspot_nbrhood_logos=False,
        min_cluster_size = 5,
        min_cluster_size_fraction = 0.001,
):

    # My hacky and probably buggy first implementation of the HotSpot method:
    #
    # "Identifying Informative Gene Modules Across Modalities of
    #  Single Cell Genomics"
    # David DeTomaso, Nir Yosef
    # https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    nbr_fracs = sorted(all_nbrs.keys())

    all_results = []
    for nbr_frac in nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        print('find_hotspot_nbrhoods for nbr_frac', nbr_frac)
        results = correlations.find_hotspot_nbrhoods(
            adata, nbrs_gex, nbrs_tcr, pval_threshold=1.0,
            also_use_cluster_graphs=False)

        # the feature_type column is already set in nbrhood_results
        #  to {tcr/gex}_nbrs_vs_graph
        results['nbr_frac'] = nbr_frac
        all_results.append(results)

    # make a plot summarizing the hotspot nbrhood pvals and also
    #  save them to a file
    if all_results:
        nbrhood_results = pd.concat(all_results, ignore_index=True)

        tcrs = preprocess.retrieve_tcrs_from_adata(adata)

        outfile = '{}_hotspot_nbrhoods.tsv'.format(outfile_prefix)
        for iab, ivj in [ (x,y) for x in range(2) for y in range(3) ]:
            key = [ 'va ja cdr3a'.split(), 'vb jb cdr3b'.split()][iab][ivj]
            nbrhood_results[key] = [tcrs[x.clone_index][iab][ivj]
                                    for x in nbrhood_results.itertuples()]
        print('making:', outfile)
        nbrhood_results.to_csv(outfile, sep='\t', index=False)

        num_clones = adata.shape[0]
        nbrhood_pvals = {'gex':np.full((num_clones,), num_clones).astype(float),
                         'tcr':np.full((num_clones,), num_clones).astype(float)}
        for l in nbrhood_results.itertuples():
            assert l.feature_type[3:] == '_nbrs_vs_graph'
            tag = l.feature_type[:3]
            nbrhood_pvals[tag][l.clone_index] = min(
                l.pvalue_adj, nbrhood_pvals[tag][l.clone_index])

        plt.figure(figsize=(12,6))
        for icol, tag in enumerate(['gex','tcr']):
            plt.subplot(1,2,icol+1)
            log10_pvals = np.log10(np.maximum(1e-100, nbrhood_pvals[tag]))
            colors = np.sqrt(np.maximum(0.0, -1*log10_pvals))
            reorder = np.argsort(colors)
            # same umap as feature nbr-type
            xy = adata.obsm['X_{}_2d'.format(tag)]
            vmax = np.sqrt(-1*np.log10(1e-5))
            plt.scatter(xy[reorder,0], xy[reorder,1], c=colors[reorder],
                        vmin=0, vmax=vmax)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.xlabel('{} UMAP1'.format(tag))
            plt.ylabel('{} UMAP2'.format(tag))
            plt.title('{} hotspot nbrhood pvalues'.format(tag))
        plt.tight_layout()
        pngfile = '{}_hotspot_nbrhoods.png'.format(outfile_prefix)
        print('making:', pngfile)
        plt.savefig(pngfile)

        if make_hotspot_nbrhood_logos:
            nbrs_gex, nbrs_tcr = all_nbrs[ max(nbr_fracs) ]
            min_cluster_size = max(
                min_cluster_size,
                int(0.5 + min_cluster_size_fraction * num_clones))
            plotting.make_hotspot_nbrhood_logo_figures(
                adata, nbrs_gex, nbrs_tcr, nbrhood_results,
                min_cluster_size, outfile_prefix,
                pvalue_threshold=1.0)


def make_hotspot_nbrhood_logo_figures(
        adata,
        nbrs_gex,
        nbrs_tcr,
        hotspot_nbrhood_results_df,
        min_cluster_size,
        outfile_prefix,
        pvalue_threshold=1.0,
        **kwargs # passed to make_logo_plots
):
    ''' this is highly experimental and not really used, fyi
    '''
    #unpacking
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    num_clones = adata.shape[0]

    #
    nbrhood_pvals = {'gex':np.full((num_clones,), num_clones).astype(float),
                     'tcr':np.full((num_clones,), num_clones).astype(float)}

    for l in hotspot_nbrhood_results_df.itertuples():
        assert l.feature_type[3:] == '_nbrs_vs_graph'
        tag = l.feature_type[:3]
        nbrhood_pvals[tag][l.clone_index] = min(
            l.pvalue_adj, nbrhood_pvals[tag][l.clone_index])

    fake_clusters = np.full( (num_clones,), 0)

    for tag, pvals in nbrhood_pvals.items():
        # are there any clusters with enough significant pvals??
        if tag == 'gex':
            fake_clusters_gex = clusters_gex
            fake_clusters_tcr = fake_clusters
        else:
            fake_clusters_gex = fake_clusters
            fake_clusters_tcr = clusters_tcr

        pngfile = f'{outfile_prefix}_hotspot_nbrhood_{tag}_clusters.png'

        if 'rank_genes_uns_tag' not in kwargs:
            kwargs['rank_genes_uns_tag']= f'rg_hotspot_{tag}_nbrhood_biclusters'

        plotting.make_cluster_logo_plots_figure(
            adata, pvals, pvalue_threshold, fake_clusters_gex,
            fake_clusters_tcr,
            nbrs_gex, nbrs_tcr, min_cluster_size, pngfile,
            **kwargs)

def make_batch_bias_plots(
        adata,
        nbrhood_results_df, # generated by correlations.py: find_batch_biases
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size_for_logos,
        pvalue_threshold_for_logos,
        outfile_prefix,
        max_color_pvalue = 1e-9, # in the UMAP figure; no brighter beyond there
):
    num_clones = adata.shape[0]
    fake_clusters_gex = np.zeros((num_clones,)).astype(int)
    fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
    batch_bias_pvals = np.full( (num_clones,), num_clones).astype(float)
    for l in nbrhood_results_df.itertuples():
        batch_bias_pvals[ l.clone_index] = min(
            l.pvalue_adj, batch_bias_pvals[l.clone_index])
        fake_clusters_tcr[l.clone_index] = l.cluster_group

    pngfile = f'{outfile_prefix}_batch_bias_logos.png'

    plotting.make_cluster_logo_plots_figure(
        adata, batch_bias_pvals, pvalue_threshold_for_logos, fake_clusters_gex,
        fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile,
        rank_genes_uns_tag='rg_batch_bias_biclusters')

    # make umaps colored by batch_bias pvals
    plt.figure(figsize=(12,6))
    for icol, tag in enumerate(['gex','tcr']):
        plt.subplot(1,2,icol+1)
        colors = np.sqrt( np.maximum(
            0.0, -1*np.log10(np.maximum(1e-100, batch_bias_pvals))))
        reorder = np.argsort(colors)
        xy = adata.obsm['X_{}_2d'.format(tag)] # same umap as feature nbr-type
        vmax = np.sqrt(-1*np.log10(max_color_pvalue))
        #vmax = np.sqrt(-1*np.log10(1e-5))
        plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder],
                     vmin=0, vmax=vmax)
        plt.colorbar()
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel('{} UMAP1'.format(tag.upper()))
        plt.ylabel('{} UMAP2'.format(tag.upper()))
        plt.title('TCR batch_bias pvalues')
    plt.tight_layout()
    pngfile = '{}_batch_bias.png'.format(outfile_prefix)
    print('making:', pngfile)
    plt.savefig(pngfile)

def analyze_proteins(
        adata,
        outfile_prefix,
        exclude_protein_prefixes = ['HTO', 'TCRV', 'TCRv'],
        n_components_prot=20,
        n_components_gex=40,

):
    ''' run pca on the protein data

    '''
    from sklearn.metrics import pairwise_distances
    from scipy.stats import describe

    # tmp hacking:
    num_clones = adata.shape[0]
    # if 'group' not in adata.uns['batch_keys'] or adata.obsm['group'].shape != (num_clones,4):
    #     print('need proper group info:', 'group' in adata.uns['batch_keys'])
    #     return

    ft_varname = util.get_feature_types_varname(adata)
    prot_mask = np.array(adata.raw.var[ft_varname] == 'Antibody Capture')
    for p in exclude_protein_prefixes:
        prot_mask &= ~adata.raw.var_names.str.startswith(p)
    print('used_protein_features:', list(adata.raw.var_names[prot_mask]))
    X = adata.raw.X.tocsc()[:,prot_mask].toarray()
    print(type(X), X.shape)
    mn = np.mean(X, axis=0)
    #print('mn:', mn)
    sd = np.std(X, axis=0)
    #print('sd:', sd)
    X = (X-mn[np.newaxis,:])/sd[np.newaxis,:]
    pca = PCA()
    X_pca = pca.fit_transform(X)
    print('pca_variance:', pca.explained_variance_ratio_[:n_components_prot+1])
    X_pca_prot = X_pca[:,:n_components_prot]
    X_pca_gex  = adata.obsm['X_pca_gex'][:,:n_components_gex]
    nrandom = 1000
    inds = np.random.permutation(adata.shape[0])[:nrandom]
    D_prot = pairwise_distances(X_pca_prot[inds,:])
    D_gex = pairwise_distances(X_pca_gex[inds,:])
    print('gex_dists: ', describe(D_gex.ravel()))
    print('prot_dists:', describe(D_prot.ravel()))
    print('dist_ratio:', np.mean(D_prot.ravel()) / np.mean(D_gex.ravel()))

    X_pca_combo = np.hstack([X_pca_gex, X_pca_prot])
    print('X_pca_combo:', X_pca_combo.shape)

    for tag, X_pca in [['prot', X_pca_prot], ['combo', X_pca_combo]]:

        adata.obsm['X_pca'] = X_pca
        adata.obsm[f'X_pca_{tag}'] = np.copy(X_pca) # force copy
        print('neighbors:', tag, X_pca.shape)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=X_pca.shape[1])

        print('umap:', tag)
        sc.tl.umap(adata)
        adata.obsm['X_umap_'+tag] = adata.obsm['X_umap']
        cluster_key_added = 'louvain_'+tag
        resolution = 1.0
        print('louvain:', tag)
        sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
        adata.obs['clusters_'+tag] = np.copy(adata.obs[cluster_key_added]).astype(int)
        adata.obsm['X_{}_2d'.format(tag)] = adata.obsm['X_umap_'+tag]




    xy_tags = ['gex','prot','combo']
    special_genes = ['ZBTB16','RTKN2','CD79A']
    trbv = 'TRBV12-5'
    genes = ['CD45RA_p']+special_genes+['min_special',trbv]

    if 'group' in adata.obs_keys():
        other = ['CD4','CD8','group_0','group_1','group_2','group_3']

        group_counts = adata.obsm['group'].astype(float)
        group_fracs = group_counts/(np.sum(group_counts,axis=1)[:,np.newaxis])
        assert np.max(group_fracs)<1.001
        assert group_fracs.shape == (num_clones, 4)
    else:
        other = ['CD4','CD8']

    nrows, ncols, plotno = len(xy_tags), 1+len(genes)+len(other), 0
    plt.figure(figsize=(ncols*2.5, nrows*3))

    var_names = list(adata.raw.var_names)
    for xy_tag in xy_tags:
        xy = adata.obsm[f'X_{xy_tag}_2d']

        # calculate nbrs for averaging
        obsm_tag = f'X_pca_{xy_tag}'
        assert obsm_tag in adata.obsm_keys()
        num_nbrs = max(1, min(num_clones//10, max(20, num_clones//500)))
        nbr_frac = (num_nbrs+0.1)/num_clones
        all_nbrs = preprocess.calc_nbrs(adata, [nbr_frac], obsm_tag_gex=obsm_tag, obsm_tag_tcr=None)
        nbrs = all_nbrs[nbr_frac][0]
        assert num_nbrs == nbrs.shape[1]
        all_nbr_avg_colors = {}
        for colortag in [xy_tag] + other + genes:
            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            cmap = 'viridis'
            title = colortag
            nbr_avg = False
            sort_order = True
            if colortag == 'min_special':
                colors = np.hstack( [all_nbr_avg_colors[x][:,np.newaxis] for x in special_genes])
                assert colors.shape == (num_clones, len(special_genes))
                colors = np.min(colors, axis=1)
            elif colortag == trbv:
                colors = np.array(adata.obs['vb'].str.startswith(trbv)).astype(float)
                nbr_avg = True
            elif colortag == 'CD4':
                colors = adata.raw.X[:,var_names.index('CD4_p')].toarray()[:,0]
            elif colortag == 'CD8':
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8_p')].toarray()[:,0]+
                                adata.raw.X[:,var_names.index('CD8a_p')].toarray()[:,0])
            elif colortag in genes:
                colors = adata.raw.X[:,var_names.index(colortag)].toarray()[:,0]
                nbr_avg = True
            elif colortag in xy_tags:
                colors = adata.obs[f'clusters_{colortag}']
                cmap = 'tab20'
                title = f'{colortag} clusters ({np.max(colors)+1})'
                sort_order = False
            else:
                assert colortag.startswith('group_')
                g = int(colortag.split('_')[1])
                colors = group_fracs[:,g]
                nbr_avg = True

            if nbr_avg:
                colors = (colors + colors[nbrs].sum(axis=1))/(num_nbrs+1)
                title = f'{title} nbr-avg-{num_nbrs}'
                all_nbr_avg_colors[colortag] = colors

            if sort_order:
                reorder = np.argsort(colors)
            else:
                reorder = np.arange(num_clones)

            plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
            plt.xticks([],[])
            plt.yticks([],[])
            #plt.xlabel(xy_tag)
            plt.ylabel(xy_tag)
            plt.title(title)
            plt.colorbar()
    pngfile = outfile_prefix+'_prot_combo_special_panels.png'
    print('making:', pngfile)
    plt.tight_layout()
    plt.savefig(pngfile)

def analyze_special_genes(
        adata,
        outfile_prefix,
        special_genes = ['ZBTB16', 'RTKN2', 'CD79A'],
        obs_tag_clusters = 'clusters_gex',
        obsm_tag_nbrs = 'X_pca_gex',
        obsm_tag_xy = 'X_gex_2d',
        trbvs = ['TRBV12-5', 'TRBV4-1'],
        other_genes = ['KLRB1'],
        nbrs = None,
        num_nbrs = None,

):
    ''' analyze special genes

    '''
    #from sklearn.metrics import pairwise_distances
    from scipy.stats import linregress, mannwhitneyu

    assert obsm_tag_xy in adata.obsm_keys()
    assert obsm_tag_nbrs in adata.obsm_keys()

    num_clones = adata.shape[0]

    genes = ['clusters', 'CD4', 'CD8'] + other_genes + special_genes + ['min_special']+trbvs

    ncols = 4
    nrows = (len(genes)+len(trbvs)-1)//ncols + 1
    plotno=0

    plt.figure(figsize=(ncols*3, nrows*3))

    var_names = list(adata.raw.var_names)
    xy = adata.obsm[obsm_tag_xy]

    # calculate nbrs for averaging
    if num_nbrs is None:
        if nbrs is None:
            num_nbrs = max(1, min(num_clones//10, max(20, num_clones//500)))
        else:
            num_nbrs = nbrs.shape[1]

    if num_nbrs and nbrs is None:
        print('computing nbrs:', adata.shape, num_nbrs)
        nbr_frac = (num_nbrs+0.1)/num_clones
        all_nbrs = preprocess.calc_nbrs(adata, [nbr_frac], obsm_tag_gex=obsm_tag_nbrs, obsm_tag_tcr=None)
        nbrs = all_nbrs[nbr_frac][0]
        assert num_nbrs == nbrs.shape[1]

    all_nbr_avg_colors = {}
    for colortag in genes:
        plotno += 1
        plt.subplot(nrows, ncols, plotno)
        cmap = 'viridis'
        title = colortag
        nbr_avg = False
        sort_order = True
        if colortag == 'min_special':
            colors = np.hstack( [all_nbr_avg_colors[x][:,np.newaxis] for x in special_genes])
            assert colors.shape == (num_clones, len(special_genes))
            colors = np.min(colors, axis=1)
            all_nbr_avg_colors[colortag] = colors # already averaged...
        elif colortag == 'clusters':
            colors = adata.obs[obs_tag_clusters]
            cmap = 'tab20'
            title = f'clusters ({np.max(colors)+1})'
            sort_order = False
        elif colortag in trbvs:
            colors = np.array(adata.obs['vb'].str.startswith(colortag)).astype(float)
            nbr_avg = True
        elif colortag == 'CD4':
            if 'CD4_p' in var_names:
                title = 'CD4_prot'
                colors = adata.raw.X[:,var_names.index('CD4_p')].toarray()[:,0]
            elif 'CD4' in var_names:
                title = 'CD4_gex'
                colors = adata.raw.X[:,var_names.index('CD4')].toarray()[:,0]
            else:
                colors = np.zeros((num_clones,))
        elif colortag == 'CD8':
            if 'CD8_p' in var_names:
                title = 'CD8_prot'
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8_p')].toarray()[:,0]+
                                 adata.raw.X[:,var_names.index('CD8a_p')].toarray()[:,0])
            else:
                title = 'CD8_gex'
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8A')].toarray()[:,0]+
                                 adata.raw.X[:,var_names.index('CD8B')].toarray()[:,0])
        else:
            colors = adata.raw.X[:,var_names.index(colortag)].toarray()[:,0]
            nbr_avg = True

        if nbr_avg:
            if num_nbrs:
                colors = (colors + colors[nbrs].sum(axis=1))/(num_nbrs+1)
            title = f'{title} nbr-avg-{num_nbrs}'
            all_nbr_avg_colors[colortag] = colors

        if sort_order:
            reorder = np.argsort(colors)
        else:
            reorder = np.arange(num_clones)

        plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.xlabel(xy_tag)
        #plt.ylabel(xy_tag)
        plt.title(title)
        plt.colorbar()

    for trbv in trbvs:
        xvals = all_nbr_avg_colors[trbv]
        yvals = all_nbr_avg_colors['min_special']
        slope, intercept, r, p, err = linregress(xvals, yvals)
        mask = np.array(adata.obs['vb'].str.startswith(trbv))
        _,mwu_p = mannwhitneyu(yvals[mask], yvals[~mask],
                               **util.mannwhitneyu_kwargs)
        print(f'trbv_vs_special: {trbv} {r:.2f} {p:.2e} {mwu_p:.2e}')
        plotno+= 1
        plt.subplot(nrows, ncols, plotno)
        plt.scatter(xvals, yvals, s=5)
        plt.title(f'{trbv} r: {r:.2f}\nP: {p:.1e} {mwu_p:.1e}')
        plt.xlabel(trbv)
        plt.ylabel('min_special')

    pngfile = outfile_prefix+'_special_genes.png'
    print('making:', pngfile)
    plt.tight_layout()
    plt.savefig(pngfile, dpi=200)
    plt.close()

    return nbrs # may be useful for caller?


def analyze_CD4_CD8(
        adata,
        nbrs_gex,
        outfile_prefix,
        cd4_genes_default = ['CD4'],
        cd8_genes_default = ['CD8A', 'CD8B'],
        cd4_proteins_default = ['CD4_p'],
        cd8_proteins_default = ['CD8_p', 'CD8a_p'],
):
    import matplotlib.pyplot as plt
    if adata.uns['organism'] != 'human':
        # only really makes sense for human tcr AB
        print('analyze_CD4_CD8 incompatible organism:', adata.uns['organism'])
        return

    clusters_gex = np.array(adata.obs['clusters_gex'])
    num_clones = adata.shape[0]

    for ptag in ['gex','prot']:
        if ptag == 'gex':
            cd4_genes = [ x for x in cd4_genes_default if x in adata.raw.var_names]
            cd8_genes = [ x for x in cd8_genes_default if x in adata.raw.var_names]
        else:
            cd4_genes = [ x for x in cd4_proteins_default if x in adata.raw.var_names]
            cd8_genes = [ x for x in cd8_proteins_default if x in adata.raw.var_names]

        if not cd4_genes or not cd8_genes:
            continue

        all_gex = {}

        for gene in cd4_genes + cd8_genes:
            index = list(adata.raw.var_names).index(gene)
            all_gex[gene] = adata.raw.X[:,index].toarray()[:,0]

        # look at neighborhoods
        xvals = np.copy(all_gex[cd8_genes[0]])
        for g in cd8_genes[1:]:
            xvals += all_gex[g]
        xvals/= len(cd8_genes)

        yvals = np.copy(all_gex[cd4_genes[0]])
        for g in cd4_genes[1:]:
            yvals += all_gex[g]
        yvals/= len(cd4_genes)

        num_nbrs = nbrs_gex.shape[1]

        xvals_nbr_avged = (xvals + xvals[nbrs_gex].sum(axis=1))/(num_nbrs+1)
        yvals_nbr_avged = (yvals + yvals[nbrs_gex].sum(axis=1))/(num_nbrs+1)

        npanels = 2 + 3 + len(cd4_genes) + len(cd8_genes)
        ncols = 4
        nrows = (npanels-1)/ncols + 1

        plt.figure(figsize=(ncols*3, nrows*3))
        plt.subplot(nrows, ncols, 1)
        plt.scatter(xvals_nbr_avged, yvals_nbr_avged, alpha=0.25, c=clusters_gex, cmap='tab20')
        if ptag=='prot':
            plt.xlim((0, 5))
            plt.ylim((0, 5))
        else:
            plt.xlim((0, plt.xlim()[1]))
            plt.ylim((0, plt.ylim()[1]))
        plt.xlabel('+'.join(cd8_genes))
        plt.ylabel('+'.join(cd4_genes))


        plt.subplot(nrows, ncols, 2)
        num_clusters = np.max(clusters_gex)+1
        vals = []
        for c in range(num_clusters):
            mask = clusters_gex==c
            x = np.sum(xvals[mask])/np.sum(mask)
            y = np.sum(yvals[mask])/np.sum(mask)
            vals.append((x,y,c))
        plt.scatter([x[0] for x in vals], [x[1] for x in vals], c=[x[2] for x in vals], cmap='tab20', s=50)
        for x,y,c in vals:
            plt.text(x,y,str(c))
        plt.xlabel('+'.join(cd8_genes))
        plt.ylabel('+'.join(cd4_genes))
        if ptag=='prot':
            plt.xlim((0, 5))
            plt.ylim((0, 5))
            plt.plot([0,5],[0,5],':k')
        else:
            plt.xlim((0, plt.xlim()[1]))
            plt.ylim((0, plt.ylim()[1]))

        plotno=2
        xy = adata.obsm['X_gex_2d']
        for color in ['cluster','n_genes','n_counts']+cd4_genes+cd8_genes:
            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            if color == 'cluster':
                colors = clusters_gex
                cmap='tab20'
                reorder = np.arange(len(colors))
            elif color in ['n_genes','n_counts']:
                colors = np.array(adata.obs[color])
                cmap = 'viridis'
                reorder = np.argsort(colors)
            else:
                colors = all_gex[color]
                cmap = 'viridis'
                reorder = np.argsort(colors)
            plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
            plt.colorbar()
            plt.title(color)

        pngfile = f'{outfile_prefix}_cd4_cd8_{ptag}.png'
        plt.savefig(pngfile)
        print('making:', pngfile)



def filter_cells_by_ribo_norm(adata):
    ''' returns  new filtered adata
    will normalize_and_log_the_raw_matrix if not already done
    '''
    from scipy.stats import gaussian_kde

    preprocess.normalize_and_log_the_raw_matrix(adata)

    X_norm = adata.raw.X

    X_norm_ribo = []
    var_names = [x.upper() for x in adata.raw.var_names ] ## NOTE UPPER-- help for mouse
    for ig,g in enumerate(var_names):
        if g.startswith('RP') and len(g)>=4 and g[2] in 'SL' and g[3].isdigit():
            X_norm_ribo.append( X_norm[:,ig].toarray() )

    X_norm_ribo = np.sum(np.hstack( X_norm_ribo ), axis=1 )/len(X_norm_ribo)

    dens = gaussian_kde(X_norm_ribo)
    # look for two maxima and a minima in between
    # or look for the closest minimum to 2.0
    step = 0.01
    x=2.0
    while True:
        y0 = dens(x)[0]
        yl = dens(x-step)[0]
        yr = dens(x+step)[0]
        print('filter_cells_by_ribo_norm:: minfind:',x,y0,yl,yr)
        if yl < yr:
            if yl < y0:
                x = x-step
            else:
                # found a local min
                ribo_norm_threshold=x
                break
        else:
            if yr < y0:
                x = x+step
            else:
                ribo_norm_threshold=x
                break
    print('filter_cells_by_ribo_norm:: ribo_norm localmin:',ribo_norm_threshold)

    mask = (X_norm_ribo>ribo_norm_threshold)
    print('filter_cells_by_ribo_norm::', adata.shape, 'downto:', np.sum(mask))
    sys.stdout.flush()
    adata = adata[mask,:].copy()
    assert np.sum(mask) == adata.shape[0]
    return adata


def _get_pvalue_from_rvalue(r, n):
    ''' hacky helper stolen from scipy.stats.linregress
    alternative='two-sided'
    '''
    from scipy import stats
    TINY = 1.0e-20
    df = n - 2  # Number of degrees of freedom
    # n-2 degrees of freedom because 2 has been used up
    # to estimate the mean and standard deviation
    t = r * np.sqrt(df / ((1.0 - r + TINY)*(1.0 + r + TINY)))
    #t, prob = _ttest_finish(df, t, alternative)
    prob = 2*stats.t.cdf(-abs(t), df)
    return prob



def run_genes_versus_features(
        adata,
        tcr_features = None,
        verbose=True,
):
    if tcr_features is None:
        tcr_features = (tcr_scoring.all_tcr_scorenames+
                        [x+'_frac' for x in tcr_scoring.amino_acids])


    print('making tcr scoretable:', adata.shape[0], len(tcr_features))
    organism = adata.uns['organism']
    num_clones = adata.shape[0]
    num_features = len(tcr_features)

    # rare genes seem to give spurious correlations
    min_cells_per_gene = max(5, min(50, 0.001*num_clones))

    tcr_scoretable = tcr_scoring.make_tcr_score_table(
        adata, tcr_features, verbose=True)

    # run with variable genes and with all genes
    # pay a factor of 2 in pvalues for that
    variable_genes = list(adata.var_names)
    all_genes = list(adata.raw.var_names)

    Y = tcr_scoretable.T.copy()
    assert Y.shape == (len(tcr_features), num_clones)

    Y_mean = Y.mean(axis=1)
    Y_std = Y.std(axis=1)
    Y = (Y - Y_mean[:,None])/Y_std[:,None]
    low_var_mask = np.abs(Y_std)<1e-6
    print('run_features_versus_features:: dropping tcr features with',
          'low variance:', [x for x,y in zip(tcr_features, low_var_mask)
                            if y])
    Y[low_var_mask,:] = 0. # dont want nans
    print('Y nans:', np.sum(np.isnan(Y)))

    dfl = []
    for is_variable in [True,False]:
        if is_variable:
            genes = variable_genes

            X = adata.X.T.copy()
            X = (X - X.mean(axis=1)[:,None])/X.std(axis=1)[:,None]

            C = X@Y.T/X.shape[1]

        else:
            genes = all_genes

            X = adata.raw.X.T
            print('dotting... sparse X and Y', X.shape, Y.shape)
            D = X.dot(Y.T)

            print('mean/std of sparse X...', X.shape)
            X_sq = X.multiply(X)

            X_mean = X.mean(axis=1).A1
            X_sq_mean = X_sq.mean(axis=1).A1
            X_mean_sq = X_mean**2
            X_var = X_sq_mean - X_mean_sq
            X_std = X_var**0.5

            X_low_var_mask = np.abs(X_std)<1e-6
            print('run_genes_versus_features:: low variance genes:',
                  [x for x,y in zip(genes, X_low_var_mask) if y])
            X_std = np.maximum(X_std, 1e-6)

            print('correlations...')
            Y_sum = Y.sum(axis=1)
            C = (D - np.outer(X_mean, Y_sum))/(num_clones * X_std[:,None])

        print('C nans:', np.sum(np.isnan(C)))
        assert np.sum(np.isnan(C))==0

        num_genes = len(genes)
        assert C.shape == (num_genes, num_features)
        pval_rescale = 2*num_genes*(num_features-np.sum(low_var_mask))

        inds = np.argsort(C.ravel())[::-1]

        for r in range(2):
            feature_counts = Counter()
            rinds = inds if r==0 else inds[::-1]
            for ind in rinds:
                top = np.unravel_index(ind, C.shape)
                rvalue = C[top]
                gene = genes[top[0]]
                if util.is_vdj_gene(gene, organism):
                    continue
                if not is_variable and gene in variable_genes:
                    continue
                feature = tcr_features[top[1]]
                tfstar = 'TF' if (organism == 'human' and
                                  gene in human_tf_gene_names) else '  '

                if is_variable:
                    gene_vals = X[top[0],:]
                else:
                    gene_vals = X[top[0],:].toarray()[0]
                feature_vals = tcr_scoretable.T[top[1],:]
                gene_nonzero = np.sum(np.abs(gene_vals)>1e-3)
                feature_nonzero = np.sum(np.abs(feature_vals)>1e-3)
                if gene_nonzero < min_cells_per_gene:
                    continue
                #pval = pval_rescale*res.pvalue
                pval = pval_rescale*_get_pvalue_from_rvalue(rvalue, num_clones)

                # res = linregress(gene_vals, feature_vals)
                # print('equal?', res.pvalue*pval_rescale, pval, res.rvalue,
                #       rvalue, gene_vals.shape, feature_vals.shape,
                #       type(gene_vals), type(feature_vals),
                #       num_clones, gene_nonzero,
                #       feature_nonzero, gene, feature, r)

                if pval>1:
                    break
                feature_counts[feature] += 1
                dfl.append(dict(
                    gene=gene,
                    feature=feature,
                    is_variable_gene=is_variable,
                    pvalue_adj=pval,
                    rvalue=rvalue,
                    is_human_tf = (gene in human_tf_gene_names),
                    num_nonzero_gene=gene_nonzero,
                    num_nonzero_feature=feature_nonzero,
                ))
                vstar = 'V' if gene in variable_genes else ' '
                if verbose and feature_counts[feature] < 3:
                    print(f'{tfstar} {vstar} {C[top]:6.3f} {rvalue:6.3f}',
                          f'{pval:9.2e} {gene_nonzero:5d} {feature_nonzero:5d}',
                          f'{gene:10s} {feature:15s}')


    results = pd.DataFrame(dfl)
    results.sort_values('pvalue_adj', inplace=True)
    return results

def filter_sorted_gene_features(
        adata,
        features,
        max_features,
):
    ''' return new_features, duplicates

    where duplicates is a dict mapping from rep features to their cluster peers
    '''

    A = make_raw_feature_scores_table(
        features, ['gex']*len(features), adata, verbose=True).T

    assert A.shape == (len(features), adata.shape[0])

    print('calc correlations:', A.shape)
    start = time.time()
    C = 1-distance.squareform(
        distance.pdist(A, metric='correlation'), force='tomatrix')
    print('took:', time.time() - start)

    new_features, duplicates = filter_sorted_features_by_correlation(
        features, C, max_features)

    return new_features, duplicates

def filter_sorted_gex_features_in_dataframe(
        results,
        adata,
        max_features,
):
    ''' returns new dataframe, and duplicates dict

    results should have the columns 'feature' and 'feature_type'

    '''

    gex_features = [x for x,y in zip(results.feature, results.feature_type)
                    if y=='gex']

    new_features, duplicates = filter_sorted_gene_features(
        adata, gex_features, max_features)

    new_features = set(new_features)

    mask = [ (y=='tcr') or (x in new_features)
             for x,y in zip(results.feature, results.feature_type)]

    return results[mask].copy(), duplicates


def find_tcr_clumping_single_chain(
        tcrs,
        organism,
        tmpfile_prefix = 'tmp',
        radii = [6, 12, 24, 48],
        num_random_samples_multiplier = 100, # ie, 100 * num_clones
        pvalue_threshold = 1.0,
        verbose=True,
        preserve_vj_pairings = False,
        pseudocount=0.25,
):
    ''' Returns a pandas dataframe with the following columns:
    - chain (A or B)
    - clone_index
    - nbr_radius
    - pvalue_adj
    - num_nbrs
    - expected_num_nbrs
    - raw_count
    - va, ja, cdr3a, vb, jb, cdr3b (ie, the 6 tcr cols for clone_index clone)
    - clumping_group: clonotypes within each other's significant nbr_radii are linked
    - clump_type: string, either 'global' or 'intra_gex_cluster' (latter only if also_find_clumps_within_gex_clusters=T)

    '''

    outprefix = f'{tmpfile_prefix}_{random.random()}_tcr_clumping'

    num_clones = len(tcrs)
    num_random_samples = num_random_samples_multiplier * num_clones

    agroups, bgroups = preprocess.setup_tcr_groups_for_tcrs(tcrs)

    tmpfiles = [] # for cleanup

    ## compute background neighbor counts at the specified radii
    tcr_clumping.estimate_background_tcrdist_distributions(
        organism, tcrs, max(radii), num_random_samples,
        tmpfile_prefix=outprefix,
        preserve_vj_pairings=preserve_vj_pairings,
        save_unpaired_dists=True,
        #nocleanup=True,
    )

    afile = f'{outprefix}_dists.txt_A.txt'
    bfile = f'{outprefix}_dists.txt_B.txt'
    if not (exists(afile) and exists(bfile)):
        print('ERROR find_tcr_clumping_single_chain::',
              'estimate_background_tcrdist_distributions failed')
        return pd.DataFrame()
    tmpfiles.extend([afile, bfile])

    acounts = np.cumsum(np.loadtxt(afile, dtype=int), axis=1)[:,radii]
    bcounts = np.cumsum(np.loadtxt(bfile, dtype=int), axis=1)[:,radii]

    # newcol = np.full((len(tcrs),1), num_random_samples)
    # acounts = np.hstack([acounts[:,radii], newcol])
    # bcounts = np.hstack([bcounts[:,radii], newcol])


    big_dfl = [] # will hold results for both chains as dataframes

    tcrs_file = outprefix +'_tcrs.tsv'
    pd.DataFrame({
        'va'   : [x[0][0] for x in tcrs],
        'cdr3a': [x[0][2] for x in tcrs],
        'vb'   : [x[1][0] for x in tcrs],
        'cdr3b': [x[1][2] for x in tcrs],
    }).to_csv(tcrs_file, sep='\t', index=False)
    tmpfiles.append(tcrs_file)

    for chain, chain_groups, bg_counts in [['A', agroups, acounts],
                                           ['B', bgroups, bcounts]]:

        # find neighbors in fg tcrs up to max(radii) ############################

        if os.name == 'posix':
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain')
        else:
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain.exe')


        db_filename = Path.joinpath(
            Path(util.path_to_tcrdist_cpp_db), f'tcrdist_info_{organism}.txt')

        tcrdist_threshold = max(radii)

        cmd = (f'{exe} -f {tcrs_file} -t {tcrdist_threshold} -d {db_filename}'
               f' -c {chain} -o {outprefix}')

        util.run_command(cmd, verbose=True)

        nbr_indices_filename = f'{outprefix}_nbr{tcrdist_threshold}_indices.txt'
        nbr_distances_filename = f'{outprefix}_nbr{tcrdist_threshold}_distances.txt'
        tmpfiles.extend([nbr_indices_filename, nbr_distances_filename])

        if not exists(nbr_indices_filename) or not exists(nbr_distances_filename):
            print('find_neighbors failed:', exists(nbr_indices_filename),
                  exists(nbr_distances_filename))
            exit(1)

        all_nbrs = []
        all_distances = []
        for group, line1, line2 in zip(chain_groups,
                                       open(nbr_indices_filename,'r'),
                                       open(nbr_distances_filename,'r')):
            nbrs = [int(x) for x in line1.split()]
            dists = [int(x) for x in line2.split()]
            assert len(nbrs) == len(dists)
            mask = [chain_groups[x] != group for x in nbrs]
            all_nbrs.append([x for x,m in zip(nbrs,mask) if m])
            all_distances.append([x for x,m in zip(dists,mask) if m])
        assert len(all_nbrs) == num_clones

        # use poisson to find nbrhoods with more tcrs than expected;
        #  have to handle agroups/bgroups
        dfl = []

        is_clumped = np.full((num_clones,), False)

        all_raw_pvalues = np.full((num_clones, len(radii)), 1.0)

        for ii in range(num_clones):
            ii_bg_counts = bg_counts[ii]
            ii_dists = all_distances[ii]
            for irad, radius in enumerate(radii):
                num_nbrs = sum(x<=radius for x in ii_dists)
                if num_nbrs<1:
                    continue
                max_nbrs = np.sum(chain_groups != chain_groups[ii])
                pval = hypergeom.sf(
                    num_nbrs-1, # total fg nbr tcrs (-1)
                    max_nbrs+num_random_samples, # total fg+bg tcrs
                    num_nbrs+ii_bg_counts[irad], # total nbr tcrs
                    max_nbrs) # total fg tcrs
                mu = max_nbrs * ii_bg_counts[irad]/num_random_samples
                all_raw_pvalues[ii, irad] = pval
                pval *= len(radii) * num_clones # simple multiple test correction
                #hgeom_pval *= len(radii) * num_clones # simple multiple test correction
                if pval <= pvalue_threshold:
                    is_clumped[ii] = True
                    # count might just be pseudocount
                    raw_count = ii_bg_counts[irad]
                    if verbose:
                        atcr_str = ' '.join(tcrs[ii][0][:3])
                        btcr_str = ' '.join(tcrs[ii][1][:3])
                        print(f'tcr_nbrs_global: {num_nbrs:2d} {mu:9.6f}',
                              f'radius: {radius:2d} pval: {pval:9.1e}',
                              f'{raw_count:9.1f} tcr: {atcr_str} {btcr_str}')

                    dfl.append( OrderedDict(
                        chain=chain,
                        clone_index=ii,
                        nbr_radius=radius,
                        pvalue_adj=pval,
                        #hgeom_pvalue_adj=hgeom_pval,
                        num_nbrs=num_nbrs,
                        expected_num_nbrs=mu,
                        raw_count=raw_count,
                        va   =tcrs[ii][0][0],
                        ja   =tcrs[ii][0][1],
                        cdr3a=tcrs[ii][0][2],
                        vb   =tcrs[ii][1][0],
                        jb   =tcrs[ii][1][1],
                        cdr3b=tcrs[ii][1][2],
                    ))

        results_df = pd.DataFrame(dfl)
        if results_df.shape[0] == 0:
            big_dfl.append(results_df)
            continue # to next chain


        # compute FDR values in addition to the simple adjusted pvalues
        _, fdr_values, _, _ = multipletests(
            all_raw_pvalues.reshape(-1), alpha=0.05, method='fdr_bh')
        fdr_values = fdr_values.reshape((num_clones, len(radii))).min(axis=1)
        # right now we don't assign fdr values for intra-gex cluster clumping
        results_df['clonotype_fdr_value'] = [
            fdr_values[x] for x in results_df.clone_index]

        # identify groups of related hits?
        all_clumped_nbrs = {}
        for l in results_df.itertuples():
            ii = l.clone_index
            radius = l.nbr_radius
            clumped_nbrs = set(x for x,y in zip(all_nbrs[ii], all_distances[ii])
                               if y<= radius and is_clumped[x])
            clumped_nbrs.add(ii)
            if ii in all_clumped_nbrs:
                all_clumped_nbrs[ii] = all_clumped_nbrs[ii] | clumped_nbrs
            else:
                all_clumped_nbrs[ii] = clumped_nbrs


        clumped_inds = sorted(all_clumped_nbrs.keys())
        assert len(clumped_inds) == np.sum(is_clumped)

        # make nbrs symmetric
        for ii in clumped_inds:
            for nbr in all_clumped_nbrs[ii]:
                assert nbr in all_clumped_nbrs
                all_clumped_nbrs[nbr].add(ii)

        all_smallest_nbr = {}
        for ii in clumped_inds:
            all_smallest_nbr[ii] = min(all_clumped_nbrs[ii])

        while True:
            updated = False
            for ii in clumped_inds:
                nbr = all_smallest_nbr[ii]
                new_nbr = min(nbr, np.min([all_smallest_nbr[x]
                                           for x in all_clumped_nbrs[ii]]))
                if nbr != new_nbr:
                    all_smallest_nbr[ii] = new_nbr
                    updated = True
            if not updated:
                break

        # define clusters, choose cluster centers
        clusters = np.array([0]*num_clones) # 0 if not clumped

        cluster_number=0
        cluster_sizes = Counter()
        for ii in clumped_inds:
            nbr = all_smallest_nbr[ii]
            if ii==nbr:
                cluster_number += 1
                members = [ x for x,y in all_smallest_nbr.items() if y==nbr]
                clusters[members] = cluster_number
                cluster_sizes[cluster_number] = len(members)

        for ii, nbrs in all_clumped_nbrs.items():
            for nbr in nbrs:
                # confirm single-linkage clusters
                assert clusters[ii] == clusters[nbr]

        assert not np.any(clusters[is_clumped]==0)
        assert np.all(clusters[~is_clumped]==0)

        # reorder the clumping groups by size
        remap = {x[0]:i+1 for i,x in enumerate(cluster_sizes.most_common())}
        remap[0] = 0

        results_df['clumping_group'] = [ remap[clusters[x.clone_index]]
                                         for x in results_df.itertuples()]

        results_df.sort_values('pvalue_adj', inplace=True)

        big_dfl.append(results_df)

    # cleanup the tmpfiles
    for tmpfile in tmpfiles:
        if exists(tmpfile):
            os.remove(tmpfile)

    results = pd.concat(big_dfl)
    return results


def make_single_chain_clumping_logos(
        results_df, # generated by find_tcr_clumping_single_chain
        adata,
        nbrs_gex,
        nbrs_tcr,
        outfile_prefix,
        min_cluster_size_for_logos=3,
        pvalue_threshold_for_logos=1.0, #pvalues are crude bonferroni corrected
        max_color_pvalue = 1e-16, # in the UMAP figure; no brighter beyond there
        **logo_plot_args,
):
    ''' Make tcr clumping logos for each chain individually, based
    on previously calculated single-chain clumping results dataframe
    '''
    num_clones = adata.shape[0]

    for chain in 'AB':
        fake_clusters_gex = np.zeros((num_clones,)).astype(int)
        fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
        clumping_pvals = np.full( (num_clones,), num_clones).astype(float)


        for l in results_df.itertuples():
            if l.chain == chain:
                clumping_pvals[ l.clone_index] = min(l.pvalue_adj,
                                                     clumping_pvals[l.clone_index])
                fake_clusters_tcr[l.clone_index] = l.clumping_group

        pngfile = f'{outfile_prefix}_{chain}_chain_clumping_logos.png'

        if 'rank_genes_uns_tag' not in logo_plot_args:
            logo_plot_args['rank_genes_uns_tag'] = f'rg_tcr_clumping_{chain}_biclusters'

        if 'conga_scores_name' not in logo_plot_args:
            logo_plot_args['conga_scores_name'] = f'TCR{chain} clumping'

        if 'show_real_clusters_gex' not in logo_plot_args:
            logo_plot_args['show_real_clusters_gex'] = True

        plotting.make_cluster_logo_plots_figure(
            adata, clumping_pvals, pvalue_threshold_for_logos,
            fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
            min_cluster_size_for_logos, pngfile,
            **logo_plot_args)





def assign_cd4_and_cd8_by_clusters(
        adata,
        key_added = 'cd4_or_cd8',
        clustering_resolution = 2.0,
        n_gex_pcs = 40,
        n_neighbors = 10,
        verbose = False,
):
    ''' adds new string column with name=key_added to adata.obs, values are
    'cd4' or 'cd8'

    Does not respect clone definitions (expanded clones may span both cd4 and cd8)

    assumes:
    * data have been preprocessed/scaled
    * data have not been reduced to a single cell per clone
    * pca may or may not have been called

    '''


    assert adata.uns['organism'] == 'human' # tmp hack

    # run pca if necessary
    if 'X_pca_gex' not in adata.obsm_keys():
        n_gex_pcs = min(adata.shape[0]-1, n_gex_pcs)
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_gex_pcs)
        adata.obsm['X_pca_gex'] = adata.obsm['X_pca']

    # run clustering with higher resolution
    adata.obsm['X_pca'] = adata.obsm['X_pca_gex']
    n_pcs = adata.obsm['X_pca'].shape[1]

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    cluster_key_added = 'leiden_gex_for_cd4_vs_cd8'

    sc.tl.leiden(adata, resolution=clustering_resolution, key_added=cluster_key_added)


    clusters_gex = np.array(adata.obs[cluster_key_added].astype(int))

    num_clones = adata.shape[0]
    all_gex = {}

    cd4_genes = ['CD4']
    cd8_genes = ['CD8A', 'CD8B']

    for gene in cd4_genes+cd8_genes:
        if gene not in adata.raw.var_names:
            print('WARNING assign_cd4_and_cd8_by_clusters: missing', gene)
            all_gex[gene] = np.zeros((num_clones,))
        else:
            index = adata.raw.var_names.get_loc(gene)
            all_gex[gene] = adata.raw.X[:,index].toarray()[:,0]

    cd8_gex = np.copy(all_gex[cd8_genes[0]])
    for g in cd8_genes[1:]:
        cd8_gex += all_gex[g]
    cd8_gex /= len(cd8_genes)

    cd4_gex = np.copy(all_gex[cd4_genes[0]])
    for g in cd4_genes[1:]:
        cd4_gex += all_gex[g]
    cd4_gex /= len(cd4_genes)

    num_clusters = np.max(clusters_gex)+1
    cd4_rescale = 1/0.6
    cd8_rescale = 1/2.0
    good_mask = np.full((num_clones,), False)
    dfl = []
    cd48s = np.full((num_clones,), 'UNK')
    for c in range(num_clusters):
        mask = clusters_gex==c
        cd4 = cd4_rescale * np.sum(cd4_gex[mask])/np.sum(mask)
        cd8 = cd8_rescale * np.sum(cd8_gex[mask])/np.sum(mask)
        tag = 'cd4' if cd4>cd8 else 'cd8'
        cd48s[mask] = tag
        if verbose:
            print(f'split_by_CD4_CD8: cluster {c:2d} {np.sum(mask):4d} {tag}',
                  f'cd4: {cd4:9.3f} cd8: {cd8:9.3f}')
        dfl.append(dict(
            cluster=c,
            size=np.sum(mask),
            cd4=cd4,
            cd8=cd8,
            ))
    adata.obs[key_added] = cd48s
    return pd.DataFrame(dfl)



def split_into_cd4_and_cd8_subsets(
        adata,
        allow_split_clones = True, # allow clones that span subsets
        max_iterations = 5,
        verbose = False,
        min_cells_for_iteration = 25,
):
    ''' returns adata_cd4, adata_cd8

    assumes:
    * data have been preprocessed/scaled
    * data have not been reduced to a single cell per clone
    * pca may or may not have been called

    takes the consensus for expanded clones, but ties contribute cells to both

    '''

    assign_cd4_and_cd8_by_clusters(adata, verbose=verbose)

    ad4 = adata[adata.obs.cd4_or_cd8 == 'cd4'].copy()
    ad8 = adata[adata.obs.cd4_or_cd8 == 'cd8'].copy()

    for r in range(max_iterations):
        print('start sizes:', r, ad4.shape[0], ad8.shape[0])
        if (ad4.shape[0] < min_cells_for_iteration or
            ad8.shape[0] < min_cells_for_iteration):
            break

        # look at clone overlap
        # tcrs_4 = preprocess.retrieve_tcrs_from_adata(ad4)
        # tcrs_8 = preprocess.retrieve_tcrs_from_adata(ad8)
        # counts_4 = Counter(tcrs_4)
        # counts_8 = Counter(tcrs_8)

        # for tcr,c4 in counts_4.items():
        #     c8 = counts_8[tcr]
        #     if c8:
        #         print(f'split_clone: {r} {c4:3d} {c8:3d} {tcr[0][2]} {tcr[1][2]}')

        # make new adatas
        dfl = []
        for ad in [ad4,ad8]:
            ad.uns['organism'] = 'human'
            assign_cd4_and_cd8_by_clusters(ad, verbose=verbose)

        ad4_new = ad4[ad4.obs.cd4_or_cd8 == 'cd4'].concatenate(
            ad8[ad8.obs.cd4_or_cd8 == 'cd4'], index_unique=None)
        ad8_new = ad8[ad8.obs.cd4_or_cd8 == 'cd8'].concatenate(
            ad4[ad4.obs.cd4_or_cd8 == 'cd8'], index_unique=None)

        # convergence?
        old_cd8_barcodes = set(ad8.obs.index)
        new_cd8_barcodes = set(ad8_new.obs.index)
        print(f'old_cd8_barcodes: {len(old_cd8_barcodes)} new_cd8_barcodes: '
              f'{len(new_cd8_barcodes)}')
        if (len(old_cd8_barcodes) == len(new_cd8_barcodes) ==
            len(old_cd8_barcodes&new_cd8_barcodes)):
            print('converged')
            break

        ad4 = ad4_new
        ad8 = ad8_new

    # now map from barcodes to cd4/cd8
    cd4_barcodes = set(ad4.obs.index)
    cd8_barcodes = set(ad8.obs.index)
    assert len(cd4_barcodes & cd8_barcodes) == 0

    # barcode_map is dict from barcode to integer (4 or 8)
    barcode_map = {x:4 for x in cd4_barcodes}
    barcode_map.update({x:8 for x in cd8_barcodes})


    # now assign clone by clone
    tcrs = preprocess.retrieve_tcrs_from_adata(adata) # has duplicates
    tcrs_sorted = sorted(set(tcrs)) # no duplicates
    tcr_to_clone_id = {x:i for i,x in enumerate(tcrs_sorted)}
    clone_ids = np.array([tcr_to_clone_id[x] for x in tcrs])

    all_barcodes = np.array(adata.obs.index)

    new_barcode_map = {}
    for clone_id, tcr in enumerate(tcrs_sorted):
        barcodes = list(all_barcodes[clone_ids==clone_id])
        counts = Counter([barcode_map[x] for x in barcodes]).most_common()
        if len(counts)>1 and verbose:
            print('split clone:', counts)

        if len(counts)>1 and counts[0][1] == counts[1][1] and allow_split_clones:
            # clone is split between cd4 and cd8, split by cell assignments
            for barcode in barcodes:
                new_barcode_map[barcode] = barcode_map[barcode]
        else:
            for barcode in barcodes:
                new_barcode_map[barcode] = counts[0][0]

    cd4_mask = np.array([new_barcode_map[x] == 4 for x in all_barcodes])
    cd8_mask = np.array([new_barcode_map[x] == 8 for x in all_barcodes])

    assert all(cd4_mask | cd8_mask)
    assert not any(cd4_mask & cd8_mask)

    adata_cd4 = adata[cd4_mask].copy()
    adata_cd8 = adata[cd8_mask].copy()

    return adata_cd4, adata_cd8


def predict_T_cell_subset (adata, model = 'gradientBoost', use_raw = True):

    assert adata.uns['organism'] == 'human', 'Model currently only supports human'
    assert model in ['logreg','MLP','gradientBoost']
    model_file = Path.joinpath( Path(models_path), f'{model}.sav')
    assert exists(model_file), 'Model not found'
    T_class_model = pd.read_pickle(model_file)
    gs_df = pd.read_csv(T_class_module_genes, sep='\t' )

    
    # check if scores logged already
    if 'Tcell_subset_scores' not in adata.uns_keys():
        subsets = gs_df.cluster.unique()
        for subset in subsets:
            gene_set = gs_df.gene[gs_df.cluster == subset].tolist()
            sc.tl.score_genes(adata, gene_set, score_name = subset, use_raw=use_raw)
        cols_keeps = subsets.tolist()
        print("logging adata.uns['Tcell_subset_scores']")
        adata.uns['Tcell_subset_scores'] = adata.obs[cols_keeps].copy()
        adata.obs = adata.obs.drop(columns = cols_keeps )

    # prediction based on module scores
    X = adata.uns['Tcell_subset_scores']
    y_pred = T_class_model.predict(X)
    adata.obs[f'T_cell_subset_{model}'] = y_pred
    print(f"Prediction logged in adata.obs['T_cell_subset_{model}']")

    return adata






## TEMPORARY CODE GRAVEYARD

# def _include_redundant_feature(count):
#     ''' count = 1 for first feature exceeding threshold, then 2, etc

#     idea:

#     '''
#     assert count>0
#     count += 2 # so 3 is the first (skip), then 4 (keep), then 5-7 (skip), etc
#     return abs(2**int(np.log2(count)) - count)<.1

    # if max_redundant_features is not None:
    #     feature_nbrs = {}

    #     # filter features by correlation
    #     # will subset: features, feature_types, feature_labels, A, nrows
    #     if verbose: print('computing feature correlations:', A.shape)
    #     C = 1-distance.squareform(distance.pdist(A, metric='correlation'),
    #                               force='tomatrix')
    #     if verbose: print('DONE computing feature correlations:', A.shape)
    #     feature_nbr_counts = [0]*len(features)
    #     feature_mask = np.full(len(features), True)
    #     for ii,f1 in enumerate(features):
    #         if verbose and ii%10==0: print('redundancy checking...', ii)
    #         # am I too close to a previous feature?
    #         for jj in range(ii-1):
    #             if ( feature_mask[jj] and
    #                  (feature_types is None or
    #                   feature_types[ii] == feature_types[jj])):
    #                 if C[ii,jj] > redundancy_threshold:
    #                     feature_nbr_counts[jj] += 1
    #                     if feature_nbr_counts[jj] > max_redundant_features:
    #                         count= feature_nbr_counts[jj]-max_redundant_features
    #                         if not _include_redundant_feature(count):
    #                             print('skip:', ii, jj, count)
    #                             feature_mask[ii] = False
    #                             feature_nbrs.setdefault(
    #                                 features[jj],[]).append(f1)
    #                         else:
    #                             print('keep:', ii, jj, count)
    #                         break
    # print('filling the score array for', len(features), 'features')
    # A = np.zeros((len(features), adata.shape[0]))
    # for ii,feature in enumerate(features):
    #     if verbose and ii%10==0:
    #         print('filling:', ii, len(features), feature)
    #     scores = get_raw_feature_scores(
    #         feature, adata,
    #         feature_types if feature_types is None else feature_types[ii])

    #     mn, std = np.mean(scores), np.std(scores)
    #     scores = (scores-mn)
    #     if std!=0:
    #         scores /= std

    #     if compute_nbr_averages:
    #         num_neighbors = nbrs.shape[1]
    #         scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

    #     if feature_types is not None:
    #         if feature_types[ii] == dist_tag:
    #             # lighten these up a bit since they will be nbr correlated
    #             scores *= rescale_factor_for_self_features

    #     A[ii,:] = scores

################################################################################
##
## this file is a place to put code that is experimental
##  or otherwise under development / not ready for prime time
##
##

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind, poisson
import scipy.sparse as sps
from scipy.sparse import issparse, csr_matrix
from statsmodels.stats.multitest import multipletests
from collections import Counter, OrderedDict
import scanpy as sc
from . import preprocess
from . import tcr_scoring
from . import util
from . import correlations
from . import plotting
from . import tcr_clumping
from .tcrdist.all_genes import all_genes
import sys
import pandas as pd
from sys import exit
import time #debugging
import random
from os.path import exists
import os
from pathlib import Path

## read a list of human tfs downloaded from
## http://humantfs.ccbr.utoronto.ca/download.php
## which are from Lambert et al (PMID:29425488)
## thanks to Matt Weirauch and collaborators
##
try:
    _tfs_listfile = util.path_to_data / 'Lambert_et_al_PMID_29425488_TF_names_v_1.01.txt'
    human_tf_gene_names = [x.split()[0] for x in open(_tfs_listfile,'r')]
except:
    print('conga.devel:: failed to read human TFs list; prob not a big deal')
    human_tf_gene_names = []

def compute_distance_correlations( adata, verbose=False ):
    ''' return pvalues, rvalues  (each 1 1d numpy array of shape (num_clones,))
    '''
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    agroups, bgroups = preprocess.setup_tcr_groups(adata)

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




def compute_cluster_interactions(
        aclusters_in,
        bclusters_in,
        barcodes_in,
        barcode2tcr,
        outlog,
        max_pval = 1.0
):

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

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

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
                        _,mwu_pval1 = mannwhitneyu(
                            nbr_scores, non_nbr_scores, alternative='greater',
                            **util.mannwhitneyu_kwargs)

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
                results = correlations.find_hotspot_features(
                    X, nbrs, features, pval_threshold)

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


def calc_nbr_zscore_threshold(nbrs, num_repeats = 50):
    num_clones, num_nbrs = nbrs.shape
    mxl = []
    for i in range(num_repeats):
        scores = np.random.randn(num_clones)
        scores = (scores-np.mean(scores))/np.std(scores) # silly
        scores = (scores + scores[nbrs].sum(axis=1))/(num_nbrs+1)
        mxl.append(np.max(scores))
    threshold = np.median(mxl)
    return threshold


def get_subset_from_scores(
        scores,
        nbrs,
        nbr_zscore_threshold=None,
        require_discrete_features_to_be_present=False,
        verbose=False,
):
    ''' returns a boolean mask

    for a binary (present/absent) feature, the nbr_zscore_threshold is not used

    '''
    num_clones, num_nbrs = nbrs.shape

    if nbr_zscore_threshold is None:
        nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)

    feature_values = set(scores)
    if len(feature_values)==1:
        return np.full((num_clones,), False)

    elif (len(feature_values)==2 and
          abs(min(feature_values)-0)<.01 and
          abs(max(feature_values)-1)<.01):
        # discrete, 0 or 1 feature
        # figure out how many of this gene we would expect in a nbrhood
        # by chance... use that to set a stringent threshold
        gene_count = int(sum(scores))
        threshold = 1
        while threshold < gene_count and threshold <= num_nbrs:
            # what are the odds of seeing a nbrhood > threshold
            pval = hypergeom.sf(threshold, num_clones, num_nbrs+1,
                                gene_count)
            if verbose:
                print('set_threshold_for_gene:', threshold, pval,
                      num_clones, num_nbrs+1, gene_count, feature)
            if pval*num_clones<1:
                break
            threshold += 1
        nbrhood_scores = scores + scores[nbrs].sum(axis=1)
        if require_discrete_features_to_be_present:
            mask = (nbrhood_scores >= threshold-0.5) & (scores>0.5)
        else:
            mask = (nbrhood_scores >= threshold-0.5)
    else:
        zscores = (scores-np.mean(scores))/np.std(scores)
        nbrhood_zscores = (zscores + zscores[nbrs].sum(axis=1))/(num_nbrs+1)

        mask = ((nbrhood_zscores >= nbr_zscore_threshold) &
                (zscores >= nbr_zscore_threshold))
    return mask



def find_distinctive_tcr_features_for_subset(
        adata,
        subset_mask, # boolean mask
        min_gene_count = 5,
):
    num_clones = adata.shape[0]
    subset_mask = np.array(subset_mask)
    assert subset_mask.shape == (adata.shape[0],)

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    organism_genes = all_genes[organism]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs for i_ab in range(2) for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    tcr_scorenames = tcr_scoring.all_tcr_scorenames[:]
    features = tcr_scorenames + count_reps
    score_table = tcr_scoring.make_tcr_score_table(adata, features)

    # from tcr_scoring.py: the feature scores:
    #all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait',
    #  'inkt', 'nndists_tcr'] + list(aa_props_df.columns)

    dfl = []

    ## use wilcoxon for the non-gene features
    for idx, feature in enumerate(features):
        if len(set(score_table[:,idx]))==1:
            continue # only one value
        t_stat, t_pval = ttest_ind(score_table[:,idx][subset_mask], #2-sided
                                   score_table[:,idx][~subset_mask])
        if feature in tcr_scorenames:
            # MWW
            _, mwu_pval = mannwhitneyu(score_table[:,idx][subset_mask],
                                       score_table[:,idx][~subset_mask],
                                       alternative='two-sided',
                                       **util.mannwhitneyu_kwargs)

            # correct for the number of tcr scores
            pvalue_adj = mwu_pval * len(tcr_scorenames)
            t_pval *= len(tcr_scorenames)
            pvalue_type='wilcoxon'
        else:
            # use hypergeometric for gene features
            gene_present = (score_table[:,idx] > 0.5)
            overlap = np.sum(gene_present & subset_mask)
            hgeom_pval = hypergeom.sf(overlap-1, num_clones,
                                      np.sum(subset_mask), np.sum(gene_present))
            pvalue_adj = hgeom_pval * len(count_reps)
            t_pval *= len(count_reps)
            pvalue_type='hypergeometric'

        if pvalue_adj < pvalue_adj_threshold:
            dfl.append(dict(
                feature=feature,
                pvalue_adj=pvalue_adj,
                ttest_stat=t_stat,
                ttest_pvalue_adj=t_pval,
                pvalue_type=pvalue_type,
            ))
    return pd.DataFrame(dfl)




def analyze_interesting_genes(
        adata,
        genes,
        nbrs, # probably TCR nbrs
        min_extreme_clonotypes=5,
        min_gene_count=5,
        pvalue_adj_threshold = 0.05, # threshold on the rank-genes pvalue_adj
        verbose = False,
):
    print('here1')
    nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)
    print('here2')

    var_names = list(adata.raw.var_names)

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    organism_genes = all_genes[adata.uns['organism']]
    counts = Counter([organism_genes[x[i_ab][j_vj]].count_rep
                      for x in tcrs for i_ab in range(2) for j_vj in range(2)])
    count_reps = [x for x,y in counts.most_common() if y >= min_gene_count ]

    tcr_scorenames = tcr_scoring.all_tcr_scorenames[:]
    tcr_features = tcr_scorenames + count_reps
    score_table = tcr_scoring.make_tcr_score_table(adata, tcr_features)
    print('here3')

    # from tcr_scoring.py: the feature scores:
    #all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait',
    #  'inkt', 'nndists_tcr'] + list(aa_props_df.columns)

    dfl = []

    for gene in genes:
        if gene not in var_names:
            print(f'ERROR: {gene} not present in adata.raw.var_names')
            continue
        gene_index = var_names.index(gene)
        gene_scores = adata.raw.X[:,gene_index].toarray()[:,0]
        mask = get_subset_from_scores(gene_scores, nbrs, nbr_zscore_threshold)

        if verbose:
            print('num_extreme:', gene, np.sum(mask))

        if np.sum(mask) < min_extreme_clonotypes:
            continue

        ## use wilcoxon for the non-gene tcr_features
        for idx, feature in enumerate(tcr_features):
            if len(set(score_table[:,idx]))==1:
                continue # only one value
            t_stat, t_pval = ttest_ind(score_table[:,idx][mask], # 2-sided
                                       score_table[:,idx][~mask])
            if feature in tcr_scorenames:
                # MWW
                _, mwu_pval = mannwhitneyu(score_table[:,idx][mask],
                                           score_table[:,idx][~mask],
                                           alternative='two-sided',
                                           **util.mannwhitneyu_kwargs)

                # correct for the number of tcr scores
                pvalue_adj = mwu_pval * len(tcr_scorenames)
                t_pval *= len(tcr_scorenames)
                pvalue_type='wilcoxon'
            else:
                # use hypergeometric for gene tcr_features
                gene_present = (score_table[:,idx] > 0.5)
                overlap = np.sum(gene_present & mask)
                hgeom_pval = hypergeom.sf(overlap-1, adata.shape[0],
                                          np.sum(mask),
                                          np.sum(gene_present))
                pvalue_adj = hgeom_pval * len(count_reps)
                t_pval *= len(count_reps)
                pvalue_type='hypergeometric'

            if pvalue_adj < pvalue_adj_threshold:
                dfl.append(dict(
                    gene=gene,
                    feature=feature,
                    pvalue_adj=pvalue_adj,
                    pvalue_type=pvalue_type,
                    ttest_stat=t_stat,
                    ttest_pvalue_adj=t_pval,
                    subset_size=np.sum(mask),
                ))
    return pd.DataFrame(dfl)




def analyze_interesting_tcr_features(
        adata,
        features,
        nbrs, # probably GEX nbrs
        min_extreme_clonotypes=5,
        only_high_scores=True,
        pvalue_adj_threshold = 0.05, # threshold on the rank-genes pvalue_adj
        verbose = False,
        require_gene_features_to_be_present = False,
):
    # calculate the feature scores
    scoretable = tcr_scoring.make_tcr_score_table(adata, features)

    # first set a heuristic threshold on the nbr-avged Z-scores
    nbr_zscore_threshold = calc_nbr_zscore_threshold(nbrs)
    if verbose:
        print('nbr_zscore_threshold:', nbr_zscore_threshold)


    num_clones, num_nbrs = nbrs.shape

    pos_tag, neg_tag = 'pos','neg' # keep same length, we put into np array
    rg_tag = 'test'
    rank_method = 'wilcoxon'
    rg_key = 'rg_feature_'

    dfl = []

    for idx,feature in enumerate(features):
        for sign in [1,-1]:
            if only_high_scores and sign==-1:
                continue

            scores = scoretable[:,idx]
            if sign==-1 and len(set(scores))<=2:
                continue
            mask = get_subset_from_scores(
                scores, nbrs, nbr_zscore_threshold,
                require_discrete_features_to_be_present=
                require_gene_features_to_be_present,
                )

            if verbose:
                print('num_extreme:', feature, sign, np.sum(mask))
            if np.sum(mask) < min_extreme_clonotypes:
                continue

            vals = np.full((num_clones,), neg_tag)
            vals[mask] = pos_tag

            adata.obs[rg_tag] = vals

            # run rank_genes on the cells with extreme values for this feature

            sc.tl.rank_genes_groups(
                adata, groupby=rg_tag, method=rank_method, groups=[pos_tag],
                reference='rest', key_added = rg_key, n_genes = 5000)

            for ii, (name, pvalue_adj, logfoldchange) in enumerate(zip(
                    adata.uns[rg_key]['names'],
                    adata.uns[rg_key]['pvals_adj'],
                    adata.uns[rg_key]['logfoldchanges'])):
                if pvalue_adj[0] > pvalue_adj_threshold:
                    break
                dfl.append(dict(
                    feature=feature,
                    sign=sign,
                    gene=name[0],
                    rank=ii,
                    pvalue_adj=pvalue_adj[0],
                    logfoldchange=logfoldchange[0],
                    subset_size=np.sum(mask),
                    ))

    results = pd.DataFrame(dfl)
    return results



def find_hotspot_nbrhoods(
        adata,
        nbrs_gex,
        nbrs_tcr,
        pval_threshold,
        also_use_cluster_graphs=False
):
    """ This is some experimental code looking to use the HotSpot method
    from the Yosef lab to find nbr-nbr interactions. Still under development.


    My hacky first implementation of the HotSpot method:
    "Identifying Informative Gene Modules Across Modalities of Single Cell Genomics"

    David DeTomaso, Nir Yosef
    https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    pvalues are crude bonferroni corrected

    """

    organism = adata.uns['organism']
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    nbrs_gex_clusters = setup_fake_nbrs_from_clusters( clusters_gex )
    nbrs_tcr_clusters = setup_fake_nbrs_from_clusters( clusters_tcr )
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
    comparisons = [['gex', nbrs_gex, 'graph', nbrs_tcr],
                   ['gex', nbrs_gex, 'clust', nbrs_tcr_clusters],
                   ['tcr', nbrs_tcr, 'graph', nbrs_gex],
                   ['tcr', nbrs_tcr, 'clust', nbrs_gex_clusters]]

    for feature_nbr_tag, feature_nbrs, graph_tag, graph_nbrs in comparisons:
        if graph_tag == 'clust' and not also_use_cluster_graphs:
            continue

        rows = [] # list of (1,num_clones) csr_matrix'es
        for ii in range(num_clones):
            ii_nbrs = feature_nbrs[ii]
            num_nbrs = len(ii_nbrs)
            data = np.full((num_nbrs,),1.0)
            row_ind = np.full((num_nbrs,),0).astype(int)
            new_row = sps.csr_matrix((data, (row_ind, ii_nbrs)),
                                     shape=(1,num_clones))
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
        df = correlations.find_hotspot_features(X, graph_nbrs, feature_names, pval_threshold)#, verbose=True)

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
        df = pd.concat(dfl, ignore_index=True)
        return df.drop(columns=['feature']) # feature is just str(clone_index), ie superfluous
    else:
        return pd.DataFrame()


def find_gex_cluster_degs(
        adata,
        outfile_prefix,
        min_cluster_size = 3,
):
    if len(set(adata.obs['clusters_gex'])) <= 1:
        print('conga.devel.find_gex_cluster_degs:: too few clusters')
        return

    # look at differentially expressed genes in gex clusters
    obs_tag = 'clusters_gex_for_degs'
    adata.obs[obs_tag] = [str(x) for x in adata.obs['clusters_gex']]
    key_added = 'degs_for_gex_clusters'
    rank_method = 'wilcoxon'
    good_clusters = sorted(
        [x for x,y in Counter(adata.obs[obs_tag]).most_common()
         if y>=min_cluster_size],
        key=int)
    sc.tl.rank_genes_groups(
        adata,
        groupby=obs_tag,
        method=rank_method,
        groups=good_clusters,
        reference='rest',
        key_added=key_added,
    )
    n_genes = 25
    sc.pl.rank_genes_groups(
        adata,
        n_genes=n_genes,
        sharey=False,
        show=False,
        key=key_added,
    )
    pngfile = outfile_prefix+'_gex_cluster_degs.png'
    plt.savefig(pngfile, bbox_inches="tight")
    print('made:', pngfile)


    new_rank_genes_genes, var_group_positions, var_group_labels = [],[],[]
    allow_gene_repeats = False
    min_rank_genes_log2fold_change = 1.0
    max_rank_genes_pval_adj=0.05
    n_genes_for_plotting = 5

    for group in good_clusters:
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
                  .format( log2fold, pval_adj,
                           adata.uns[key_added]['scores'][group][igene],
                           gene, group ) )
            my_genes.append( gene )
        if my_genes:
            var_group_positions.append((len(new_rank_genes_genes),
                                        len(new_rank_genes_genes)+
                                        len(my_genes)-1 ) )
            var_group_labels.append( group )
            new_rank_genes_genes.extend( my_genes )

    if new_rank_genes_genes:
        sc.pl.stacked_violin(
            adata, var_names = new_rank_genes_genes, groupby=obs_tag,
            figsize=(10,n_genes_for_plotting*10),
            use_raw = True, stripplot=True, show=False, swap_axes=True,
            var_group_positions = var_group_positions,
            var_group_labels = var_group_labels,
            var_group_rotation = 1.0 )
        pngfile = outfile_prefix+'_gex_cluster_degs_violin.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:',pngfile)

        sc.pl.dotplot(
            adata, var_names=new_rank_genes_genes, groupby=obs_tag, show=False,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions)
        pngfile = outfile_prefix+'_gex_cluster_degs_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

        # this plot_scatter seems to have moved in scanpy; need to update
        #sc.pl._tools.plot_scatter(adata, 'gex_2d', ncols = 6,
        #                          color = new_rank_genes_genes, show=False,
        #                          use_raw = True, s=40)
        #pngfile = args.outfile_prefix+'_gex_cluster_degs_tsne.png'
        #plt.savefig(pngfile, bbox_inches="tight")
        #print('made:', pngfile)


    if adata.uns['organism'] == 'human_ig':
        # list of B cell marker genes from "Human germinal centres engage
        # memory and naive B cells after influenza vaccination"
        # Turner...Ellebedy, Nature 2020:
        # https://doi.org/10.1038/s41586-020-2711-0
        # note that they say activated B cells are distinguished by lack of CR2
        genes_lines = """GC-Bs BCL6, RGS13, MEF2B, STMN1, ELL3, SERPINA9
        PBs XBP1, IRF4, SEC11C, FKBP11, JCHAIN, PRDM1
        naive TCL1A, IL4R, CCR7, IGHM, IGHD
        act-Bs TBX21, FCRL5, ITGAX, NKG7, ZEB2, CR2
        rest TNFRSF13B, CD27, CD24
        misc IGHA1 IGHA2 IGHG1 IGHG2 IGHG3 IGHG4 IGHE"""\
            .replace(',',' ').split('\n')
        genes, var_group_positions, var_group_labels = [], [], []
        for line in genes_lines:
            my_genes = [x for x in line.split()[1:] if x in adata.raw.var_names]
            print(len(my_genes), line.split())
            if my_genes:
                var_group_positions.append((len(genes),
                                            len(genes)+len(my_genes)-1))
                var_group_labels.append( line.split()[0])
                genes.extend(my_genes)
        sc.pl.dotplot(
            adata, var_names=genes, groupby=obs_tag, show=False,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions)
        pngfile = outfile_prefix+'_gex_cluster_bcell_genes_dotplot.png'
        plt.savefig(pngfile, bbox_inches="tight")
        print('made:', pngfile)

    # show some of our marker genes
    organism = adata.uns['organism']
    genes = (plotting.default_logo_genes[organism] +
             plotting.default_gex_header_genes[organism])
    genes = sorted(set(x for x in genes if x in adata.raw.var_names))
    sc.pl.dotplot(adata, var_names=genes, groupby=obs_tag, show=False)
    pngfile = outfile_prefix+'_gex_cluster_marker_genes_dotplot.png'
    plt.savefig(pngfile, bbox_inches="tight")
    print('made:', pngfile)




def find_hotspot_nbrhoods(
        adata,
        all_nbrs,
        outfile_prefix,
        make_hotspot_nbrhood_logos=False,
        min_cluster_size = 5,
        min_cluster_size_fraction = 0.001,
):

    # My hacky and probably buggy first implementation of the HotSpot method:
    #
    # "Identifying Informative Gene Modules Across Modalities of
    #  Single Cell Genomics"
    # David DeTomaso, Nir Yosef
    # https://www.biorxiv.org/content/10.1101/2020.02.06.937805v1

    nbr_fracs = sorted(all_nbrs.keys())

    all_results = []
    for nbr_frac in nbr_fracs:
        nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]
        print('find_hotspot_nbrhoods for nbr_frac', nbr_frac)
        results = correlations.find_hotspot_nbrhoods(
            adata, nbrs_gex, nbrs_tcr, pval_threshold=1.0,
            also_use_cluster_graphs=False)

        # the feature_type column is already set in nbrhood_results
        #  to {tcr/gex}_nbrs_vs_graph
        results['nbr_frac'] = nbr_frac
        all_results.append(results)

    # make a plot summarizing the hotspot nbrhood pvals and also
    #  save them to a file
    if all_results:
        nbrhood_results = pd.concat(all_results, ignore_index=True)

        tcrs = preprocess.retrieve_tcrs_from_adata(adata)

        outfile = '{}_hotspot_nbrhoods.tsv'.format(outfile_prefix)
        for iab, ivj in [ (x,y) for x in range(2) for y in range(3) ]:
            key = [ 'va ja cdr3a'.split(), 'vb jb cdr3b'.split()][iab][ivj]
            nbrhood_results[key] = [tcrs[x.clone_index][iab][ivj]
                                    for x in nbrhood_results.itertuples()]
        print('making:', outfile)
        nbrhood_results.to_csv(outfile, sep='\t', index=False)

        num_clones = adata.shape[0]
        nbrhood_pvals = {'gex':np.full((num_clones,), num_clones).astype(float),
                         'tcr':np.full((num_clones,), num_clones).astype(float)}
        for l in nbrhood_results.itertuples():
            assert l.feature_type[3:] == '_nbrs_vs_graph'
            tag = l.feature_type[:3]
            nbrhood_pvals[tag][l.clone_index] = min(
                l.pvalue_adj, nbrhood_pvals[tag][l.clone_index])

        plt.figure(figsize=(12,6))
        for icol, tag in enumerate(['gex','tcr']):
            plt.subplot(1,2,icol+1)
            log10_pvals = np.log10(np.maximum(1e-100, nbrhood_pvals[tag]))
            colors = np.sqrt(np.maximum(0.0, -1*log10_pvals))
            reorder = np.argsort(colors)
            # same umap as feature nbr-type
            xy = adata.obsm['X_{}_2d'.format(tag)]
            vmax = np.sqrt(-1*np.log10(1e-5))
            plt.scatter(xy[reorder,0], xy[reorder,1], c=colors[reorder],
                        vmin=0, vmax=vmax)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.xlabel('{} UMAP1'.format(tag))
            plt.ylabel('{} UMAP2'.format(tag))
            plt.title('{} hotspot nbrhood pvalues'.format(tag))
        plt.tight_layout()
        pngfile = '{}_hotspot_nbrhoods.png'.format(outfile_prefix)
        print('making:', pngfile)
        plt.savefig(pngfile)

        if make_hotspot_nbrhood_logos:
            nbrs_gex, nbrs_tcr = all_nbrs[ max(nbr_fracs) ]
            min_cluster_size = max(
                min_cluster_size,
                int(0.5 + min_cluster_size_fraction * num_clones))
            plotting.make_hotspot_nbrhood_logo_figures(
                adata, nbrs_gex, nbrs_tcr, nbrhood_results,
                min_cluster_size, outfile_prefix,
                pvalue_threshold=1.0)


def make_hotspot_nbrhood_logo_figures(
        adata,
        nbrs_gex,
        nbrs_tcr,
        hotspot_nbrhood_results_df,
        min_cluster_size,
        outfile_prefix,
        pvalue_threshold=1.0,
        **kwargs # passed to make_logo_plots
):
    ''' this is highly experimental and not really used, fyi
    '''
    #unpacking
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    num_clones = adata.shape[0]

    #
    nbrhood_pvals = {'gex':np.full((num_clones,), num_clones).astype(float),
                     'tcr':np.full((num_clones,), num_clones).astype(float)}

    for l in hotspot_nbrhood_results_df.itertuples():
        assert l.feature_type[3:] == '_nbrs_vs_graph'
        tag = l.feature_type[:3]
        nbrhood_pvals[tag][l.clone_index] = min(
            l.pvalue_adj, nbrhood_pvals[tag][l.clone_index])

    fake_clusters = np.full( (num_clones,), 0)

    for tag, pvals in nbrhood_pvals.items():
        # are there any clusters with enough significant pvals??
        if tag == 'gex':
            fake_clusters_gex = clusters_gex
            fake_clusters_tcr = fake_clusters
        else:
            fake_clusters_gex = fake_clusters
            fake_clusters_tcr = clusters_tcr

        pngfile = f'{outfile_prefix}_hotspot_nbrhood_{tag}_clusters.png'

        if 'rank_genes_uns_tag' not in kwargs:
            kwargs['rank_genes_uns_tag']= f'rg_hotspot_{tag}_nbrhood_biclusters'

        plotting.make_cluster_logo_plots_figure(
            adata, pvals, pvalue_threshold, fake_clusters_gex,
            fake_clusters_tcr,
            nbrs_gex, nbrs_tcr, min_cluster_size, pngfile,
            **kwargs)

def make_batch_bias_plots(
        adata,
        nbrhood_results_df, # generated by correlations.py: find_batch_biases
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size_for_logos,
        pvalue_threshold_for_logos,
        outfile_prefix,
        max_color_pvalue = 1e-9, # in the UMAP figure; no brighter beyond there
):
    num_clones = adata.shape[0]
    fake_clusters_gex = np.zeros((num_clones,)).astype(int)
    fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
    batch_bias_pvals = np.full( (num_clones,), num_clones).astype(float)
    for l in nbrhood_results_df.itertuples():
        batch_bias_pvals[ l.clone_index] = min(
            l.pvalue_adj, batch_bias_pvals[l.clone_index])
        fake_clusters_tcr[l.clone_index] = l.cluster_group

    pngfile = f'{outfile_prefix}_batch_bias_logos.png'

    plotting.make_cluster_logo_plots_figure(
        adata, batch_bias_pvals, pvalue_threshold_for_logos, fake_clusters_gex,
        fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile,
        rank_genes_uns_tag='rg_batch_bias_biclusters')

    # make umaps colored by batch_bias pvals
    plt.figure(figsize=(12,6))
    for icol, tag in enumerate(['gex','tcr']):
        plt.subplot(1,2,icol+1)
        colors = np.sqrt( np.maximum(
            0.0, -1*np.log10(np.maximum(1e-100, batch_bias_pvals))))
        reorder = np.argsort(colors)
        xy = adata.obsm['X_{}_2d'.format(tag)] # same umap as feature nbr-type
        vmax = np.sqrt(-1*np.log10(max_color_pvalue))
        #vmax = np.sqrt(-1*np.log10(1e-5))
        plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder],
                     vmin=0, vmax=vmax)
        plt.colorbar()
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel('{} UMAP1'.format(tag.upper()))
        plt.ylabel('{} UMAP2'.format(tag.upper()))
        plt.title('TCR batch_bias pvalues')
    plt.tight_layout()
    pngfile = '{}_batch_bias.png'.format(outfile_prefix)
    print('making:', pngfile)
    plt.savefig(pngfile)

def analyze_proteins(
        adata,
        outfile_prefix,
        exclude_protein_prefixes = ['HTO', 'TCRV', 'TCRv'],
        n_components_prot=20,
        n_components_gex=40,

):
    ''' run pca on the protein data

    '''
    from sklearn.metrics import pairwise_distances
    from scipy.stats import describe

    # tmp hacking:
    num_clones = adata.shape[0]
    # if 'group' not in adata.uns['batch_keys'] or adata.obsm['group'].shape != (num_clones,4):
    #     print('need proper group info:', 'group' in adata.uns['batch_keys'])
    #     return

    ft_varname = util.get_feature_types_varname(adata)
    prot_mask = np.array(adata.raw.var[ft_varname] == 'Antibody Capture')
    for p in exclude_protein_prefixes:
        prot_mask &= ~adata.raw.var_names.str.startswith(p)
    print('used_protein_features:', list(adata.raw.var_names[prot_mask]))
    X = adata.raw.X.tocsc()[:,prot_mask].toarray()
    print(type(X), X.shape)
    mn = np.mean(X, axis=0)
    #print('mn:', mn)
    sd = np.std(X, axis=0)
    #print('sd:', sd)
    X = (X-mn[np.newaxis,:])/sd[np.newaxis,:]
    pca = PCA()
    X_pca = pca.fit_transform(X)
    print('pca_variance:', pca.explained_variance_ratio_[:n_components_prot+1])
    X_pca_prot = X_pca[:,:n_components_prot]
    X_pca_gex  = adata.obsm['X_pca_gex'][:,:n_components_gex]
    nrandom = 1000
    inds = np.random.permutation(adata.shape[0])[:nrandom]
    D_prot = pairwise_distances(X_pca_prot[inds,:])
    D_gex = pairwise_distances(X_pca_gex[inds,:])
    print('gex_dists: ', describe(D_gex.ravel()))
    print('prot_dists:', describe(D_prot.ravel()))
    print('dist_ratio:', np.mean(D_prot.ravel()) / np.mean(D_gex.ravel()))

    X_pca_combo = np.hstack([X_pca_gex, X_pca_prot])
    print('X_pca_combo:', X_pca_combo.shape)

    for tag, X_pca in [['prot', X_pca_prot], ['combo', X_pca_combo]]:

        adata.obsm['X_pca'] = X_pca
        adata.obsm[f'X_pca_{tag}'] = np.copy(X_pca) # force copy
        print('neighbors:', tag, X_pca.shape)
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=X_pca.shape[1])

        print('umap:', tag)
        sc.tl.umap(adata)
        adata.obsm['X_umap_'+tag] = adata.obsm['X_umap']
        cluster_key_added = 'louvain_'+tag
        resolution = 1.0
        print('louvain:', tag)
        sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
        adata.obs['clusters_'+tag] = np.copy(adata.obs[cluster_key_added]).astype(int)
        adata.obsm['X_{}_2d'.format(tag)] = adata.obsm['X_umap_'+tag]




    xy_tags = ['gex','prot','combo']
    special_genes = ['ZBTB16','RTKN2','CD79A']
    trbv = 'TRBV12-5'
    genes = ['CD45RA_p']+special_genes+['min_special',trbv]

    if 'group' in adata.obs_keys():
        other = ['CD4','CD8','group_0','group_1','group_2','group_3']

        group_counts = adata.obsm['group'].astype(float)
        group_fracs = group_counts/(np.sum(group_counts,axis=1)[:,np.newaxis])
        assert np.max(group_fracs)<1.001
        assert group_fracs.shape == (num_clones, 4)
    else:
        other = ['CD4','CD8']

    nrows, ncols, plotno = len(xy_tags), 1+len(genes)+len(other), 0
    plt.figure(figsize=(ncols*2.5, nrows*3))

    var_names = list(adata.raw.var_names)
    for xy_tag in xy_tags:
        xy = adata.obsm[f'X_{xy_tag}_2d']

        # calculate nbrs for averaging
        obsm_tag = f'X_pca_{xy_tag}'
        assert obsm_tag in adata.obsm_keys()
        num_nbrs = max(1, min(num_clones//10, max(20, num_clones//500)))
        nbr_frac = (num_nbrs+0.1)/num_clones
        all_nbrs = preprocess.calc_nbrs(adata, [nbr_frac], obsm_tag_gex=obsm_tag, obsm_tag_tcr=None)
        nbrs = all_nbrs[nbr_frac][0]
        assert num_nbrs == nbrs.shape[1]
        all_nbr_avg_colors = {}
        for colortag in [xy_tag] + other + genes:
            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            cmap = 'viridis'
            title = colortag
            nbr_avg = False
            sort_order = True
            if colortag == 'min_special':
                colors = np.hstack( [all_nbr_avg_colors[x][:,np.newaxis] for x in special_genes])
                assert colors.shape == (num_clones, len(special_genes))
                colors = np.min(colors, axis=1)
            elif colortag == trbv:
                colors = np.array(adata.obs['vb'].str.startswith(trbv)).astype(float)
                nbr_avg = True
            elif colortag == 'CD4':
                colors = adata.raw.X[:,var_names.index('CD4_p')].toarray()[:,0]
            elif colortag == 'CD8':
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8_p')].toarray()[:,0]+
                                adata.raw.X[:,var_names.index('CD8a_p')].toarray()[:,0])
            elif colortag in genes:
                colors = adata.raw.X[:,var_names.index(colortag)].toarray()[:,0]
                nbr_avg = True
            elif colortag in xy_tags:
                colors = adata.obs[f'clusters_{colortag}']
                cmap = 'tab20'
                title = f'{colortag} clusters ({np.max(colors)+1})'
                sort_order = False
            else:
                assert colortag.startswith('group_')
                g = int(colortag.split('_')[1])
                colors = group_fracs[:,g]
                nbr_avg = True

            if nbr_avg:
                colors = (colors + colors[nbrs].sum(axis=1))/(num_nbrs+1)
                title = f'{title} nbr-avg-{num_nbrs}'
                all_nbr_avg_colors[colortag] = colors

            if sort_order:
                reorder = np.argsort(colors)
            else:
                reorder = np.arange(num_clones)

            plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
            plt.xticks([],[])
            plt.yticks([],[])
            #plt.xlabel(xy_tag)
            plt.ylabel(xy_tag)
            plt.title(title)
            plt.colorbar()
    pngfile = outfile_prefix+'_prot_combo_special_panels.png'
    print('making:', pngfile)
    plt.tight_layout()
    plt.savefig(pngfile)

def analyze_special_genes(
        adata,
        outfile_prefix,
        special_genes = ['ZBTB16', 'RTKN2', 'CD79A'],
        obs_tag_clusters = 'clusters_gex',
        obsm_tag_nbrs = 'X_pca_gex',
        obsm_tag_xy = 'X_gex_2d',
        trbvs = ['TRBV12-5', 'TRBV4-1'],
        other_genes = ['KLRB1'],
        nbrs = None,
        num_nbrs = None,

):
    ''' analyze special genes

    '''
    #from sklearn.metrics import pairwise_distances
    from scipy.stats import linregress, mannwhitneyu

    assert obsm_tag_xy in adata.obsm_keys()
    assert obsm_tag_nbrs in adata.obsm_keys()

    num_clones = adata.shape[0]

    genes = ['clusters', 'CD4', 'CD8'] + other_genes + special_genes + ['min_special']+trbvs

    ncols = 4
    nrows = (len(genes)+len(trbvs)-1)//ncols + 1
    plotno=0

    plt.figure(figsize=(ncols*3, nrows*3))

    var_names = list(adata.raw.var_names)
    xy = adata.obsm[obsm_tag_xy]

    # calculate nbrs for averaging
    if num_nbrs is None:
        if nbrs is None:
            num_nbrs = max(1, min(num_clones//10, max(20, num_clones//500)))
        else:
            num_nbrs = nbrs.shape[1]

    if num_nbrs and nbrs is None:
        print('computing nbrs:', adata.shape, num_nbrs)
        nbr_frac = (num_nbrs+0.1)/num_clones
        all_nbrs = preprocess.calc_nbrs(adata, [nbr_frac], obsm_tag_gex=obsm_tag_nbrs, obsm_tag_tcr=None)
        nbrs = all_nbrs[nbr_frac][0]
        assert num_nbrs == nbrs.shape[1]

    all_nbr_avg_colors = {}
    for colortag in genes:
        plotno += 1
        plt.subplot(nrows, ncols, plotno)
        cmap = 'viridis'
        title = colortag
        nbr_avg = False
        sort_order = True
        if colortag == 'min_special':
            colors = np.hstack( [all_nbr_avg_colors[x][:,np.newaxis] for x in special_genes])
            assert colors.shape == (num_clones, len(special_genes))
            colors = np.min(colors, axis=1)
            all_nbr_avg_colors[colortag] = colors # already averaged...
        elif colortag == 'clusters':
            colors = adata.obs[obs_tag_clusters]
            cmap = 'tab20'
            title = f'clusters ({np.max(colors)+1})'
            sort_order = False
        elif colortag in trbvs:
            colors = np.array(adata.obs['vb'].str.startswith(colortag)).astype(float)
            nbr_avg = True
        elif colortag == 'CD4':
            if 'CD4_p' in var_names:
                title = 'CD4_prot'
                colors = adata.raw.X[:,var_names.index('CD4_p')].toarray()[:,0]
            elif 'CD4' in var_names:
                title = 'CD4_gex'
                colors = adata.raw.X[:,var_names.index('CD4')].toarray()[:,0]
            else:
                colors = np.zeros((num_clones,))
        elif colortag == 'CD8':
            if 'CD8_p' in var_names:
                title = 'CD8_prot'
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8_p')].toarray()[:,0]+
                                 adata.raw.X[:,var_names.index('CD8a_p')].toarray()[:,0])
            else:
                title = 'CD8_gex'
                colors = 0.5 * ( adata.raw.X[:,var_names.index('CD8A')].toarray()[:,0]+
                                 adata.raw.X[:,var_names.index('CD8B')].toarray()[:,0])
        else:
            colors = adata.raw.X[:,var_names.index(colortag)].toarray()[:,0]
            nbr_avg = True

        if nbr_avg:
            if num_nbrs:
                colors = (colors + colors[nbrs].sum(axis=1))/(num_nbrs+1)
            title = f'{title} nbr-avg-{num_nbrs}'
            all_nbr_avg_colors[colortag] = colors

        if sort_order:
            reorder = np.argsort(colors)
        else:
            reorder = np.arange(num_clones)

        plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.xlabel(xy_tag)
        #plt.ylabel(xy_tag)
        plt.title(title)
        plt.colorbar()

    for trbv in trbvs:
        xvals = all_nbr_avg_colors[trbv]
        yvals = all_nbr_avg_colors['min_special']
        slope, intercept, r, p, err = linregress(xvals, yvals)
        mask = np.array(adata.obs['vb'].str.startswith(trbv))
        _,mwu_p = mannwhitneyu(yvals[mask], yvals[~mask],
                               **util.mannwhitneyu_kwargs)
        print(f'trbv_vs_special: {trbv} {r:.2f} {p:.2e} {mwu_p:.2e}')
        plotno+= 1
        plt.subplot(nrows, ncols, plotno)
        plt.scatter(xvals, yvals, s=5)
        plt.title(f'{trbv} r: {r:.2f}\nP: {p:.1e} {mwu_p:.1e}')
        plt.xlabel(trbv)
        plt.ylabel('min_special')

    pngfile = outfile_prefix+'_special_genes.png'
    print('making:', pngfile)
    plt.tight_layout()
    plt.savefig(pngfile, dpi=200)
    plt.close()

    return nbrs # may be useful for caller?


def analyze_CD4_CD8(
        adata,
        nbrs_gex,
        outfile_prefix,
        cd4_genes_default = ['CD4'],
        cd8_genes_default = ['CD8A', 'CD8B'],
        cd4_proteins_default = ['CD4_p'],
        cd8_proteins_default = ['CD8_p', 'CD8a_p'],
):
    import matplotlib.pyplot as plt
    if adata.uns['organism'] != 'human':
        # only really makes sense for human tcr AB
        print('analyze_CD4_CD8 incompatible organism:', adata.uns['organism'])
        return

    clusters_gex = np.array(adata.obs['clusters_gex'])
    num_clones = adata.shape[0]

    for ptag in ['gex','prot']:
        if ptag == 'gex':
            cd4_genes = [ x for x in cd4_genes_default if x in adata.raw.var_names]
            cd8_genes = [ x for x in cd8_genes_default if x in adata.raw.var_names]
        else:
            cd4_genes = [ x for x in cd4_proteins_default if x in adata.raw.var_names]
            cd8_genes = [ x for x in cd8_proteins_default if x in adata.raw.var_names]

        if not cd4_genes or not cd8_genes:
            continue

        all_gex = {}

        for gene in cd4_genes + cd8_genes:
            index = list(adata.raw.var_names).index(gene)
            all_gex[gene] = adata.raw.X[:,index].toarray()[:,0]

        # look at neighborhoods
        xvals = np.copy(all_gex[cd8_genes[0]])
        for g in cd8_genes[1:]:
            xvals += all_gex[g]
        xvals/= len(cd8_genes)

        yvals = np.copy(all_gex[cd4_genes[0]])
        for g in cd4_genes[1:]:
            yvals += all_gex[g]
        yvals/= len(cd4_genes)

        num_nbrs = nbrs_gex.shape[1]

        xvals_nbr_avged = (xvals + xvals[nbrs_gex].sum(axis=1))/(num_nbrs+1)
        yvals_nbr_avged = (yvals + yvals[nbrs_gex].sum(axis=1))/(num_nbrs+1)

        npanels = 2 + 3 + len(cd4_genes) + len(cd8_genes)
        ncols = 4
        nrows = (npanels-1)/ncols + 1

        plt.figure(figsize=(ncols*3, nrows*3))
        plt.subplot(nrows, ncols, 1)
        plt.scatter(xvals_nbr_avged, yvals_nbr_avged, alpha=0.25, c=clusters_gex, cmap='tab20')
        if ptag=='prot':
            plt.xlim((0, 5))
            plt.ylim((0, 5))
        else:
            plt.xlim((0, plt.xlim()[1]))
            plt.ylim((0, plt.ylim()[1]))
        plt.xlabel('+'.join(cd8_genes))
        plt.ylabel('+'.join(cd4_genes))


        plt.subplot(nrows, ncols, 2)
        num_clusters = np.max(clusters_gex)+1
        vals = []
        for c in range(num_clusters):
            mask = clusters_gex==c
            x = np.sum(xvals[mask])/np.sum(mask)
            y = np.sum(yvals[mask])/np.sum(mask)
            vals.append((x,y,c))
        plt.scatter([x[0] for x in vals], [x[1] for x in vals], c=[x[2] for x in vals], cmap='tab20', s=50)
        for x,y,c in vals:
            plt.text(x,y,str(c))
        plt.xlabel('+'.join(cd8_genes))
        plt.ylabel('+'.join(cd4_genes))
        if ptag=='prot':
            plt.xlim((0, 5))
            plt.ylim((0, 5))
            plt.plot([0,5],[0,5],':k')
        else:
            plt.xlim((0, plt.xlim()[1]))
            plt.ylim((0, plt.ylim()[1]))

        plotno=2
        xy = adata.obsm['X_gex_2d']
        for color in ['cluster','n_genes','n_counts']+cd4_genes+cd8_genes:
            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            if color == 'cluster':
                colors = clusters_gex
                cmap='tab20'
                reorder = np.arange(len(colors))
            elif color in ['n_genes','n_counts']:
                colors = np.array(adata.obs[color])
                cmap = 'viridis'
                reorder = np.argsort(colors)
            else:
                colors = all_gex[color]
                cmap = 'viridis'
                reorder = np.argsort(colors)
            plt.scatter(xy[reorder,0], xy[reorder,1], s=5, c=colors[reorder], cmap=cmap)
            plt.colorbar()
            plt.title(color)

        pngfile = f'{outfile_prefix}_cd4_cd8_{ptag}.png'
        plt.savefig(pngfile)
        print('making:', pngfile)



def filter_cells_by_ribo_norm(adata):
    ''' returns  new filtered adata
    will normalize_and_log_the_raw_matrix if not already done
    '''
    from scipy.stats import gaussian_kde

    preprocess.normalize_and_log_the_raw_matrix(adata)

    X_norm = adata.raw.X

    X_norm_ribo = []
    var_names = [x.upper() for x in adata.raw.var_names ] ## NOTE UPPER-- help for mouse
    for ig,g in enumerate(var_names):
        if g.startswith('RP') and len(g)>=4 and g[2] in 'SL' and g[3].isdigit():
            X_norm_ribo.append( X_norm[:,ig].toarray() )

    X_norm_ribo = np.sum(np.hstack( X_norm_ribo ), axis=1 )/len(X_norm_ribo)

    dens = gaussian_kde(X_norm_ribo)
    # look for two maxima and a minima in between
    # or look for the closest minimum to 2.0
    step = 0.01
    x=2.0
    while True:
        y0 = dens(x)[0]
        yl = dens(x-step)[0]
        yr = dens(x+step)[0]
        print('filter_cells_by_ribo_norm:: minfind:',x,y0,yl,yr)
        if yl < yr:
            if yl < y0:
                x = x-step
            else:
                # found a local min
                ribo_norm_threshold=x
                break
        else:
            if yr < y0:
                x = x+step
            else:
                ribo_norm_threshold=x
                break
    print('filter_cells_by_ribo_norm:: ribo_norm localmin:',ribo_norm_threshold)

    mask = (X_norm_ribo>ribo_norm_threshold)
    print('filter_cells_by_ribo_norm::', adata.shape, 'downto:', np.sum(mask))
    sys.stdout.flush()
    adata = adata[mask,:].copy()
    assert np.sum(mask) == adata.shape[0]
    return adata


def _get_pvalue_from_rvalue(r, n):
    ''' hacky helper stolen from scipy.stats.linregress
    alternative='two-sided'
    '''
    from scipy import stats
    TINY = 1.0e-20
    df = n - 2  # Number of degrees of freedom
    # n-2 degrees of freedom because 2 has been used up
    # to estimate the mean and standard deviation
    t = r * np.sqrt(df / ((1.0 - r + TINY)*(1.0 + r + TINY)))
    #t, prob = _ttest_finish(df, t, alternative)
    prob = 2*stats.t.cdf(-abs(t), df)
    return prob



def run_genes_versus_features(
        adata,
        tcr_features = None,
        verbose=True,
):
    if tcr_features is None:
        tcr_features = (tcr_scoring.all_tcr_scorenames+
                        [x+'_frac' for x in tcr_scoring.amino_acids])


    print('making tcr scoretable:', adata.shape[0], len(tcr_features))
    organism = adata.uns['organism']
    num_clones = adata.shape[0]
    num_features = len(tcr_features)

    # rare genes seem to give spurious correlations
    min_cells_per_gene = max(5, min(50, 0.001*num_clones))

    tcr_scoretable = tcr_scoring.make_tcr_score_table(
        adata, tcr_features, verbose=True)

    # run with variable genes and with all genes
    # pay a factor of 2 in pvalues for that
    variable_genes = list(adata.var_names)
    all_genes = list(adata.raw.var_names)

    Y = tcr_scoretable.T.copy()
    assert Y.shape == (len(tcr_features), num_clones)

    Y_mean = Y.mean(axis=1)
    Y_std = Y.std(axis=1)
    Y = (Y - Y_mean[:,None])/Y_std[:,None]
    low_var_mask = np.abs(Y_std)<1e-6
    print('run_features_versus_features:: dropping tcr features with',
          'low variance:', [x for x,y in zip(tcr_features, low_var_mask)
                            if y])
    Y[low_var_mask,:] = 0. # dont want nans
    print('Y nans:', np.sum(np.isnan(Y)))

    dfl = []
    for is_variable in [True,False]:
        if is_variable:
            genes = variable_genes

            X = adata.X.T.copy()
            X = (X - X.mean(axis=1)[:,None])/X.std(axis=1)[:,None]

            C = X@Y.T/X.shape[1]

        else:
            genes = all_genes

            X = adata.raw.X.T
            print('dotting... sparse X and Y', X.shape, Y.shape)
            D = X.dot(Y.T)

            print('mean/std of sparse X...', X.shape)
            X_sq = X.multiply(X)

            X_mean = X.mean(axis=1).A1
            X_sq_mean = X_sq.mean(axis=1).A1
            X_mean_sq = X_mean**2
            X_var = X_sq_mean - X_mean_sq
            X_std = X_var**0.5

            X_low_var_mask = np.abs(X_std)<1e-6
            print('run_genes_versus_features:: low variance genes:',
                  [x for x,y in zip(genes, X_low_var_mask) if y])
            X_std = np.maximum(X_std, 1e-6)

            print('correlations...')
            Y_sum = Y.sum(axis=1)
            C = (D - np.outer(X_mean, Y_sum))/(num_clones * X_std[:,None])

        print('C nans:', np.sum(np.isnan(C)))
        assert np.sum(np.isnan(C))==0

        num_genes = len(genes)
        assert C.shape == (num_genes, num_features)
        pval_rescale = 2*num_genes*(num_features-np.sum(low_var_mask))

        inds = np.argsort(C.ravel())[::-1]

        for r in range(2):
            feature_counts = Counter()
            rinds = inds if r==0 else inds[::-1]
            for ind in rinds:
                top = np.unravel_index(ind, C.shape)
                rvalue = C[top]
                gene = genes[top[0]]
                if util.is_vdj_gene(gene, organism):
                    continue
                if not is_variable and gene in variable_genes:
                    continue
                feature = tcr_features[top[1]]
                tfstar = 'TF' if (organism == 'human' and
                                  gene in human_tf_gene_names) else '  '

                if is_variable:
                    gene_vals = X[top[0],:]
                else:
                    gene_vals = X[top[0],:].toarray()[0]
                feature_vals = tcr_scoretable.T[top[1],:]
                gene_nonzero = np.sum(np.abs(gene_vals)>1e-3)
                feature_nonzero = np.sum(np.abs(feature_vals)>1e-3)
                if gene_nonzero < min_cells_per_gene:
                    continue
                #pval = pval_rescale*res.pvalue
                pval = pval_rescale*_get_pvalue_from_rvalue(rvalue, num_clones)

                # res = linregress(gene_vals, feature_vals)
                # print('equal?', res.pvalue*pval_rescale, pval, res.rvalue,
                #       rvalue, gene_vals.shape, feature_vals.shape,
                #       type(gene_vals), type(feature_vals),
                #       num_clones, gene_nonzero,
                #       feature_nonzero, gene, feature, r)

                if pval>1:
                    break
                feature_counts[feature] += 1
                dfl.append(dict(
                    gene=gene,
                    feature=feature,
                    is_variable_gene=is_variable,
                    pvalue_adj=pval,
                    rvalue=rvalue,
                    is_human_tf = (gene in human_tf_gene_names),
                    num_nonzero_gene=gene_nonzero,
                    num_nonzero_feature=feature_nonzero,
                ))
                vstar = 'V' if gene in variable_genes else ' '
                if verbose and feature_counts[feature] < 3:
                    print(f'{tfstar} {vstar} {C[top]:6.3f} {rvalue:6.3f}',
                          f'{pval:9.2e} {gene_nonzero:5d} {feature_nonzero:5d}',
                          f'{gene:10s} {feature:15s}')


    results = pd.DataFrame(dfl)
    results.sort_values('pvalue_adj', inplace=True)
    return results

def filter_sorted_gene_features(
        adata,
        features,
        max_features,
):
    ''' return new_features, duplicates

    where duplicates is a dict mapping from rep features to their cluster peers
    '''

    A = make_raw_feature_scores_table(
        features, ['gex']*len(features), adata, verbose=True).T

    assert A.shape == (len(features), adata.shape[0])

    print('calc correlations:', A.shape)
    start = time.time()
    C = 1-distance.squareform(
        distance.pdist(A, metric='correlation'), force='tomatrix')
    print('took:', time.time() - start)

    new_features, duplicates = filter_sorted_features_by_correlation(
        features, C, max_features)

    return new_features, duplicates

def filter_sorted_gex_features_in_dataframe(
        results,
        adata,
        max_features,
):
    ''' returns new dataframe, and duplicates dict

    results should have the columns 'feature' and 'feature_type'

    '''

    gex_features = [x for x,y in zip(results.feature, results.feature_type)
                    if y=='gex']

    new_features, duplicates = filter_sorted_gene_features(
        adata, gex_features, max_features)

    new_features = set(new_features)

    mask = [ (y=='tcr') or (x in new_features)
             for x,y in zip(results.feature, results.feature_type)]

    return results[mask].copy(), duplicates


def find_tcr_clumping_single_chain(
        tcrs,
        organism,
        tmpfile_prefix = 'tmp',
        radii = [6, 12, 24, 48],
        num_random_samples_multiplier = 100, # ie, 100 * num_clones
        pvalue_threshold = 1.0,
        verbose=True,
        preserve_vj_pairings = False,
        bg_tcrs = None, # usually better to leave this None
):
    ''' Returns a pandas dataframe with the following columns:
    - chain (A or B)
    - clone_index
    - nbr_radius
    - pvalue_adj
    - num_nbrs
    - expected_num_nbrs
    - raw_count
    - va, ja, cdr3a, vb, jb, cdr3b (ie, the 6 tcr cols for clone_index clone)
    - clumping_group: clonotypes within each other's significant nbr_radii are linked
    - clump_type: string, either 'global' or 'intra_gex_cluster' (latter only if also_find_clumps_within_gex_clusters=T)

    '''

    outprefix = f'{tmpfile_prefix}_{random.random()}_tcr_clumping'

    num_clones = len(tcrs)
    num_random_samples = num_random_samples_multiplier * num_clones

    if bg_tcrs is None:
        bg_tcrs = tcrs

    agroups, bgroups = preprocess.setup_tcr_groups_for_tcrs(tcrs)

    tmpfiles = [] # for cleanup

    ## compute background neighbor counts at the specified radii
    tcr_clumping.estimate_background_tcrdist_distributions(
        organism, tcrs, max(radii), num_random_samples,
        tmpfile_prefix=outprefix,
        preserve_vj_pairings=preserve_vj_pairings,
        save_unpaired_dists=True,
        tcrs_for_background_generation=bg_tcrs,
        #nocleanup=True,
    )

    afile = f'{outprefix}_dists.txt_A.txt'
    bfile = f'{outprefix}_dists.txt_B.txt'
    if not (exists(afile) and exists(bfile)):
        print('ERROR find_tcr_clumping_single_chain::',
              'estimate_background_tcrdist_distributions failed')
        return pd.DataFrame()
    tmpfiles.extend([afile, bfile])

    acounts = np.cumsum(np.loadtxt(afile, dtype=int), axis=1)[:,radii]
    bcounts = np.cumsum(np.loadtxt(bfile, dtype=int), axis=1)[:,radii]

    # newcol = np.full((len(tcrs),1), num_random_samples)
    # acounts = np.hstack([acounts[:,radii], newcol])
    # bcounts = np.hstack([bcounts[:,radii], newcol])


    big_dfl = [] # will hold results for both chains as dataframes

    tcrs_file = outprefix +'_tcrs.tsv'
    pd.DataFrame({
        'va'   : [x[0][0] for x in tcrs],
        'cdr3a': [x[0][2] for x in tcrs],
        'vb'   : [x[1][0] for x in tcrs],
        'cdr3b': [x[1][2] for x in tcrs],
    }).to_csv(tcrs_file, sep='\t', index=False)
    tmpfiles.append(tcrs_file)

    for chain, chain_groups, bg_counts in [['A', agroups, acounts],
                                           ['B', bgroups, bcounts]]:

        # find neighbors in fg tcrs up to max(radii) ############################

        if os.name == 'posix':
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain')
        else:
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain.exe')


        db_filename = Path.joinpath(
            Path(util.path_to_tcrdist_cpp_db), f'tcrdist_info_{organism}.txt')

        tcrdist_threshold = max(radii)

        cmd = (f'{exe} -f {tcrs_file} -t {tcrdist_threshold} -d {db_filename}'
               f' -c {chain} -o {outprefix}')

        util.run_command(cmd, verbose=True)

        nbr_indices_filename = f'{outprefix}_nbr{tcrdist_threshold}_indices.txt'
        nbr_distances_filename = f'{outprefix}_nbr{tcrdist_threshold}_distances.txt'
        tmpfiles.extend([nbr_indices_filename, nbr_distances_filename])

        if not exists(nbr_indices_filename) or not exists(nbr_distances_filename):
            print('find_neighbors failed:', exists(nbr_indices_filename),
                  exists(nbr_distances_filename))
            exit(1)

        all_nbrs = []
        all_distances = []
        for group, line1, line2 in zip(chain_groups,
                                       open(nbr_indices_filename,'r'),
                                       open(nbr_distances_filename,'r')):
            nbrs = [int(x) for x in line1.split()]
            dists = [int(x) for x in line2.split()]
            assert len(nbrs) == len(dists)
            mask = [chain_groups[x] != group for x in nbrs]
            all_nbrs.append([x for x,m in zip(nbrs,mask) if m])
            all_distances.append([x for x,m in zip(dists,mask) if m])
        assert len(all_nbrs) == num_clones

        # use poisson to find nbrhoods with more tcrs than expected;
        #  have to handle agroups/bgroups
        dfl = []

        is_clumped = np.full((num_clones,), False)

        all_raw_pvalues = np.full((num_clones, len(radii)), 1.0)

        for ii in range(num_clones):
            ii_bg_counts = bg_counts[ii]
            ii_dists = all_distances[ii]
            for irad, radius in enumerate(radii):
                num_nbrs = sum(x<=radius for x in ii_dists)
                if num_nbrs<1:
                    continue
                max_nbrs = np.sum(chain_groups != chain_groups[ii])
                pval = hypergeom.sf(
                    num_nbrs-1, # total fg nbr tcrs (-1)
                    max_nbrs+num_random_samples, # total fg+bg tcrs
                    num_nbrs+ii_bg_counts[irad], # total nbr tcrs
                    max_nbrs) # total fg tcrs
                mu = max_nbrs * ii_bg_counts[irad]/num_random_samples
                all_raw_pvalues[ii, irad] = pval
                pval *= len(radii) * num_clones # simple multiple test correction
                #hgeom_pval *= len(radii) * num_clones # simple multiple test correction
                if pval <= pvalue_threshold:
                    is_clumped[ii] = True
                    raw_count = ii_bg_counts[irad]
                    if verbose:
                        atcr_str = ' '.join(tcrs[ii][0][:3])
                        btcr_str = ' '.join(tcrs[ii][1][:3])
                        print(f'tcr_nbrs_global: {num_nbrs:2d} {mu:9.6f}',
                              f'radius: {radius:2d} pval: {pval:9.1e}',
                              f'{raw_count:9.1f} tcr: {atcr_str} {btcr_str}')

                    dfl.append( OrderedDict(
                        chain=chain,
                        clone_index=ii,
                        nbr_radius=radius,
                        pvalue_adj=pval,
                        #hgeom_pvalue_adj=hgeom_pval,
                        num_nbrs=num_nbrs,
                        expected_num_nbrs=mu,
                        raw_count=raw_count,
                        va   =tcrs[ii][0][0],
                        ja   =tcrs[ii][0][1],
                        cdr3a=tcrs[ii][0][2],
                        vb   =tcrs[ii][1][0],
                        jb   =tcrs[ii][1][1],
                        cdr3b=tcrs[ii][1][2],
                    ))

        results_df = pd.DataFrame(dfl)
        if results_df.shape[0] == 0:
            big_dfl.append(results_df)
            continue # to next chain


        # compute FDR values in addition to the simple adjusted pvalues
        _, fdr_values, _, _ = multipletests(
            all_raw_pvalues.reshape(-1), alpha=0.05, method='fdr_bh')
        fdr_values = fdr_values.reshape((num_clones, len(radii))).min(axis=1)
        # right now we don't assign fdr values for intra-gex cluster clumping
        results_df['clonotype_fdr_value'] = [
            fdr_values[x] for x in results_df.clone_index]

        # identify groups of related hits?
        all_clumped_nbrs = {}
        for l in results_df.itertuples():
            ii = l.clone_index
            radius = l.nbr_radius
            clumped_nbrs = set(x for x,y in zip(all_nbrs[ii], all_distances[ii])
                               if y<= radius and is_clumped[x])
            clumped_nbrs.add(ii)
            if ii in all_clumped_nbrs:
                all_clumped_nbrs[ii] = all_clumped_nbrs[ii] | clumped_nbrs
            else:
                all_clumped_nbrs[ii] = clumped_nbrs


        clumped_inds = sorted(all_clumped_nbrs.keys())
        assert len(clumped_inds) == np.sum(is_clumped)

        # make nbrs symmetric
        for ii in clumped_inds:
            for nbr in all_clumped_nbrs[ii]:
                assert nbr in all_clumped_nbrs
                all_clumped_nbrs[nbr].add(ii)

        all_smallest_nbr = {}
        for ii in clumped_inds:
            all_smallest_nbr[ii] = min(all_clumped_nbrs[ii])

        while True:
            updated = False
            for ii in clumped_inds:
                nbr = all_smallest_nbr[ii]
                new_nbr = min(nbr, np.min([all_smallest_nbr[x]
                                           for x in all_clumped_nbrs[ii]]))
                if nbr != new_nbr:
                    all_smallest_nbr[ii] = new_nbr
                    updated = True
            if not updated:
                break

        # define clusters, choose cluster centers
        clusters = np.array([0]*num_clones) # 0 if not clumped

        cluster_number=0
        cluster_sizes = Counter()
        for ii in clumped_inds:
            nbr = all_smallest_nbr[ii]
            if ii==nbr:
                cluster_number += 1
                members = [ x for x,y in all_smallest_nbr.items() if y==nbr]
                clusters[members] = cluster_number
                cluster_sizes[cluster_number] = len(members)

        for ii, nbrs in all_clumped_nbrs.items():
            for nbr in nbrs:
                # confirm single-linkage clusters
                assert clusters[ii] == clusters[nbr]

        assert not np.any(clusters[is_clumped]==0)
        assert np.all(clusters[~is_clumped]==0)

        # reorder the clumping groups by size
        remap = {x[0]:i+1 for i,x in enumerate(cluster_sizes.most_common())}
        remap[0] = 0

        results_df['clumping_group'] = [ remap[clusters[x.clone_index]]
                                         for x in results_df.itertuples()]

        results_df.sort_values('pvalue_adj', inplace=True)

        big_dfl.append(results_df)

    # cleanup the tmpfiles
    for tmpfile in tmpfiles:
        if exists(tmpfile):
            os.remove(tmpfile)

    results = pd.concat(big_dfl)
    return results


def make_single_chain_clumping_logos(
        results_df, # generated by find_tcr_clumping_single_chain
        adata,
        nbrs_gex,
        nbrs_tcr,
        outfile_prefix,
        min_cluster_size_for_logos=3,
        pvalue_threshold_for_logos=1.0, #pvalues are crude bonferroni corrected
        max_color_pvalue = 1e-16, # in the UMAP figure; no brighter beyond there
        **logo_plot_args,
):
    ''' Make tcr clumping logos for each chain individually, based
    on previously calculated single-chain clumping results dataframe
    '''
    num_clones = adata.shape[0]

    for chain in 'AB':
        fake_clusters_gex = np.zeros((num_clones,)).astype(int)
        fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
        clumping_pvals = np.full( (num_clones,), num_clones).astype(float)


        for l in results_df.itertuples():
            if l.chain == chain:
                clumping_pvals[ l.clone_index] = min(l.pvalue_adj,
                                                     clumping_pvals[l.clone_index])
                fake_clusters_tcr[l.clone_index] = l.clumping_group

        pngfile = f'{outfile_prefix}_{chain}_chain_clumping_logos.png'

        if 'rank_genes_uns_tag' not in logo_plot_args:
            logo_plot_args['rank_genes_uns_tag'] = f'rg_tcr_clumping_{chain}_biclusters'

        if 'conga_scores_name' not in logo_plot_args:
            logo_plot_args['conga_scores_name'] = f'TCR{chain} clumping'

        if 'show_real_clusters_gex' not in logo_plot_args:
            logo_plot_args['show_real_clusters_gex'] = True

        plotting.make_cluster_logo_plots_figure(
            adata, clumping_pvals, pvalue_threshold_for_logos,
            fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
            min_cluster_size_for_logos, pngfile,
            **logo_plot_args)





def assign_cd4_and_cd8_by_clusters(
        adata,
        key_added = 'cd4_or_cd8',
        clustering_resolution = 2.0,
        n_gex_pcs = 40,
        n_neighbors = 10,
        verbose = False,
):
    ''' adds new string column with name=key_added to adata.obs, values are
    'cd4' or 'cd8'

    Does not respect clone definitions (expanded clones may span both cd4 and cd8)

    assumes:
    * data have been preprocessed/scaled
    * data have not been reduced to a single cell per clone
    * pca may or may not have been called

    '''


    assert adata.uns['organism'] == 'human' # tmp hack

    # run pca if necessary
    if 'X_pca_gex' not in adata.obsm_keys():
        n_gex_pcs = min(adata.shape[0]-1, n_gex_pcs)
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_gex_pcs)
        adata.obsm['X_pca_gex'] = adata.obsm['X_pca']

    # run clustering with higher resolution
    adata.obsm['X_pca'] = adata.obsm['X_pca_gex']
    n_pcs = adata.obsm['X_pca'].shape[1]

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    cluster_key_added = 'leiden_gex_for_cd4_vs_cd8'

    sc.tl.leiden(adata, resolution=clustering_resolution, key_added=cluster_key_added)


    clusters_gex = np.array(adata.obs[cluster_key_added].astype(int))

    num_clones = adata.shape[0]
    all_gex = {}

    cd4_genes = ['CD4']
    cd8_genes = ['CD8A', 'CD8B']

    for gene in cd4_genes+cd8_genes:
        if gene not in adata.raw.var_names:
            print('WARNING assign_cd4_and_cd8_by_clusters: missing', gene)
            all_gex[gene] = np.zeros((num_clones,))
        else:
            index = adata.raw.var_names.get_loc(gene)
            all_gex[gene] = adata.raw.X[:,index].toarray()[:,0]

    cd8_gex = np.copy(all_gex[cd8_genes[0]])
    for g in cd8_genes[1:]:
        cd8_gex += all_gex[g]
    cd8_gex /= len(cd8_genes)

    cd4_gex = np.copy(all_gex[cd4_genes[0]])
    for g in cd4_genes[1:]:
        cd4_gex += all_gex[g]
    cd4_gex /= len(cd4_genes)

    num_clusters = np.max(clusters_gex)+1
    cd4_rescale = 1/0.6
    cd8_rescale = 1/2.0
    good_mask = np.full((num_clones,), False)
    dfl = []
    cd48s = np.full((num_clones,), 'UNK')
    for c in range(num_clusters):
        mask = clusters_gex==c
        cd4 = cd4_rescale * np.sum(cd4_gex[mask])/np.sum(mask)
        cd8 = cd8_rescale * np.sum(cd8_gex[mask])/np.sum(mask)
        tag = 'cd4' if cd4>cd8 else 'cd8'
        cd48s[mask] = tag
        if verbose:
            print(f'split_by_CD4_CD8: cluster {c:2d} {np.sum(mask):4d} {tag}',
                  f'cd4: {cd4:9.3f} cd8: {cd8:9.3f}')
        dfl.append(dict(
            cluster=c,
            size=np.sum(mask),
            cd4=cd4,
            cd8=cd8,
            ))
    adata.obs[key_added] = cd48s
    return pd.DataFrame(dfl)



def split_into_cd4_and_cd8_subsets(
        adata,
        allow_split_clones = True, # allow clones that span subsets
        max_iterations = 5,
        verbose = False,
        min_cells_for_iteration = 25,
):
    ''' returns adata_cd4, adata_cd8

    assumes:
    * data have been preprocessed/scaled
    * data have not been reduced to a single cell per clone
    * pca may or may not have been called

    takes the consensus for expanded clones, but ties contribute cells to both

    '''

    assign_cd4_and_cd8_by_clusters(adata, verbose=verbose)

    ad4 = adata[adata.obs.cd4_or_cd8 == 'cd4'].copy()
    ad8 = adata[adata.obs.cd4_or_cd8 == 'cd8'].copy()

    for r in range(max_iterations):
        print('start sizes:', r, ad4.shape[0], ad8.shape[0])
        if (ad4.shape[0] < min_cells_for_iteration or
            ad8.shape[0] < min_cells_for_iteration):
            break

        # look at clone overlap
        # tcrs_4 = preprocess.retrieve_tcrs_from_adata(ad4)
        # tcrs_8 = preprocess.retrieve_tcrs_from_adata(ad8)
        # counts_4 = Counter(tcrs_4)
        # counts_8 = Counter(tcrs_8)

        # for tcr,c4 in counts_4.items():
        #     c8 = counts_8[tcr]
        #     if c8:
        #         print(f'split_clone: {r} {c4:3d} {c8:3d} {tcr[0][2]} {tcr[1][2]}')

        # make new adatas
        dfl = []
        for ad in [ad4,ad8]:
            ad.uns['organism'] = 'human'
            assign_cd4_and_cd8_by_clusters(ad, verbose=verbose)

        ad4_new = ad4[ad4.obs.cd4_or_cd8 == 'cd4'].concatenate(
            ad8[ad8.obs.cd4_or_cd8 == 'cd4'], index_unique=None)
        ad8_new = ad8[ad8.obs.cd4_or_cd8 == 'cd8'].concatenate(
            ad4[ad4.obs.cd4_or_cd8 == 'cd8'], index_unique=None)

        # convergence?
        old_cd8_barcodes = set(ad8.obs.index)
        new_cd8_barcodes = set(ad8_new.obs.index)
        print(f'old_cd8_barcodes: {len(old_cd8_barcodes)} new_cd8_barcodes: '
              f'{len(new_cd8_barcodes)}')
        if (len(old_cd8_barcodes) == len(new_cd8_barcodes) ==
            len(old_cd8_barcodes&new_cd8_barcodes)):
            print('converged')
            break

        ad4 = ad4_new
        ad8 = ad8_new

    # now map from barcodes to cd4/cd8
    cd4_barcodes = set(ad4.obs.index)
    cd8_barcodes = set(ad8.obs.index)
    assert len(cd4_barcodes & cd8_barcodes) == 0

    # barcode_map is dict from barcode to integer (4 or 8)
    barcode_map = {x:4 for x in cd4_barcodes}
    barcode_map.update({x:8 for x in cd8_barcodes})


    # now assign clone by clone
    tcrs = preprocess.retrieve_tcrs_from_adata(adata) # has duplicates
    tcrs_sorted = sorted(set(tcrs)) # no duplicates
    tcr_to_clone_id = {x:i for i,x in enumerate(tcrs_sorted)}
    clone_ids = np.array([tcr_to_clone_id[x] for x in tcrs])

    all_barcodes = np.array(adata.obs.index)

    new_barcode_map = {}
    for clone_id, tcr in enumerate(tcrs_sorted):
        barcodes = list(all_barcodes[clone_ids==clone_id])
        counts = Counter([barcode_map[x] for x in barcodes]).most_common()
        if len(counts)>1 and verbose:
            print('split clone:', counts)

        if len(counts)>1 and counts[0][1] == counts[1][1] and allow_split_clones:
            # clone is split between cd4 and cd8, split by cell assignments
            for barcode in barcodes:
                new_barcode_map[barcode] = barcode_map[barcode]
        else:
            for barcode in barcodes:
                new_barcode_map[barcode] = counts[0][0]

    cd4_mask = np.array([new_barcode_map[x] == 4 for x in all_barcodes])
    cd8_mask = np.array([new_barcode_map[x] == 8 for x in all_barcodes])

    assert all(cd4_mask | cd8_mask)
    assert not any(cd4_mask & cd8_mask)

    adata_cd4 = adata[cd4_mask].copy()
    adata_cd8 = adata[cd8_mask].copy()

    return adata_cd4, adata_cd8



def run_umap_and_clustering_from_indices_distances(
        adata,
        knn_indices,
        knn_distances,
        umap_key_added,
        cluster_key_added,
        clustering_resolution = None,
        clustering_method=None,
        n_components_umap = 2,
):
    ''' code is borrowed from conga.preprocess.calc_tcrdist_nbrs_umap_clusters_cpp

    return TRUE on success, FALSE on failure (umap graph is too disconnected)

    '''
    assert knn_indices.shape == knn_distances.shape
    num_nbrs = knn_indices.shape[1]

    try: # HACK: the naming of this function changes across scanpy versions
        distances, connectivities= sc.neighbors.compute_connectivities_umap(
            knn_indices, knn_distances, adata.shape[0], num_nbrs)
    except:
        print('try new name for compute_connectivities_umap')
        distances, connectivities= sc.neighbors._compute_connectivities_umap(
            knn_indices, knn_distances, adata.shape[0], num_nbrs)

    if issparse(connectivities): # I think this is always true
        from scipy.sparse.csgraph import connected_components
        connected_components = connected_components(connectivities)
        number_connected_components = connected_components[0]
        print('number_connected_components:', number_connected_components)
        if number_connected_components > 2*n_components_umap:
            print('run_umap_and_clustering_from_indices_distances:',
                  'too many connected components in the',
                  'neighbor graph:', number_connected_components)
            return False # signal failure

    ################
    # stash the stuff in adata, stolen from scanpy/neighbors/__init__.py
    #
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params']={'n_neighbors': num_nbrs, 'method': 'umap'}
    adata.uns['neighbors']['params']['metric'] = 'tcrdist'# fake metric
    # if metric_kwds:
    #     adata.uns['neighbors']['params']['metric_kwds'] = metric_kwds
    # if use_rep is not None:
    #     adata.uns['neighbors']['params']['use_rep'] = use_rep
    # if n_pcs is not None:
    #     adata.uns['neighbors']['params']['n_pcs'] = n_pcs
    adata.uns['neighbors']['distances'] = distances
    adata.uns['neighbors']['connectivities'] = connectivities

    # as far as I can tell, these are only used if there are too many connected
    # components in the nbr graph... see the infinite while loop up above.
    print('temporarily putting random pca vectors into adata...')
    fake_pca = np.random.randn(adata.shape[0], 10)
    adata.obsm['X_pca'] = fake_pca

    print('running umap', adata.shape)
    sc.tl.umap(adata, n_components=n_components_umap)
    print('DONE running umap')
    adata.obsm[umap_key_added] = adata.obsm['X_umap']


    resolution = 1.0 if clustering_resolution is None else clustering_resolution
    if clustering_method=='louvain':
        sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
        print('ran louvain clustering:', resolution, cluster_key_added)
    elif clustering_method=='leiden':
        sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key_added)
        print('ran leiden clustering:', resolution, cluster_key_added)
    else: # try both (hacky)
        try:
            sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key_added)
            print('ran leiden clustering:', resolution, cluster_key_added)
        except ImportError: # hacky
            sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
            print('ran louvain clustering:', resolution, cluster_key_added)

    adata.obs[cluster_key_added] = np.copy(adata.obs[cluster_key_added]).astype(int)
    print('DONE running louvain', cluster_key_added)

    del adata.obsm['X_pca'] # delete the fake pcas
    del adata.obsm['X_umap'] # delete the extra umap copy

    return True


def create_conga_hits_adata(
        adata,
        outfile_prefix, # for plots
        small_nbr_frac = 0.1,
        big_nbr_frac = 0.25,
        num_nbrs = 10, # for UMAP/clustering
        verbose = True,
):
    ''' Try UMAP/clustering just the conga hits, using a measure that combines GEX and
    TCR

    returns a new adata which just contains conga hits, and has new fields

    adata.obs.clusters_combo
    adata.obsm.X_combo_2d

    '''

    # subset to the conga hits
    mask = adata.obs.conga_scores < 1.0
    adata = adata[mask].copy()

    all_nbrs = preprocess.calc_nbrs(
        adata, [small_nbr_frac, big_nbr_frac], sort_nbrs=True)

    nbrs_gex, nbrs_tcr= all_nbrs[big_nbr_frac]

    N = adata.shape[0]
    overlaps = []

    knn_indices = np.zeros((N,num_nbrs), dtype=np.int32)
    knn_distances = np.zeros((N,num_nbrs), dtype=float)

    big_num_nbrs = nbrs_gex.shape[1]
    scores = np.zeros((big_num_nbrs,))

    for ii in range(N):
        if ii and ii%1000==0 and verbose:
            print(ii)
        gnbrs = nbrs_gex[ii]
        tnbrs_index = {x:i for i,x in enumerate(nbrs_tcr[ii])}
        for j,jj in enumerate(gnbrs):
            k = tnbrs_index.get(jj,0.5*N)
            scores[j] = 0.5*(j+k)/N
        inds = np.argpartition(scores, num_nbrs-1)[:num_nbrs]
        knn_distances[ii,:] = scores[inds]
        knn_indices[ii,:] = gnbrs[inds]


    run_umap_and_clustering_from_indices_distances(
        adata,
        knn_indices,
        knn_distances,
        'X_combo_2d',
        'clusters_combo',
        #clustering_resolution=2.0,
    )

    # show new UMAP space colored in various ways
    pngfile = outfile_prefix+'_combo_umaps.png'
    xy = adata.obsm['X_combo_2d']

    genes = ('conga_scores clone_sizes KLRB1 GZMK GZMH CCL5 SELL '
             'CCR7 IKZF2 ZNF683 KIR TCF7'.split())

    kir_genes = ['KIR2DL1', 'KIR2DL3', 'KIR2DL4', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3']

    nrows, ncols = 3,4
    plt.figure(figsize=(ncols*4,nrows*4))
    for ig,gene in enumerate(genes):
        plt.subplot(nrows, ncols, ig+1)
        if gene == 'conga_scores':
            colors = np.sqrt(-1*np.log10(adata.obs.conga_scores))
        elif gene == 'clone_sizes':
            colors = np.log10(adata.obs.clone_sizes)
        elif gene == 'KIR':
            colors = np.sum(np.vstack(
                [adata.raw.X[:,adata.raw.var_names.get_loc(g)].toarray()[:,0]
                 for g in kir_genes if g in adata.raw.var_names]),axis=0)
        else:
            colors = adata.raw.X[:,adata.raw.var_names.get_loc(gene)].toarray()[:,0]
        reorder = np.argsort(colors)
        plt.scatter(xy[reorder,0], xy[reorder,1], c=colors[reorder], s=5)
        plt.title(gene)
        plt.xticks([],[])
        plt.yticks([],[])
    plt.tight_layout()
    plt.savefig(pngfile)
    print('made:', pngfile)

    #borrowed from  plotting.make_tcr_clumping_plots
    fake_clusters_gex = np.zeros((N,)).astype(int)
    fake_clusters_tcr = np.array(adata.obs.clusters_combo)

    pngfile = outfile_prefix+'_combo_clusters_old_2d.png'

    nbrs_gex, nbrs_tcr = all_nbrs[small_nbr_frac]

    min_cluster_size_for_logos=3

    plotting.make_cluster_logo_plots_figure(
        adata, adata.obs.conga_scores, 10.0, # 10.0 is congathreshold (anything above 1)
        fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile, show_real_clusters_gex=True,
    )
    print('made:', pngfile)

    X_gex_2d_save = np.array(adata.obsm['X_gex_2d']).copy()
    X_tcr_2d_save = np.array(adata.obsm['X_tcr_2d']).copy()
    adata.obsm['X_gex_2d'] = adata.obsm['X_combo_2d']
    adata.obsm['X_tcr_2d'] = adata.obsm['X_combo_2d']

    pngfile = outfile_prefix+'_combo_clusters_new_2d.png'
    plotting.make_cluster_logo_plots_figure(
        adata, adata.obs.conga_scores, 10.0, # 10.0 is congathreshold (anything above 1)
        fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile, show_real_clusters_gex=True,
    )
    print('made:', pngfile)

    adata.obsm['X_gex_2d'] = X_gex_2d_save
    adata.obsm['X_tcr_2d'] = X_tcr_2d_save

    return adata


## TEMPORARY CODE GRAVEYARD

# def _include_redundant_feature(count):
#     ''' count = 1 for first feature exceeding threshold, then 2, etc

#     idea:

#     '''
#     assert count>0
#     count += 2 # so 3 is the first (skip), then 4 (keep), then 5-7 (skip), etc
#     return abs(2**int(np.log2(count)) - count)<.1

    # if max_redundant_features is not None:
    #     feature_nbrs = {}

    #     # filter features by correlation
    #     # will subset: features, feature_types, feature_labels, A, nrows
    #     if verbose: print('computing feature correlations:', A.shape)
    #     C = 1-distance.squareform(distance.pdist(A, metric='correlation'),
    #                               force='tomatrix')
    #     if verbose: print('DONE computing feature correlations:', A.shape)
    #     feature_nbr_counts = [0]*len(features)
    #     feature_mask = np.full(len(features), True)
    #     for ii,f1 in enumerate(features):
    #         if verbose and ii%10==0: print('redundancy checking...', ii)
    #         # am I too close to a previous feature?
    #         for jj in range(ii-1):
    #             if ( feature_mask[jj] and
    #                  (feature_types is None or
    #                   feature_types[ii] == feature_types[jj])):
    #                 if C[ii,jj] > redundancy_threshold:
    #                     feature_nbr_counts[jj] += 1
    #                     if feature_nbr_counts[jj] > max_redundant_features:
    #                         count= feature_nbr_counts[jj]-max_redundant_features
    #                         if not _include_redundant_feature(count):
    #                             print('skip:', ii, jj, count)
    #                             feature_mask[ii] = False
    #                             feature_nbrs.setdefault(
    #                                 features[jj],[]).append(f1)
    #                         else:
    #                             print('keep:', ii, jj, count)
    #                         break
    # print('filling the score array for', len(features), 'features')
    # A = np.zeros((len(features), adata.shape[0]))
    # for ii,feature in enumerate(features):
    #     if verbose and ii%10==0:
    #         print('filling:', ii, len(features), feature)
    #     scores = get_raw_feature_scores(
    #         feature, adata,
    #         feature_types if feature_types is None else feature_types[ii])

    #     mn, std = np.mean(scores), np.std(scores)
    #     scores = (scores-mn)
    #     if std!=0:
    #         scores /= std

    #     if compute_nbr_averages:
    #         num_neighbors = nbrs.shape[1]
    #         scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

    #     if feature_types is not None:
    #         if feature_types[ii] == dist_tag:
    #             # lighten these up a bit since they will be nbr correlated
    #             scores *= rescale_factor_for_self_features

    #     A[ii,:] = scores
