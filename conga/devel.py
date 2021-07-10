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
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind
#from scipy.sparse import issparse, csr_matrix
import scipy.sparse as sps
from statsmodels.stats.multitest import multipletests
from collections import Counter, OrderedDict
import scanpy as sc
from . import preprocess
from . import tcr_scoring
from . import util
from . import correlations
from . import plotting
from .tcrdist.all_genes import all_genes
import sys
import pandas as pd
from sys import exit
import time #debugging



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
                                       alternative='two-sided')

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
                                           alternative='two-sided')

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
        _,mwu_p = mannwhitneyu(yvals[mask], yvals[~mask])
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

