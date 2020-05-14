import numpy as np
from scipy import stats
from scipy.stats import hypergeom, mannwhitneyu
from scipy.sparse import issparse
from collections import Counter
import scanpy as sc
from . import preprocess as pp
from . import tcr_scores
from . import util
import sys
import pandas as pd
from sys import exit

def compute_cluster_interactions( aclusters_in, bclusters_in, barcodes_in, barcode2tcr, outlog, max_pval = 1.0 ):

    ''' Compute the cluster-cluster intxn (positive) with the lowest hypergeometric pval, and report.
    Then eliminate one of the two interacting clusters: first try taking the one with the fewest cells
    outside the pairwise interaction... Then iterate.

    pval includes a rescaling factor = num_a_clusters * num_b_clusters, max_pval applies to this rescaled pval
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

    #all_mait_fractions = {} # going to be a dictionary from clusters and cluster pairs to mait fractions
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

        # if not all_mait_fractions:
        #     for a in range(num_a):
        #         a_mask = (aclusters==a)
        #         all_mait_fractions[ (a,num_b) ] \
        #             = sum( is_mait( barcode2tcr[x][0] ) for x in barcodes[a_mask] ) / float(max(1,np.sum(a_mask)))
        #         for b in range(num_b):
        #             ab_mask = a_mask & ( bclusters==b )
        #             all_mait_fractions[ (a,b) ] \
        #                 = sum( is_mait( barcode2tcr[x][0] ) for x in barcodes[ab_mask] ) / float(max(1,np.sum(ab_mask)))
        #     for b in range(num_b):
        #         b_mask = (bclusters==b)
        #         all_mait_fractions[ (num_a,b) ] \
        #             = sum( is_mait( barcode2tcr[x][0] ) for x in barcodes[b_mask] ) / float(max(1,np.sum(b_mask)))


        # most significant interaction:
        a,b = np.unravel_index( np.argmin( pvals ), pvals.shape )
        overlap = abcounts[a,b]
        if not overlap:
            break
        overlap_barcodes = barcodes[ ( aclusters==a ) & ( bclusters==b ) ]
        # mait_frac = sum( is_mait( barcode2tcr[x][0] ) for x in overlap_barcodes ) / float(overlap)
        # a_mait_frac = sum( is_mait( barcode2tcr[x][0] ) for x in barcodes[ aclusters==a ] ) / float(acounts[a])
        # b_mait_frac = sum( is_mait( barcode2tcr[x][0] ) for x in barcodes[ bclusters==b ] ) / float(bcounts[b])
        expected = acounts[a] * bcounts[b] / float(total)
        pval = pval_rescale * pvals[a,b]

        if pval > max_pval:
            break

        outlog.write('clusclus2_intxn: {:2d} {:2d} {:8.1e} {:3d} {:6.1f} {:4d} {:4d} {:2d} {:2d} {:5d} ex {}\n'\
                     .format( a,b,pval,overlap,expected,acounts[a],bcounts[b],
                              len(set(aclusters)),len(set(bclusters)),len(barcodes),
                              len(overlap_barcodes)-overlap ))
        # print('clusclus2_intxn: {:2d} {:2d} {:8.1e} {:3d} {:6.1f} {:4d} {:4d} {:2d} {:2d} {:5d} {:5.3f} {:5.3f} {:5.3f} ex {}'\
        #       .format( a,b,pval,overlap,expected,acounts[a],bcounts[b],
        #                len(set(aclusters)),len(set(bclusters)),len(barcodes),
        #                all_mait_fractions[(a,b)], all_mait_fractions[(a,num_b)], all_mait_fractions[(num_a,b)],
        #                len(overlap_barcodes)-overlap ))

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


def find_neighbor_neighbor_interactions(
        adata,
        nbrs_gex,
        nbrs_tcr,
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True
):
    ''' Returns a pandas dataframe with results
    of all clones with nbr-nbr overlaps having pvalues <= pval_threshold

    '''
    pp.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_mait']

    num_clones = len(nbrs_gex)

    pval_rescale = num_clones if scale_pvals_by_num_clones else 1.0

    results = []
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

        if nbr_pval > pval_threshold:
            continue ## NOTE

        double_nbrs = [ x for x in ii_nbrs_tcr if x in ii_nbrs_gex_set ]
        assert overlap == len(double_nbrs)

        overlap_corrected = overlap
        if correct_overlaps_for_groups:
            overlap_corrected = min( len(set(agroups[double_nbrs])), len(set(bgroups[double_nbrs])) )

            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                nbr_pval = pval_rescale * hypergeom.sf( overlap_corrected-1, actual_num_clones, num_neighbors_gex-delta,
                                                        num_neighbors_tcr-delta)
                if nbr_pval > pval_threshold:
                    continue ## NOTE

        results.append( dict( conga_score=nbr_pval,
                              num_neighbors_gex=num_neighbors_gex,
                              num_neighbors_tcr=num_neighbors_tcr,
                              overlap=overlap,
                              overlap_corrected=overlap_corrected,
                              mait_fraction=np.sum(is_mait[double_nbrs])/overlap,
                              clone_index=ii ))#double_nbrs ] )

    return pd.DataFrame(results)


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
    '''
    pp.add_mait_info_to_adata_obs(adata) # for annotation of overlaps
    is_mait = adata.obs['is_mait']

    num_clones = len(nbrs)

    pval_rescale = num_clones if scale_pvals_by_num_clones else 1.0

    results = []
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

        if nbr_pval > pval_threshold:
            continue ## NOTE

        same_cluster_nbrs = ii_nbrs[ clusters[ii_nbrs] == ii_cluster ]
        assert len(same_cluster_nbrs)== overlap

        overlap_corrected = overlap
        if correct_overlaps_for_groups:
            overlap_corrected = min( len(set(agroups[same_cluster_nbrs])), len(set(bgroups[same_cluster_nbrs])) )
            if overlap_corrected < overlap:
                delta = overlap-overlap_corrected
                new_nbr_pval = pval_rescale * hypergeom.sf( overlap_corrected-1, actual_num_clones, num_neighbors-delta,
                                                            ii_cluster_clustersize-delta )


        results.append( dict( conga_score=nbr_pval,
                              num_neighbors=num_neighbors,
                              cluster_size=ii_cluster_clustersize,
                              overlap=overlap,
                              overlap_corrected=overlap_corrected,
                              mait_fraction=np.sum(is_mait[same_cluster_nbrs])/overlap,
                              clone_index=ii ))#double_nbrs ] )

    return pd.DataFrame(results)

def run_rank_genes_on_good_cluster_pairs(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        rank_method='wilcoxon',
        rg_tag = 'test',
        neg_tag='none',
        min_count=5,
        key_added = 'rank_genes_good_cluster_pairs'
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
        print('run_rank_genes_on_good_cluster_pairs: no good cluster pairs')
        return

    adata.obs[rg_tag] = vals

    print('run rank_genes_groups', Counter(vals).most_common() )

    sc.tl.rank_genes_groups( adata, groupby=rg_tag, method=rank_method, groups=pos_tags, reference='rest',
                             key_added = key_added )

    # retvals = []

    # for igene,gene in enumerate( adata.uns['rank_genes_groups']['names'][pos_tag] ):
    #     log2fold = adata.uns['rank_genes_groups']['logfoldchanges'][pos_tag][igene]
    #     pval_adj = adata.uns['rank_genes_groups']['pvals_adj'][pos_tag][igene]
    #     score = adata.uns['rank_genes_groups']['scores'][pos_tag][igene]
    #     if pval_adj <= max_pval_for_output:
    #         print('cellslog2fold: {:3d} {:.2f} pval_adj: {:9.1e} score: {:.1f} {} {} {}'\
    #               .format( igene, log2fold, pval_adj, score, gene, numpos, pos_tag ) )
    #     retvals.append( [ igene, gene, log2fold, pval_adj, score ] )

    # sys.stdout.flush()
    # return retvals


def run_rank_genes_on_cells(
        adata,
        mask,
        key_added='rank_genes_on_cells',
        pos_tag='pos',
        #rank_method='t-test',
        rank_method='wilcoxon',
        rg_tag = 'test',
        neg_tag='neg',
        min_pos=3,
        max_pval_for_output=10,
        n_genes=100
):
    ''' Returns list of the top n_genes [rank(0-indexed), gene, log2fold, pval_adj, score ]
    '''
    assert mask.shape[0] == adata.shape[0]
    if np.sum(mask) < min_pos:
        return []

    if not pp.check_if_raw_matrix_is_logged(adata):
        print('ERROR need to log the raw matrix before calling run_rank_genes_on_cells')
        exit()

    vals = [ pos_tag if a else neg_tag for a in mask ]

    adata.obs[rg_tag] = vals

    sc.tl.rank_genes_groups(adata, groupby=rg_tag, method=rank_method, groups=[pos_tag], reference=neg_tag,
                            key_added=key_added, n_genes=n_genes)

    retvals = []

    for igene,gene in enumerate( adata.uns[key_added]['names'][pos_tag] ):
        log2fold = adata.uns[key_added]['logfoldchanges'][pos_tag][igene]
        pval_adj = adata.uns[key_added]['pvals_adj'][pos_tag][igene]
        score = adata.uns[key_added]['scores'][pos_tag][igene]
        retvals.append( [ igene, gene, log2fold, pval_adj, score ] )

    return retvals

def calc_good_cluster_tcr_features(
        adata,
        good_mask,
        clusters_gex,
        clusters_tcr,
        tcr_score_names,
        min_count=5,
):
    num_clones = adata.shape[0]

    clp_counts = Counter( (x,y) for x,y,z in zip( clusters_gex, clusters_tcr, good_mask ) if z )
    print( clp_counts.most_common())

    good_clps = [ x for x,y in clp_counts.items() if y>=min_count]

    fake_nbrs_gex = []
    seen = set()
    for cl_gex, cl_tcr, m in zip( clusters_gex, clusters_tcr, good_mask ):
        clp=(cl_gex, cl_tcr)
        if m and clp in good_clps and clp not in seen:
            seen.add(clp)
            fake_nbrs_gex.append( np.nonzero( (clusters_gex==cl_gex) & (clusters_tcr==cl_tcr) )[0] )
        else:
            fake_nbrs_gex.append([])

    pval_threshold = 1.
    results_df = gex_nbrhood_rank_tcr_scores( adata, fake_nbrs_gex, tcr_score_names, pval_threshold,
                                              prefix_tag = 'good_clp' )

    all_tcr_features = {}
    for row in results_df.itertuples():
        clp = (row.gex_cluster, row.tcr_cluster)
        assert clp in good_clps
        all_tcr_features.setdefault(clp,[]).append( ( row.mwu_pvalue_adj, row.ttest_stat, row.score_name))

    for clp in all_tcr_features:
        # plotting code expects name then stat then pvalue
        # but we want them sorted by significance
        all_tcr_features[clp] = [ (x[2], x[1], x[0]) for x in sorted(all_tcr_features[clp]) ]

    return all_tcr_features


# def get_nbrhood_infostrings( adata, nbrs ):
#     ''' Returns a list of infostrings of len= adata.shape[0]
#     clusters and consensus clusters
#     mait fraction
#     '''

#     tcrs = pp.retrieve_tcrs_from_adata(adata)

#     clusters_gex = adata.obs['clusters_gex']
#     clusters_tcr = adata.obs['clusters_tcr']

#     infostrings = []
#     for ii in range(adata.shape[0]):
#         tcr_string = '{} {}'.format(' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3]) )
#         is_mait = tcr_scores.is_human_mait_alpha_chain(tcrs[ii][0])
#         mait_count = sum( tcr_scores.is_human_mait_alpha_chain(tcrs[x][0]) for x in nbrs[ii]) + is_mait
#         mait_frac = mait_count/(1+len(nbrs[ii])) ## including ii too
#         clp = ( clusters_gex[ii], clusters_tcr[ii])
#         clp_counts = Counter( zip( (clusters_gex[x] for x in nbrs[ii]), (clusters_tcr[x] for x in nbrs[ii]) ) )
#         clp_counts[clp] += 1
#         top_clp = clp_counts.most_common(1)[0][0]
#         info = 'clp {:2d} {:2d} top {:2d} {:2d} {} m: {:d} {:.3f} {}'\
#                .format( clp[0], clp[1], top_clp[0], top_clp[1], tcr_string, is_mait, mait_frac, ii)
#         infostrings.append(info)
#     return infostrings



# def tcr_nbrhood_rank_genes( adata, nbrs_tcr, pval_threshold, rank_method=None):
#     assert False # not using this anymore?
#     if rank_method is None:
#         rank_method = 'wilcoxon'

#     num_clones = adata.shape[0]
#     tcrs = pp.retrieve_tcrs_from_adata(adata)

#     for ii in range(num_clones):
#         nbrhood_mask = np.full( (num_clones,), False)
#         nbrhood_mask[ nbrs_tcr[ii] ] = True
#         nbrhood_mask[ ii ] = True

#         print('run_rank_genes_on_cells:', ii, num_clones)
#         results = run_rank_genes_on_cells(adata, nbrhood_mask, rank_method=rank_method )
#         for igene, gene, log2fold, pval_adj, score in results:
#             if pval_adj < pval_threshold and gene.lower()[:4] not in ['trav','trbv']:
#                 print('tcr_nbrhood_rank_genes: {:4d} {:2d} {:9.2e} {:7.3f} {} {} {}'\
#                       .format( ii, igene, pval_adj, log2fold, gene, ' '.join(tcrs[ii][0][:3]),
#                                ' '.join(tcrs[ii][1][:3]) ) )
#         sys.stdout.flush()

# def get_split_mean_var( X, mask, mean, mean_sq ):
#     # use: var = (mean_sq - mean**2)
#     #
#     N = mask.shape[0]
#     assert X.shape[0] == N
#     fg_wt = np.sum(mask) / N
#     bg_wt = 1. - fg_wt
#     fg_X = X[mask]
#     fg_mean = fg_X.mean(axis=0)
#     fg_mean_sq = np.multiply(fg_X, fg_X).mean(axis=0)
#     bg_mean = (mean - fg_wt*fg_mean)/bg_wt
#     bg_mean_sq = (mean_sq - fg_wt*fg_mean_sq)/bg_wt
#     fg_var = (fg_mean_sq - fg_mean**2)
#     bg_var = (bg_mean_sq - bg_mean**2)
#     return fg_mean, fg_var, bg_mean, bg_var
def get_split_mean_var( X, X_sq, mask, mean, mean_sq ):
    # use: var = (mean_sq - mean**2)
    #
    N = mask.shape[0]
    assert X.shape[0] == N
    wt_fg = np.sum(mask) / N
    wt_bg = 1. - wt_fg
    mean_fg = X[mask].mean(axis=0)
    mean_sq_fg = X_sq[mask].mean(axis=0)
    if issparse(X):
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
    clusters_gex = adata.obs['clusters_gex']
    clusters_tcr = adata.obs['clusters_tcr']

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold*10

    #nbrhood_infos = get_nbrhood_infostrings(adata, nbrs_gex)
    print('making tcr score table:', tcr_score_names)
    score_table = tcr_scores.make_tcr_score_table(adata, tcr_score_names)
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
            mean1=mean_fg, std1=np.sqrt(var_fg), nobs1=num_fg,
            mean2=mean_bg, std2=np.sqrt(var_bg), nobs2=num_clones-num_fg, # scanpy ttest-over-estim-var uses nobs1 here
            equal_var=False  # Welch's
        )

        scores[np.isnan(scores)] = 0  # I think it's only nan when means are the same and vars are 0
        pvals[np.isnan(pvals)] = 1  # This also has to happen for Benjamini Hochberg

        # crude bonferroni
        pvals *= pval_rescale

        nbrhood_clusters_gex, nbrhood_clusters_tcr = None,None # lazy

        for ind in np.argsort(pvals):
            pval = pvals[ind]

            if pval>ttest_pval_threshold_for_mwu_calc:
                continue

            _,mwu_pval = mannwhitneyu( score_table[:,ind][nbrhood_mask], score_table[:,ind][~nbrhood_mask] )
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

                if verbose:
                    print('gex_{}_score: {:9.2e} {:9.2e} {:7.2f} clp {:2d} {:2d} {:7.3f} {:7.3f} {:15s} {} {} mf: {:.3f} {}'\
                          .format( prefix_tag, pval, mwu_pval_adj, score, gex_cluster, tcr_cluster, mean_fg[ind],
                                   mean_bg[ind], score_name,
                                   ' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3]), mait_fraction, ii ))

                results.append( dict(ttest_pvalue_adj=pval,
                                     ttest_stat=score,
                                     mwu_pvalue_adj=mwu_pval_adj,
                                     gex_cluster=gex_cluster,
                                     tcr_cluster=tcr_cluster,
                                     score_name=score_name,
                                     mait_fraction=mait_fraction,
                                     clone_index=ii) )

        sys.stdout.flush()
        # if corr_method == 'benjamini-hochberg':
        #     _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        # elif corr_method == 'bonferroni':
        #     pvals_adj = np.minimum(pvals * n_genes, 1.0)
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
    clusters_gex = adata.obs['clusters_gex']
    clusters_tcr = adata.obs['clusters_tcr']
    tcrs = pp.retrieve_tcrs_from_adata(adata)
    pp.add_mait_info_to_adata_obs(adata)
    is_mait = adata.obs['is_mait']
    ## done unpacking

    if clone_display_names is None:
        clone_display_names = [ '{} {}'.format(' '.join(x[0][:3]), ' '.join(x[1][:3])) for x in tcrs ]

    if ttest_pval_threshold_for_mwu_calc is None:
        ttest_pval_threshold_for_mwu_calc = pval_threshold * 10

    #corr_method = 'benjamini-hochberg'

    #from scipy import stats
    from statsmodels.stats.multitest import multipletests
    assert pp.check_if_raw_matrix_is_logged(adata)

    rankby_abs = False

    num_clones = adata.shape[0]
    #tcrs = pp.retrieve_tcrs_from_adata(adata)
    #nbrhood_infos = get_nbrhood_infostrings(adata, nbrs_tcr)

    genes = adata.raw.var_names
    X = adata.raw.X
    assert issparse(X)
    assert X.shape[1] == len(genes)

    print('to csc')
    X_csc = adata.raw.X.tocsc()

    print('square X')
    X2 = X.multiply(X)
    print('done squaring X')

    mean = X.mean(axis=0)
    mean_sq = X2.mean(axis=0)
    print('done taking big means')

    reference_indices = np.arange(X.shape[1], dtype=int)
    # need this?
    mean = mean.A1
    mean_sq = mean_sq.A1

    num_nonempty_nbrhoods = sum(1 for x in nbrs_tcr if len(x)>0)

    # len(genes) is probably too hard since lots of the genes are all zeros
    min_nonzero_cells = 3
    gene_nonzero_counts = Counter( X.nonzero()[1] )
    bad_gene_mask = np.array([ gene_nonzero_counts[x] < min_nonzero_cells for x in range(len(genes)) ])
    n_genes_eff = np.sum(~bad_gene_mask)
    pval_rescale = num_nonempty_nbrhoods * n_genes_eff

    results = []

    for ii in range(num_clones):
        if len(nbrs_tcr[ii])==0:
            continue
        nbrhood_mask = np.full( (num_clones,), False)
        nbrhood_mask[ nbrs_tcr[ii] ] = True
        nbrhood_mask[ ii ] = True

        mean_fg, var_fg, mean_bg, var_bg = get_split_mean_var(X, X2, nbrhood_mask, mean, mean_sq)
        num_fg = np.sum(nbrhood_mask)
        if num_fg < min_num_fg:
            continue
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean_fg, std1=np.sqrt(var_fg), nobs1=num_fg,
            mean2=mean_bg, std2=np.sqrt(var_bg), nobs2=num_clones-num_fg,
            equal_var=False  # Welch's
        )

        # scanpy code:
        scores[np.isnan(scores)] = 0.
        pvals [np.isnan(pvals)] = 1.
        pvals [bad_gene_mask] = 1.
        logfoldchanges = np.log2((np.expm1(mean_fg) + 1e-9) / (np.expm1(mean_bg) + 1e-9))

        pvals_adj = pvals * pval_rescale
        # if corr_method == 'benjamini-hochberg':
        #     _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
        # elif corr_method == 'bonferroni':
        #     pvals_adj = np.minimum(pvals * len(genes), 1.0)
        # else:
        #     print('unrecognized corr_method:', corr_method)
        #     exit()

        scores_sort = np.abs(scores) if rankby_abs else scores
        partition = np.argpartition(scores_sort, -top_n)[-top_n:]
        partial_indices = np.argsort(scores_sort[partition])[::-1]
        global_indices = reference_indices[partition][partial_indices]

        nbrhood_clusters_gex, nbrhood_clusters_tcr = None,None

        for igene, ind in enumerate(global_indices):
            gene = genes[ind]
            if gene.lower()[:4] in ['trav','trbv']:
                continue
            pval_adj = pvals_adj[ind]
            log2fold= logfoldchanges[ind]

            if pval_adj > ttest_pval_threshold_for_mwu_calc:
                continue

            col = X_csc[:,ind][nbrhood_mask]
            noncol = X_csc[:,ind][~nbrhood_mask]
            _,mwu_pval = mannwhitneyu( col.todense(), noncol.todense() )
            mwu_pval_adj = mwu_pval * pval_rescale

            if min(mwu_pval_adj, pval_adj) < pval_threshold:
                if nbrhood_clusters_gex is None: # lazy
                    nbrhood_clusters_gex = clusters_gex[nbrhood_mask]
                    nbrhood_clusters_tcr = clusters_tcr[nbrhood_mask]
                    nbrhood_is_mait = is_mait[nbrhood_mask]


                # better annotation of the enriched tcrs...
                num_top = num_fg//4
                if len(col.data)>num_top: # more than a quarter non-zero
                    top_indices = col.indices[ np.argpartition(col.data, -num_top)[-num_top:] ]
                    #bot_indices = col.indices[ np.argpartition(col.data, num_top-1)[:num_top] ]
                    #assert np.mean(col[top_indices]) > np.mean(col[bot_indices])
                else:
                    top_indices = col.indices
                #top_indices = np.nonzero(nbrhood_mask)[0][col_top_inds]

                gex_cluster = Counter( nbrhood_clusters_gex[ top_indices ]).most_common(1)[0][0]
                tcr_cluster = Counter( nbrhood_clusters_tcr[ top_indices ]).most_common(1)[0][0]
                mait_fraction = np.sum(nbrhood_is_mait[ top_indices ] )/len(top_indices)


                if verbose:
                    print('tcr_{}_gene: {:9.2e} {:9.2e} {:7.3f} clp {:2d} {:2d} {:8s} {:.4f} {:.4f} {:4d} {} mf: {:.3f} {} {}'\
                          .format( prefix_tag, pval_adj, mwu_pval_adj, log2fold, gex_cluster, tcr_cluster, gene,
                                   mean_fg[ind], mean_bg[ind], num_fg, clone_display_names[ii], mait_fraction,
                                   ii, igene ))
                results.append( { 'ttest_pvalue_adj': pval_adj,
                                  'mwu_pvalue_adj': mwu_pval_adj,
                                  'log2enr': log2fold,
                                  'gex_cluster': gex_cluster,
                                  'tcr_cluster': tcr_cluster,
                                  'gene': gene,
                                  'mean_fg': mean_fg[ind],
                                  'mean_bg': mean_bg[ind],
                                  'num_fg': num_fg,
                                  'clone_index': ii,
                                  'mait_fraction': mait_fraction } )

        sys.stdout.flush()

    return pd.DataFrame(results)
