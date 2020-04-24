import numpy as np
from scipy.stats import hypergeom
from collections import Counter
import scanpy as sc

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
        nbrs_gex,
        nbrs_tcr,
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True
):
    ''' Returns a list: [ [ pvalue, index, overlapping_indices ], ...]
    of all clones with nbr-nbr overlaps having pvalues <= pval_threshold

    '''
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

        if correct_overlaps_for_groups:
            new_overlap = min( len(set(agroups[double_nbrs])), len(set(bgroups[double_nbrs])) )
            if new_overlap < overlap:
                delta = overlap-new_overlap
                nbr_pval = pval_rescale * hypergeom.sf( new_overlap-1, actual_num_clones, num_neighbors_gex-delta,
                                                        num_neighbors_tcr-delta)
                if nbr_pval > pval_threshold:
                    continue ## NOTE


        results.append( [ nbr_pval, ii, double_nbrs ] )
    return results


def find_neighbor_cluster_interactions(
        nbrs,
        clusters, # for example, or could be swapped: tcr/gex
        agroups,
        bgroups,
        pval_threshold,
        counts_correction=0, # if there are a lot of maits, for example; this is the number of them
        correct_overlaps_for_groups=True,
        scale_pvals_by_num_clones=True
):
    ''' Returns a list: [ [ pvalue, index, overlapping_indices ], ...]
    of all clones with nbr-cluster overlaps having pvalues <= pval_threshold
    '''
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

        if correct_overlaps_for_groups:
            new_overlap = min( len(set(agroups[same_cluster_nbrs])), len(set(bgroups[same_cluster_nbrs])) )
            if new_overlap < overlap:
                delta = overlap-new_overlap
                new_nbr_pval = pval_rescale * hypergeom.sf( new_overlap-1, actual_num_clones, num_neighbors-delta,
                                                            ii_cluster_clustersize-delta )


        results.append( [ nbr_pval, ii, same_cluster_nbrs ] )
    return results

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
