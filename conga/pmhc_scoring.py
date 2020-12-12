from . import preprocess as pp
from . import util
import numpy as np
from scipy.stats import hypergeom
from scipy.special import binom
import sys
from collections import Counter
import pandas as pd

def get_gene_ids_varname( adata ):
    ''' Figure out the correct "gene_ids" varname, e.g. gene_ids-0 or gene_ids-0-0
    '''
    for name in adata.var:
        if name.startswith('gene_ids'):
            return name
    print('unable to find feature_types varname')
    print(adata.var_names)
    return None

def get_tenx_agbt_pmhc_var_names( adata ):
    raw = adata.raw
    if raw is None:
        raw = adata
    ab_capture_mask = ( raw.var[ util.get_feature_types_varname( adata )] == 'Antibody Capture' )
    not_totalseq_mask = ~raw.var_names.str.contains('TotalSeq')
    pmhc_mask = ab_capture_mask & not_totalseq_mask

    pmhc_var_names = list( raw.var_names[ pmhc_mask ] )
    return pmhc_var_names

def get_pmhc_short_and_long_names_dicts( adata ):
    ''' Returns dictionaries mapping from the feature names to short and long names that include HLA info
    '''
    raw = adata.raw
    if raw is None:
        raw = adata

    ab_capture_mask = ( raw.var[ util.get_feature_types_varname( adata )] == 'Antibody Capture' )
    not_totalseq_mask = ~raw.var_names.str.contains('TotalSeq')
    pmhc_mask = ab_capture_mask & not_totalseq_mask

    pmhc_gene_ids = list( raw.var[ get_gene_ids_varname(adata) ][ pmhc_mask ])

    pmhc_feature_names = list( raw.var_names[ pmhc_mask ] )

    pmhc_short_names = {}
    pmhc_long_names = {}

    for f, g in zip( pmhc_feature_names, pmhc_gene_ids ):
        long_name = g[:]
        if long_name.startswith('NR('):
            long_name = long_name[3:].replace(')','')
        assert long_name[0] in 'AB'
        hla,pep = long_name.split('_')[:2]
        short_name = '{}_{}{}'.format( hla[:3], pep[:4], len(pep))

        pmhc_short_names[f] = short_name
        pmhc_long_names[f]  = long_name

    return pmhc_short_names, pmhc_long_names


def shorten_pmhc_var_names(adata):
    pmhc_short_names ,_ = get_pmhc_short_and_long_names_dicts(adata)

    for ad in [adata, adata.raw]:
        if ad is None:
            continue

        names = list(ad.var.index)
        changed=False

        for ii in range(len(names)):
            if names[ii] in pmhc_short_names:
                changed=True
                print('change var_name from {} to {}'.format( names[ii], pmhc_short_names[names[ii]] ))
                names[ii] = pmhc_short_names[names[ii]]
        if changed:
            ad.var.index = names



def _get_X_pmhc( adata, pmhc_var_names ):
    ''' Returns:  X_pmhc
    Does not log the raw array so do this beforehand if you need it logged

    Not for use after conga adata has been setup; at that point use adata.obsm['X_pmhc']
    '''

    raw = adata.raw
    if raw is None:
        raw = adata

    var_names = list( raw.var_names )
    pmhc_indices = np.array( [ var_names.index(x) for x in pmhc_var_names ] )

    X_pmhc = raw[:,pmhc_indices].X.toarray()

    assert X_pmhc.shape == ( adata.shape[0], len(pmhc_var_names ) )

    return X_pmhc


# def get_pmhc_mask( ipmhc, threshold_type, umi_threshold, X_pmhc, X_pmhc_sorted, X_pmhc_argsorted ):

#     if threshold_type=='a': # absolute threshold
#         threshold = np.log( 1 + umi_threshold - 0.5 )
#         pmhc_mask = ( X_pmhc[:,ipmhc] >= threshold )
#     else:
#         top_pmhc_index = X_pmhc_argsorted[:,0]
#         log1p_delta = X_pmhc_sorted[:,0] - X_pmhc_sorted[:,1]
#         actual_delta = np.expm1(X_pmhc_sorted[:,0]) - np.expm1(X_pmhc_sorted[:,1])
#         pmhc_mask = ( ( top_pmhc_index == ipmhc ) & ( log1p_delta >= float(umi_threshold)/10 ) &
#                       ( actual_delta >= min_actual_delta_for_r_thresholds ) )
#     return pmhc_mask



# def setup_top_pmhcs( threshold_type, umi_threshold, pmhcs, tcrs, tcr2X_pmhc ):
#     ''' Return an np.array of strings, either 'none' or the pmhc for that clonotype (fullname, not shortname)
#     '''
#     assert threshold_type=='r' # otherwise it's not well defined

#     X_pmhc = np.array( [ tcr2X_pmhc[ x ] for x in tcrs ] )
#     X_pmhc_argsorted = np.argsort( -1*X_pmhc, axis=1 ) # negative since we want the biggest first
#     X_pmhc_sorted = -1*np.sort( -1*X_pmhc, axis=1 )
#     assert X_pmhc.shape == (len(tcrs),len(pmhcs))

#     top_pmhcs = ['none'] * len(tcrs)

#     for ipmhc, pmhc in enumerate( pmhcs ):
#         pmhc_mask = get_pmhc_mask( ipmhc, threshold_type, umi_threshold, X_pmhc, X_pmhc_sorted, X_pmhc_argsorted )

#         for ii, (tcr,m) in enumerate( zip( tcrs, pmhc_mask ) ):
#             if m:
#                 assert top_pmhcs[ii] == 'none'
#                 top_pmhcs[ii] = pmhc

#     return np.array( top_pmhcs )

def product_cdf(p1,p2):
    x = p1*p2
    if x<=0:
        return 0
    else:
        try:
            return max(0, x - x*np.log(x) )
        except:
            print('ERROR in product_cdf!',x)
            return x

def calc_sf_max( m, P ): # ~= 1 - ( 1-P)^m
    ''' If P is the sf of a random var X, return the sf of the max of m independent Xs
    looks like  1 - ( 1-P)^m = m*P - binom(m,2)*P^2 + binom(m,3)*P^3 ... +(-P)^m
    '''
    tot = 0.0
    for i in range(1,m+1):
        term = -1 * binom(m,i) * (-1*P)**i
        tot += term
        if tot and abs(term)<1e-3*abs(tot):
            break
    return tot

def calc_pmhc_nbrs_total_pval( pmhc_mask_in, nbrs, agroups, bgroups, verbose=False ):
    ''' at each step, eliminate nbrs of max-overlap cell as well as same-group cells
    '''
    pmhc_mask = np.copy( pmhc_mask_in )
    num_neighbors = nbrs.shape[1]

    min_combo_pval = 1.0
    combo_pval = None
    counter=0
    while True:
        counter+=1
        num_pos_cells = np.sum( pmhc_mask )
        if not num_pos_cells:
            break
        # nbrs = nbrs_in[pmhc_mask,:]
        # agroups = agroups_in[pmhc_mask,:]
        # bgroups = bgroups_in[pmhc_mask,:]

        # find the cell with the greatest number of nbrs
        max_overlap, ii_max = 0, None
        for ii,m in enumerate(pmhc_mask):
            if m:
                overlap = np.sum( pmhc_mask[ nbrs[ii,:] ] )
                if overlap>max_overlap:
                    max_overlap = overlap
                    ii_max = ii

        if max_overlap==0:
            # should we add a contribution here? like another product_cdf with 1.0? otherwise seems like there
            # might be bias...
            break

        # what are the odds of seeing this many nbrs?
        ii = ii_max

        not_same_group = ~( ( agroups==agroups[ii] ) | ( bgroups==bgroups[ii] ) )
        possible_pmhc_pos_nbrs = np.sum( pmhc_mask[ not_same_group ] )
        possible_nbrs = np.sum( not_same_group )

        expected = float(possible_pmhc_pos_nbrs*num_neighbors)/possible_nbrs
        if max_overlap<expected:
            break

        pval = hypergeom.sf( max_overlap-1, possible_nbrs, possible_pmhc_pos_nbrs, num_neighbors )
        if pval*num_pos_cells>1.0:
            break

        if pval==0:
            return pval

        sf_pval = calc_sf_max( num_pos_cells, pval )
        if sf_pval<0:
            print('ERROR lo_sf_pval1:', sf_pval, pval, num_pos_cells, max_overlap,
                  possible_nbrs, possible_pmhc_pos_nbrs, num_neighbors )
            break

        if combo_pval is None: ## first time through
            combo_pval = sf_pval
        else:
            combo_pval = product_cdf( combo_pval, sf_pval )

        if combo_pval==0:
            return combo_pval

        min_combo_pval = min(min_combo_pval, combo_pval)

        if verbose:
            print('calc_pmhc_nbrs_total_pval: {:2d} N= {} k= {} m= {} max_nbrs= {} pval= {:.1e} {:.1e} total_pval= {:.1e}'\
                  .format( counter, pmhc_mask.shape[0], num_neighbors, num_pos_cells, max_overlap, pval, sf_pval,
                           combo_pval ) )

        # remove the cell with the most nbrs, continue looping
        pmhc_mask[ nbrs[ii,:] ] = False
        pmhc_mask[ agroups==agroups[ii] ] = False
        pmhc_mask[ bgroups==bgroups[ii] ] = False


    sys.stdout.flush()
    return min_combo_pval


def compute_pmhc_versus_nbrs(
        adata,
        nbrs,
        agroups,
        bgroups,
        min_log1p_delta=2.0,
        min_actual_delta=3,
        min_positive_clones=3
):
    pmhc_var_names = adata.uns['pmhc_var_names']
    num_clones = adata.shape[0]

    X_pmhc = adata.obsm['X_pmhc']

    X_pmhc_sorted = -1 * np.sort( -1 * X_pmhc, axis=1 ) # in decreasing order
    X_pmhc_argsorted = np.argsort( -1 * X_pmhc, axis=1 ) # ditto

    results = []
    for ip, pmhc in enumerate(pmhc_var_names):
        top_pmhc_index = X_pmhc_argsorted[:,0]
        log1p_delta = X_pmhc_sorted[:,0] - X_pmhc_sorted[:,1]
        actual_delta = np.expm1(X_pmhc_sorted[:,0]) - np.expm1(X_pmhc_sorted[:,1])
        pmhc_mask = ( ( top_pmhc_index == ip ) & ( log1p_delta >= min_log1p_delta ) &
                      ( actual_delta >= min_actual_delta ) )

        num_positive_clones = np.sum( pmhc_mask )
        if num_positive_clones < min_positive_clones:
            continue ############## NOTE

        # find the total number of nbrs and the max nbrs per cell
        total_nbrs, max_nbrs, expected_total_nbrs = 0, 0, 0
        for ind1 in ( x for x,y in enumerate( pmhc_mask ) if y ):
            ind1_nbrs=0#pmhc pos nbrs that is
            diff_group_mask = (agroups!=agroups[ind1]) & (bgroups!=bgroups[ind1])
            # if len(nbrs[ind11]) nbrs are chosen at random from among the diff_group_mask clones, how many
            #  would we expect to be pmhc positive by chance?
            expected_total_nbrs += len(nbrs[ind1]) * np.sum( pmhc_mask & diff_group_mask) / np.sum(diff_group_mask)
            for ind2 in nbrs[ind1]:
                if pmhc_mask[ind2]:
                    assert ind2 != ind1
                    ind1_nbrs += 1
            max_nbrs = max( max_nbrs, ind1_nbrs )
            total_nbrs += ind1_nbrs

        expected_total_nbrs = max(1e-6, expected_total_nbrs) # no div by zero

        total_pval = calc_pmhc_nbrs_total_pval( pmhc_mask, nbrs, agroups, bgroups )

        if total_nbrs==0:
            if expected_total_nbrs <1.0:
                log2ratio = 0.0
            else:
                # give a pseudocount of 0.5
                log2ratio = np.log2(0.5/expected_total_nbrs)
        else:
            log2ratio = np.log2(float(total_nbrs)/expected_total_nbrs)

        results.append( dict(total_nbrs=total_nbrs,
                             expected_total_nbrs=expected_total_nbrs,
                             max_nbrs=max_nbrs,
                             log2_enrich=log2ratio,
                             pvalue=total_pval,
                             num_positive_clones=num_positive_clones,
                             pmhc=pmhc ) )
        # print('pmhc_{}_nbrs {:7d} {:3d} l2r {:7.3f} totP {:9.1e} Npos {:3d} {}'\
        #       .format(prefix_tag, total_nbrs, max_nbrs, log2ratio, total_pval,
        #               num_positive_clones, pmhc ))

    results_df = pd.DataFrame(results)
    return results_df

def calc_clone_pmhc_pvals(adata, min_log1p_delta=2.0, min_actual_delta=3 ):
    ''' This needs to be called before we subset to a single cell per clone
    '''
    pmhc_var_names = adata.uns['pmhc_var_names']

    pp.normalize_and_log_the_raw_matrix(adata) # just in case

    X_pmhc = get_X_pmhc( adata, pmhc_var_names )

    N = adata.shape[0] # total num cells, bigger than num clones
    assert X_pmhc.shape[0] == N
    X_pmhc_sorted = -1 * np.sort( -1 * X_pmhc )
    X_pmhc_argsorted = np.argsort( -1 * X_pmhc )
    log1p_delta = X_pmhc_sorted[:,0] - X_pmhc_sorted[:,1]
    actual_delta = np.expm1(X_pmhc_sorted[:,0]) - np.expm1(X_pmhc_sorted[:,1])

    is_pmhc_pos = ( actual_delta >= min_actual_delta ) & ( log1p_delta >= min_log1p_delta )
    top_pmhcs = np.array( [ pmhc_var_names[ X_pmhc_argsorted[ii,0] ] if m else 'none'
                            for ii,m in enumerate(is_pmhc_pos) ] )
    all_pmhc_counts = Counter( top_pmhcs )

    tcrs = pp.retrieve_tcrs_from_adata(adata, include_subject_id_if_present=True) # may contain duplicates
    unique_tcrs = sorted(set(tcrs))
    num_clones = len(unique_tcrs)
    tcr2clone_id = { y:x for x,y in enumerate(unique_tcrs)}
    clone_ids = np.array( [ tcr2clone_id[x] for x in tcrs ] )

    results = []

    print('calc_clone_pmhc_pvals: num_clones=', len(unique_tcrs))

    for c, tcr in enumerate(unique_tcrs):
        #clone_cells = np.array( [ x for x,y in enumerate(clone_ids) if y==c] )
        clone_mask = (clone_ids==c)
        clone_size = np.sum( clone_mask )
        X_pmhc_clone = X_pmhc[clone_mask,:]
        X_pmhc_clone_avg = np.sum(X_pmhc_clone, axis=0 )/clone_size
        X_pmhc_clone_avg_sorted = -1 * np.sort( -1 * X_pmhc_clone_avg )
        pmhc_counts = Counter( top_pmhcs[ clone_mask ] )
        for pmhc, count in pmhc_counts.items():
            if pmhc=='none': continue
            # surprise at seeing this many?
            # this is not be quite right since there are multiple possible pmhcs...
            pval = N*hypergeom.sf( count-1, N, clone_size, all_pmhc_counts[pmhc] )
            expected = float(clone_size*all_pmhc_counts[pmhc])/N
            l2r = np.log2( count / expected )
            # how would this pmhc look with an "r20" clone theshold??
            avglog1p_delta = X_pmhc_clone_avg[pmhc_var_names.index(pmhc)] - X_pmhc_clone_avg_sorted[1]

            if pval<1 or avglog1p_delta>= min_log1p_delta:
                results.append( {'clone_index':c,
                                 'adjusted_pvalue': pval,
                                 'pmhc_positive_fraction': float(count)/clone_size,
                                 'avglog1p_delta': avglog1p_delta,
                                 'num_pmhc_positive_in_clone': count,
                                 'clone_size': clone_size,
                                 'pmhc': pmhc,
                                 'num_pmhc_positive_overall': all_pmhc_counts[pmhc],
                                 'num_clones': num_clones,
                                 'va': tcr[0][0],
                                 'ja': tcr[0][1],
                                 'cdr3a': tcr[0][2],
                                 'cdr3a_nucseq': tcr[0][3],
                                 'vb': tcr[1][0],
                                 'jb': tcr[1][1],
                                 'cdr3b': tcr[1][2],
                                 'cdr3b_nucseq': tcr[1][3],
                             } )


                # print('clone_pval: {:9.2e} {:.3f} al1pd {:6.2f} l2r {:6.2f} {:4d} {:4d} {:5d} {:5d} {}'\
                #       .format(pval, float(count)/clone_size,  avglog1p_delta, l2r,
                #               count, clone_size, all_pmhc_counts[pmhc], N, pmhc))

    results_df = pd.DataFrame( results )
    return results_df
