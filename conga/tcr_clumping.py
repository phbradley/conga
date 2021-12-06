######################################################################################88
# import scanpy as sc
# import random
import pandas as pd
from os.path import exists
from pathlib import Path
from collections import OrderedDict, Counter
# from sklearn.metrics import pairwise_distances
# from sklearn.utils import sparsefuncs
# from sklearn.decomposition import KernelPCA
import numpy as np
# import scipy
# from scipy.cluster import hierarchy
# from scipy.spatial.distance import squareform, cdist
# from scipy.sparse import issparse#, csr_matrix
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
# from anndata import AnnData
import sys
import os
from sys import exit
from . import util
from . import preprocess
from .tags import *
from .tcrdist import tcr_sampler
import random


def estimate_background_tcrdist_distributions(
        organism,
        tcrs,
        max_dist,
        num_random_samples = 50000,
        pseudocount = 0.25,
        tmpfile_prefix = None,
        background_alpha_chains = None, # default is to get these by shuffling tcrs_for_background_generation
        background_beta_chains = None, #  -- ditto --
        tcrs_for_background_generation = None, # default is to use 'tcrs'
        preserve_vj_pairings = False,
        save_unpaired_dists = False, # tell C++ code to save to files (dev feature)
        nocleanup = False,
):
    ''' Returns a numpy float matrix P of shape (len(tcrs), max_dist+1)

    the (i,j) entry in P is an estimate of the probability of seeing a paired tcrdist
    score <= j for tcr #i.

    '''
    if not util.tcrdist_cpp_available():
        print('conga.tcr_clumping.estimate_background_tcrdist_distributions:: need to compile the C++ tcrdist executables')
        exit(1)

    if tmpfile_prefix is None:
        tmpfile_prefix = Path('./tmp_nbrs{}'.format(random.randrange(1,10000)))
    else:
        tmpfile_prefix = Path(tmpfile_prefix)


    if tcrs_for_background_generation is None:
        # only used when background_alpha_chains and/or background_beta_chains is None
        #tcrs_for_background_generation = tcrs
        # since 10x doesn't always provide allele information, we need to try out alternate
        # alleles to get the best parses...
        tcrs_for_background_generation = tcr_sampler.find_alternate_alleles_for_tcrs(
            organism, tcrs, verbose=False)

    max_dist = int(0.1+max_dist) ## need an integer

    if background_alpha_chains is None or background_beta_chains is None:
        # parse the V(D)J junction regions of the tcrs to define split-points for shuffling
        junctions_df = tcr_sampler.parse_tcr_junctions(organism, tcrs_for_background_generation)

        # resample shuffled single-chain tcrs
        if background_alpha_chains is None:
            background_alpha_chains = tcr_sampler.resample_shuffled_tcr_chains(
                organism, num_random_samples, 'A', junctions_df,
                preserve_vj_pairings = preserve_vj_pairings)
        if background_beta_chains is None:
            background_beta_chains  = tcr_sampler.resample_shuffled_tcr_chains(
                organism, num_random_samples, 'B', junctions_df,
                preserve_vj_pairings = preserve_vj_pairings)

    # save all tcrs to files
    achains_file = str(tmpfile_prefix) + '_bg_achains.tsv'
    bchains_file = str(tmpfile_prefix) + '_bg_bchains.tsv'
    tcrs_file = str(tmpfile_prefix) + '_tcrs.tsv'

    pd.DataFrame({'va'   :[x[0] for x in background_alpha_chains],
                  'cdr3a':[x[2] for x in background_alpha_chains]}).to_csv(
                      achains_file, sep='\t', index=False)

    pd.DataFrame({'vb'   :[x[0] for x in background_beta_chains ],
                  'cdr3b':[x[2] for x in background_beta_chains ]}).to_csv(
                      bchains_file, sep='\t', index=False)

    pd.DataFrame({'va':[x[0][0] for x in tcrs], 'cdr3a':[x[0][2] for x in tcrs],
                  'vb':[x[1][0] for x in tcrs], 'cdr3b':[x[1][2] for x in tcrs]})\
      .to_csv(tcrs_file, sep='\t', index=False)


    # compute distributions vs background chains
    if os.name == 'posix':
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) ,
                             'calc_distributions')
    else:
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) ,
                             'calc_distributions.exe')

    outfile = str(tmpfile_prefix) + '_dists.txt'

    db_filename = Path.joinpath( Path(util.path_to_tcrdist_cpp_db) ,
                                 'tcrdist_info_{}.txt'.format( organism))

    cmd = '{} -f {} -m {} -d {} -a {} -b {} -o {}'\
          .format(exe, tcrs_file, max_dist, db_filename,
                  achains_file, bchains_file, outfile)

    if save_unpaired_dists:
        cmd += ' -u'

    util.run_command(cmd, verbose=True)

    if not exists(outfile):
        print('tcr_clumping:: calc_distributions failed: missing', outfile)
        exit(1)

    counts = np.loadtxt(outfile, dtype=int)
    counts = np.cumsum(counts, axis=1)
    assert counts.shape == (len(tcrs), max_dist+1)
    n_bg_pairs = len(background_alpha_chains) * len(background_beta_chains)
    tcrdist_freqs = np.maximum(pseudocount, counts.astype(float))/n_bg_pairs

    if not nocleanup:
        for filename in [achains_file, bchains_file, tcrs_file, outfile]:
            os.remove(filename)

    return tcrdist_freqs

# not to be confused with assess_tcr_clumping which takes in an adata
# and is basically a wrapper around this function
def find_tcr_clumping(
        tcrs,
        organism,
        tmpfile_prefix = 'tmp',
        radii = [24, 48, 72, 96],
        num_random_samples = 50000, # higher ==> more sensitive, but slower
        pvalue_threshold = 1.0,
        verbose=True,
        clusters_gex=None, # if passed, will look for clumps within clusters
        preserve_vj_pairings = False,
        bg_tcrs = None, # usually better to leave this None
):
    ''' Returns a pandas dataframe with the following columns:
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
    if not util.tcrdist_cpp_available():
        print('conga.tcr_clumping.find_tcr_clumping::',
              'need to compile the C++ tcrdist executables')
        exit(1)

    num_clones = len(tcrs)

    if bg_tcrs is None:
        bg_tcrs = tcrs

    radii = [int(x+0.1) for x in radii] #ensure integers

    outprefix = f'{tmpfile_prefix}_{random.random()}_tcr_clumping'
    tmpfiles = []

    # this will optimize V/J alleles to fit the cdr3_nucseqs better
    # in case 10x didn't give us the alleles
    bg_freqs = estimate_background_tcrdist_distributions(
        organism, tcrs, max(radii),
        num_random_samples=num_random_samples,
        tmpfile_prefix=outprefix,
        preserve_vj_pairings=preserve_vj_pairings,
        tcrs_for_background_generation=bg_tcrs,
    )

    tcrs_file = outprefix +'_tcrs.tsv'
    pd.DataFrame({
        'va'   : [x[0][0] for x in tcrs],
        'cdr3a': [x[0][2] for x in tcrs],
        'vb'   : [x[1][0] for x in tcrs],
        'cdr3b': [x[1][2] for x in tcrs],
    }).to_csv(tcrs_file, sep='\t', index=False)
    tmpfiles.append(tcrs_file)

    # find neighbors in fg tcrs up to max(radii) ############################

    if os.name == 'posix':
        exe = Path.joinpath(
            Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors')
    else:
        exe = Path.joinpath(
            Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors.exe')

    agroups, bgroups = preprocess.setup_tcr_groups_for_tcrs(tcrs)
    agroups_filename = outprefix+'_agroups.txt'
    bgroups_filename = outprefix+'_bgroups.txt'
    np.savetxt(agroups_filename, agroups, fmt='%d')
    np.savetxt(bgroups_filename, bgroups, fmt='%d')
    tmpfiles.extend([agroups_filename, bgroups_filename])

    db_filename = Path.joinpath(
        Path(util.path_to_tcrdist_cpp_db), f'tcrdist_info_{organism}.txt')

    tcrdist_threshold = max(radii)

    cmd = '{} -f {} -t {} -d {} -o {} -a {} -b {}'\
          .format(exe, tcrs_file, tcrdist_threshold, db_filename, outprefix,
                  agroups_filename, bgroups_filename)

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
    for line1, line2 in zip(open(nbr_indices_filename,'r'),
                            open(nbr_distances_filename,'r')):
        l1 = line1.split()
        l2 = line2.split()
        assert len(l1) == len(l2)
        #ii = len(all_nbrs)
        all_nbrs.append([int(x) for x in l1])
        all_distances.append([int(x) for x in l2])
    assert len(all_nbrs) == num_clones

    # cleanup the tmpfiles
    for tmpfile in tmpfiles:
        if exists(tmpfile):
            os.remove(tmpfile)

    # we were printing this out in verbose mode...
    #clone_sizes = adata.obs['clone_sizes']

    # use poisson to find nbrhoods with more tcrs than expected;
    #  have to handle agroups/bgroups
    dfl = []

    is_clumped = np.full((num_clones,), False)

    n_bg_pairs = num_random_samples * num_random_samples

    all_raw_pvalues = np.full((num_clones, len(radii)), 1.0)

    for ii in range(num_clones):
        ii_freqs = bg_freqs[ii]
        ii_dists = all_distances[ii]
        for irad, radius in enumerate(radii):
            num_nbrs = sum(x<=radius for x in ii_dists)
            if num_nbrs<1:
                continue # NOTE: OK since wont have intra-cluster nbrs either
            max_nbrs = np.sum((agroups!=agroups[ii]) & (bgroups!=bgroups[ii]))
            # adjust for number of tests
            mu = max_nbrs * ii_freqs[radius]
            pval = poisson.sf(num_nbrs-1, mu)
            all_raw_pvalues[ii, irad] = pval
            pval *= len(radii) * num_clones # simple multiple test correction
            if pval <= pvalue_threshold:
                is_clumped[ii] = True
                # count might just be pseudocount
                raw_count = ii_freqs[radius]*n_bg_pairs
                if verbose:
                    atcr_str = ' '.join(tcrs[ii][0][:3])
                    btcr_str = ' '.join(tcrs[ii][1][:3])
                    print(f'tcr_nbrs_global: {num_nbrs:2d} {mu:9.6f}',
                          f'radius: {radius:2d} pval: {pval:9.1e}',
                          f'{raw_count:9.1f} tcr: {atcr_str} {btcr_str}')

                dfl.append( OrderedDict(
                    clump_type='global',
                    clone_index=ii,
                    nbr_radius=radius,
                    pvalue_adj=pval,
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
            if clusters_gex is not None:
                # look for clumping within the GEX cluster containing ii
                ii_nbrs = all_nbrs[ii]
                ii_cluster = clusters_gex[ii]
                ii_cluster_mask = (clusters_gex==ii_cluster)
                num_nbrs = sum((x<=radius and
                                clusters_gex[y]==ii_cluster)
                               for x,y in zip(ii_dists, ii_nbrs))
                if num_nbrs<1:
                    continue ## NOTE-- continue
                max_nbrs = np.sum((agroups!=agroups[ii]) &
                                  (bgroups!=bgroups[ii]) &
                                  ii_cluster_mask)
                mu = max_nbrs * ii_freqs[radius]
                pval = (len(radii) * num_clones *
                        poisson.sf(num_nbrs-1, mu ))
                if pval <= pvalue_threshold:
                    is_clumped[ii] = True
                    # if count was 0, will be pseudocount
                    raw_count = ii_freqs[radius]*n_bg_pairs
                    if verbose:
                        atcr_str = ' '.join(tcrs[ii][0][:3])
                        btcr_str = ' '.join(tcrs[ii][1][:3])
                        print(f'tcr_nbrs_intra: {num_nbrs:2d} {mu:9.6f}',
                              f'radius: {radius:2d} pval: {pval:9.1e}',
                              f'{raw_count:9.1f} tcr: {atcr_str} {btcr_str}')
                    dfl.append( OrderedDict(
                        clump_type='intra_gex_cluster',
                        clone_index=ii,
                        nbr_radius=radius,
                        pvalue_adj=pval,
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
        return results_df ## NOTE EARLY RETURN #########################

    # compute FDR values in addition to the simple adjusted pvalues
    _, fdr_values, _, _ = multipletests(
        all_raw_pvalues.reshape(-1), alpha=0.05, method='fdr_bh')
    fdr_values = fdr_values.reshape((num_clones, len(radii))).min(axis=1)
    # right now we don't assign fdr values for intra-gex cluster clumping
    fdr_column = [fdr_values[x.clone_index] if x.clump_type=='global'
                  else np.nan for x in results_df.itertuples()]
    results_df['clonotype_fdr_value'] = fdr_column

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

    return results_df ## there is an early return if results_df is empty



def assess_tcr_clumping(
        adata,
        also_find_clumps_within_gex_clusters = False,
        outfile_prefix = None, # for saving the table to tsv
        verbose = False,
        **kwargs,
        # these are the possible kwargs: see find_tcr_clumping above
        #tmpfile_prefix = 'tmp'
        #radii = [24, 48, 72, 96],
        #num_random_samples = 50000,
        #pvalue_threshold = 1.0,
):
    ''' This is a wrapper around find_tcr_clumping that takes an AnnData

    Returns a pandas dataframe with the following columns:
    - clone_index
    - nbr_radius
    - pvalue_adj
    - num_nbrs
    - expected_num_nbrs
    - raw_count
    - va, ja, cdr3a, vb, jb, cdr3b (ie, the 6 tcr cols for clone_index clone)
    - clumping_group: clonotypes within each other's significant nbr_radii are linked
    - clump_type: string, either 'global' or 'intra_gex_cluster' (latter only if also_find_clumps_within_gex_clusters=T)

    Wrapper around find_tcr_clumping; here we take adata and unpack the info

    returns the results dataframe and also stores them in the

    adata.uns['conga_results'][conga.tags.TCR_CLUMPING]

    '''
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    if also_find_clumps_within_gex_clusters:
        clusters_gex = np.array(adata.obs['clusters_gex'])
    else:
        clusters_gex = None

    organism = adata.uns['organism']

    results = find_tcr_clumping(
        tcrs,
        organism,
        clusters_gex=clusters_gex,
        verbose=verbose,
        **kwargs
    )

    ## store best pvalues in adata.obs:
    num_clones = adata.shape[0]
    tcr_clumping_pvalues = np.full((num_clones,), num_clones).astype(float)
    for l in results.itertuples():
        tcr_clumping_pvalues[l.clone_index] = min(
            tcr_clumping_pvalues[l.clone_index], l.pvalue_adj)
    # stash in adata.obs
    adata.obs['tcr_clumping_pvalues'] = tcr_clumping_pvalues


    if not results.empty: ## add cluster info
        clusters_gex = np.array(adata.obs['clusters_gex']) # over-write old
        clusters_tcr = np.array(adata.obs['clusters_tcr'])
        results['clusters_gex'] = [clusters_gex[x] for x in results.clone_index]
        results['clusters_tcr'] = [clusters_tcr[x] for x in results.clone_index]
        results.sort_values('pvalue_adj', inplace=True)

    ## store results in adata.uns
    adata.uns.setdefault('conga_results', {})[TCR_CLUMPING] = results

    #
    help_message = """ This table stores the results of the TCR "clumping"
analysis, which looks for neighborhoods in TCR space with more TCRs than
expected by chance under a simple null model of VDJ rearrangement.

For each TCR in the dataset, we count how many TCRs are within a set of
fixed TCRdist radii (defaults: 24,48,72,96), and compare that number
to the expected number given the size of the dataset using the poisson
model. Inspired by the ALICE and TCRnet methods.

Columns:
clump_type='global' unless we are optionally looking for TCR clumps within
   the individual GEX clusters
num_nbrs = neighborhood size (number of other TCRs with TCRdist<nbr_radius)
pvalue_adj= raw Poisson P value of the neighborhood size multiplied by
   the number of TCRs in the dataset and the number of radii
clonotype_fdr_value= Benjamin-Hochberg FDR value for the per-TCR
   significance (ie, looking at all the radii).

    """

    adata.uns['conga_results'][TCR_CLUMPING+HELP_SUFFIX] = help_message

    if outfile_prefix is not None:
        util.save_table_and_helpfile(
            TCR_CLUMPING, adata, outfile_prefix)

    return results

def tcrs_from_dataframe_helper(df, add_j_and_nucseq=False):
    ''' Helper function that creates a tcrs list from a dataframe
    the dataframe should have columns va, cdr3a, vb, cdr3b

    if add_j_and_nucseq is True, then the dataframe should also have the
    columns ja, cdr3a_nucseq, jb, cdr3b_nucseq
    '''
    if add_j_and_nucseq:
        return [ ( (x.va, x.ja, x.cdr3a, x.cdr3a_nucseq.lower()),
                   (x.vb, x.jb, x.cdr3b, x.cdr3b_nucseq.lower()) )
                 for x in df.itertuples() ]
    else:
        return [ ( (x.va, None, x.cdr3a), (x.vb, None, x.cdr3b) )
                 for x in df.itertuples() ]

def find_significant_tcrdist_matches(
        query_tcrs_df,
        db_tcrs_df,
        organism,
        tmpfile_prefix = None,
        adjusted_pvalue_threshold = 1.0, # adjusted for query and db sizes
        background_tcrs_df = None, # default is to use query_tcrs_df
        num_random_samples_for_bg_freqs = 50000,
        nocleanup=False,
        fixup_allele_assignments_in_background_tcrs_df=True,
):
    ''' Computes paired tcrdist distances between query_tcrs_df and
    db_tcrs_df and converts them to
    to p-values, adjusted for both the number of query and db tcrs
    (ie, pvalue_adj = num_query_tcrs * num_db_tcrs * pvalue_raw)

    Anything ending in _df is a pandas DataFrame

    Required columns in query_tcrs_df and db_tcrs_df: va, cdr3a, vb, cdr3b

    Required columns in background_tcrs_df:
       va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq

    returns results pd.DataFrame with columns:

    * tcrdist
    * pvalue_adj
    * va, {ja}, cdr3a, vb, {jb}, cdr3b {}:if present in query_tcrs_df
    * PLUS all columns in db_tcrs_df prepended with 'db_' string

    NOTE: reported pvalues are adjusted for sizes of both query_tcrs_df and
       db_tcrs_df
       db_tcrs_tsvfile

    '''
    if tmpfile_prefix is None:
        tmpfile_prefix = 'tmpfile{}'.format(random.randrange(1,10000))

    if background_tcrs_df is None:
        background_tcrs_df = query_tcrs_df

    for tag, df in [ ['query', query_tcrs_df],
                     ['db', db_tcrs_df],
                     ['background', background_tcrs_df]]:
        for ab in 'ab':
            required_cols = f'cdr3{ab} v{ab}'.split()
            if tag == 'background':
                required_cols += f'cdr3{ab}_nucseq j{ab}'.split()
            for col in required_cols:
                if col not in df.columns:
                    print('ERROR find_significant_tcrdist_matches::',
                          f'{tag} df is missing {col} column')
                    return

    query_tcrs = tcrs_from_dataframe_helper(query_tcrs_df)
    db_tcrs = tcrs_from_dataframe_helper(db_tcrs_df)
    background_tcrs = tcrs_from_dataframe_helper(
        background_tcrs_df, add_j_and_nucseq=True)

    num_comparisons = len(query_tcrs) * len(db_tcrs)
    print('num_comparisons:', num_comparisons, len(query_tcrs), len(db_tcrs))

    if fixup_allele_assignments_in_background_tcrs_df:
        background_tcrs = tcr_sampler.find_alternate_alleles_for_tcrs(
            organism, background_tcrs, verbose=False)


    max_dist = 200

    bg_freqs = estimate_background_tcrdist_distributions(
        organism, query_tcrs, max_dist,
        num_random_samples= num_random_samples_for_bg_freqs,
        tcrs_for_background_generation= background_tcrs)
    assert bg_freqs.shape == (len(query_tcrs), max_dist+1)

    adjusted_bg_freqs = num_comparisons * bg_freqs

    could_match = np.any( adjusted_bg_freqs<= adjusted_pvalue_threshold, axis=0)
    assert could_match.shape == (max_dist+1,)

    max_dist_for_matching = 0 # must be some numpy way of doing this
    while could_match[max_dist_for_matching] and max_dist_for_matching<max_dist:
        max_dist_for_matching += 1

    print(f'find_significant_tcrdist_matches:: max_dist: {max_dist}',
          f'max_dist_for_matching: {max_dist_for_matching}')

    # now run C++ matching code
    query_tcrs_file = tmpfile_prefix+'temp_query_tcrs.tsv'
    db_tcrs_file = tmpfile_prefix+'temp_db_tcrs.tsv'
    query_tcrs_df['va cdr3a vb cdr3b'.split()].to_csv(query_tcrs_file, sep='\t',
                                                      index=False)
    db_tcrs_df['va cdr3a vb cdr3b'.split()].to_csv(db_tcrs_file, sep='\t',
                                                   index=False)

    if os.name == 'posix':
        exe = Path(util.path_to_tcrdist_cpp_bin) / 'find_paired_matches'
    else:
        exe = Path(util.path_to_tcrdist_cpp_bin) / 'find_paired_matches.exe'

    if not exists(exe):
        print(f'ERROR: find_paired_matches:: tcrdist_cpp executable {exe}',
              'is missing')
        print('ERROR: see instructions in github repository README for',
              'compiling tcrdist_cpp (tldr: type make in conga/tcrdist_cpp/')
        return

    db_filename = (Path(util.path_to_tcrdist_cpp_db) /
                   f'tcrdist_info_{organism}.txt')

    outfilename = tmpfile_prefix+'temp_tcr_matching.tsv'

    cmd = '{} -i {} -j {} -t {} -d {} -o {}'\
          .format(exe, query_tcrs_file, db_tcrs_file, max_dist_for_matching,
                  db_filename, outfilename)

    util.run_command(cmd, verbose=True)

    df = pd.read_csv(outfilename, sep='\t')

    num_matches = df.shape[0]
    raw_pvalues = [bg_freqs[x.index1][x.tcrdist] for x in df.itertuples()]
    raw_pvalues_argsort = np.argsort(raw_pvalues)
    raw_pvalues_sorted = np.concatenate(
        [np.sort(raw_pvalues), np.full((num_comparisons - num_matches,), 1.0)])
    assert raw_pvalues_sorted.shape[0] == num_comparisons #BIG
    # the alpha setting here does not matter for method='fdr_bh'
    _, fdr_values_sorted, _, _ = multipletests(
        raw_pvalues_sorted, alpha=0.05, is_sorted=True, method='fdr_bh')
    fdr_values = fdr_values_sorted[np.argsort(raw_pvalues_argsort)]

    dfl = []
    for l, fdr_value in zip(df.itertuples(), fdr_values):
        i = l.index1
        j = l.index2
        pvalue_adj = num_comparisons * bg_freqs[i][l.tcrdist]
        if pvalue_adj > adjusted_pvalue_threshold:
            continue
        query_row = query_tcrs_df.iloc[i]
        db_row = db_tcrs_df.iloc[j]
        assert query_row.cdr3b == l.cdr3b1 # sanity check
        assert db_row.cdr3b == l.cdr3b2 # ditto
        D = OrderedDict(tcrdist= l.tcrdist,
                        pvalue_adj=pvalue_adj,
                        fdr_value=fdr_value, # Benjamini-Hochberg
                        query_index= i,
                        db_index= j,
                        va= query_row.va,
                        #ja= query_row.ja,
                        cdr3a= query_row.cdr3a,
                        vb= query_row.vb,
                        #jb= query_row.jb,
                        cdr3b= query_row.cdr3b)
        if 'ja' in query_row:
            D['ja'] = query_row.ja
        if 'jb' in query_row:
            D['jb'] = query_row.jb
        for tag in db_row.index:
            D['db_'+tag] = db_row[tag]
        dfl.append(D)

    if not nocleanup:
        os.remove(outfilename)
        os.remove(query_tcrs_file)
        os.remove(db_tcrs_file)

    results_df = pd.DataFrame(dfl)
    if dfl:
        results_df.sort_values('pvalue_adj', inplace=True)
        for col in results_df.columns:
            # object columns with missing values cause problems in h5 writing
            if (np.sum(results_df[col].isna()) and
                results_df[col].dtype is np.dtype('O')):
                print('replacing nans in object column that has missing values',
                      col)
                results_df[col] = results_df[col].fillna('').astype(str)


    return results_df


def match_adata_tcrs_to_db_tcrs(
        adata,
        outfile_prefix=None,
        db_tcrs_tsvfile=None,
        tmpfile_prefix=None,
        adjusted_pvalue_threshold = 1.0, # adjusted for size of adata AND db_tcrs_tsvfile
        tcrs_for_background_generation = None, # default is to use tcrs from adata
        num_random_samples_for_bg_freqs = 50000,
        nocleanup=False,
        fixup_allele_assignments_in_background_tcrs_df=True,
):
    ''' Find significant tcrdist matches between tcrs in adata and tcrs
    in db_tcrs_tsvfile

    by calling the find_significant_tcrdist_matches function above (see that
    docstring too)

    returns results pd.DataFrame with columns:
       tcrdist, pvalue_adj, va, ja, cdr3a, vb, jb, cdr3b,
       PLUS all columns in db_tcrs_tsvfile prepended with 'db_' string

    db_tcrs_tsvfile has at a minimum the columns:
       va (or va_gene) cdr3a vb (or vb_gene) cdr3b

    NOTE: reported pvalues are adjusted for sizes of both adata and
       db_tcrs_tsvfile

    '''

    if db_tcrs_tsvfile is None:
        if adata.uns['organism'] != 'human':
            print('ERROR: match_adata_tcrs_to_db_tcrs db_tcrs_tsvfile is None')
            print('but we only have built-in database for organism=human')
            return pd.DataFrame() ##### NOTE EARLY RETURN HERE ################

        print('tcr_clumping.match_adata_tcrs_to_db_tcrs: Matching to',
              'default literature TCR database; for more info see',
              'conga/data/new_paired_tcr_db_for_matching_nr_README.txt')
        db_tcrs_tsvfile = Path.joinpath(
            util.path_to_data, 'new_paired_tcr_db_for_matching_nr.tsv')

    print('Matching to paired tcrs in', db_tcrs_tsvfile)

    query_tcrs_df = adata.obs['va ja cdr3a vb jb cdr3b'.split()].copy()
    db_tcrs_df = pd.read_csv(db_tcrs_tsvfile, sep='\t')

    # possibly swap legacy column names
    if 'va' not in db_tcrs_df.columns and 'va_gene' in db_tcrs_df.columns:
        db_tcrs_df['va'] = db_tcrs_df['va_gene']
    if 'vb' not in db_tcrs_df.columns and 'vb_gene' in db_tcrs_df.columns:
        db_tcrs_df['vb'] = db_tcrs_df['vb_gene']

    if tcrs_for_background_generation is None:
        background_tcrs_df = adata.obs['va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq'.split()]\
                                  .copy()
    else:
        background_tcrs_df = pd.DataFrame(
            [ dict(va=x[0], ja=x[1], cdr3a=x[2], cdr3a_nucseq=x[3],
                   vb=y[0], jb=y[1], cdr3b=y[2], cdr3b_nucseq=y[3])
              for x,y in tcrs_for_background_generation ])

    results = find_significant_tcrdist_matches(
        query_tcrs_df,
        db_tcrs_df,
        adata.uns['organism'],
        tmpfile_prefix=tmpfile_prefix,
        adjusted_pvalue_threshold=adjusted_pvalue_threshold,
        background_tcrs_df=background_tcrs_df,
        num_random_samples_for_bg_freqs=num_random_samples_for_bg_freqs,
        fixup_allele_assignments_in_background_tcrs_df=fixup_allele_assignments_in_background_tcrs_df,
        nocleanup=nocleanup,
        )

    if not results.empty:
        results.rename(columns={'query_index':'clone_index'}, inplace=True)
        barcodes = list(adata.obs.index)
        results['barcode'] = [barcodes[x] for x in results.clone_index]
        results.sort_values(['pvalue_adj','tcrdist'], inplace=True)


    # store results in adata.uns
    table_tag = TCR_DB_MATCH
    adata.uns.setdefault('conga_results', {})[table_tag] = results

    help_message = f"""This table stores significant matches between
TCRs in adata and TCRs in the file {db_tcrs_tsvfile}

P values of matches are assigned by turning the raw TCRdist
score into a P value based on a model of the V(D)J rearrangement
process, so matches between TCRs that are very far from germline
(for example) are assigned a higher significance.

Columns:

tcrdist: TCRdist distance between the two TCRs (adata query and db hit)

pvalue_adj: raw P value of the match * num query TCRs * num db TCRs

fdr_value: Benjamini-Hochberg FDR value for match

clone_index: index within adata of the query TCR clonotype

db_index: index of the hit in the database being matched

va,ja,cdr3a,vb,jb,cdr3b

db_XXX: where XXX is a field in the literature database

    """

    # store help info in adata.uns
    adata.uns['conga_results'][table_tag+HELP_SUFFIX] = help_message

    if outfile_prefix is not None:
        util.save_table_and_helpfile(table_tag, adata, outfile_prefix)

    return results



def strict_single_chain_match_adata_tcrs_to_db_tcrs(
        adata,
        outfile_prefix=None,
        db_tcrs_tsvfile=None,
):
    ''' Find CDR3a and CDR3b matches between adata tcrs and tcrs in db_tcrs_tsvfile
    based on strictly matching amino acid sequences only

    returns a list containing two DataFrames; 0 is alpha chain hits, 1 is beta chain hits.
    DataFrames contain a columns with the UMI barcodes, GEX clusters, and TCR clusters of matching clonotypes

    db_tcrs_tsvfile has at a minimum the columns: cdr3a and cdr3b

    '''

    if db_tcrs_tsvfile is None:
        if adata.uns['organism'] == 'human':
            db_tcrs_tsvfile = Path.joinpath(util.path_to_data, 'human_tcr_db_for_matching.tsv')
        elif adata.uns['organism'] == 'mouse':
            db_tcrs_tsvfile = Path.joinpath(util.path_to_data, 'mouse_tcr_db_for_matching.tsv')
        print('tcr_clumping.strict_single_chain_match_adata_tcrs_to_db_tcrs: Matching to default literature TCR database; for more info see conga/data/human_and_mouse_tcr_db_for_matching_README.txt')

    print('Matching to CDR3a and CDR3b sequences in', db_tcrs_tsvfile)

    query_tcrs_df = adata.obs['va ja cdr3a vb jb cdr3b'.split()].copy()
    db_tcrs_df = pd.read_csv(db_tcrs_tsvfile, sep='\t')

    # possibly swap legacy column names.
    if 'va' not in db_tcrs_df.columns and 'va_gene' in db_tcrs_df.columns:
        db_tcrs_df['va'] = db_tcrs_df['va_gene']
    if 'vb' not in db_tcrs_df.columns and 'vb_gene' in db_tcrs_df.columns:
        db_tcrs_df['vb'] = db_tcrs_df['vb_gene']

    # generate a list of df with database CDR3a or CDR3b matching to adata.obs
    matched_dfs = []
    chains = ('cdr3a', 'cdr3b')
    for i in range(len(chains)):

        chain = chains[i]

        matched_dfs.append( db_tcrs_df[ db_tcrs_df[chain].isin(query_tcrs_df[chain]) ].copy() )

        matched_dfs[i] = matched_dfs[i].rename(columns = {'Unnamed: 0': 'db_index'})

        if matched_dfs[i].empty:
            print(f'No {chain} matches detected')
        else:
            clones = []
            gex_clusters = []
            tcr_clusters = []

            for (idx, row) in matched_dfs[i].iterrows():

                clone_bc = adata.obs.index[adata.obs[chain] == row.loc[chain] ].to_list()
                clone_gex = adata.obs.clusters_gex[adata.obs[chain] == row.loc[chain] ].astype(str).unique().tolist()
                clone_tcr = adata.obs.clusters_tcr[adata.obs[chain] == row.loc[chain] ].astype(str).unique().tolist()

                clones.append( ",".join(clone_bc) )
                gex_clusters.append( ",".join(clone_gex) )
                tcr_clusters.append( ",".join(clone_tcr) )

            matched_dfs[i][f'{chain}_match_UMI'] = clones
            matched_dfs[i][f'{chain}_match_gex_clusters'] = gex_clusters
            matched_dfs[i][f'{chain}_match_tcr_clusters'] = tcr_clusters

        if outfile_prefix is not None:
            csv_file = outfile_prefix + f'_single_chain_db_matches_{chain}.csv'
            matched_dfs[i].to_csv(csv_file, index=False)

    return matched_dfs
