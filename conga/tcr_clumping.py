# import scanpy as sc
# import random
import pandas as pd
from os.path import exists
from pathlib import Path
from collections import OrderedDict#, Counter
# from sklearn.metrics import pairwise_distances
# from sklearn.utils import sparsefuncs
# from sklearn.decomposition import KernelPCA
import numpy as np
# import scipy
# from scipy.cluster import hierarchy
# from scipy.spatial.distance import squareform, cdist
# from scipy.sparse import issparse#, csr_matrix
from scipy.stats import poisson
# from anndata import AnnData
import sys
import os
from sys import exit
from . import util
from . import preprocess
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
):
    if not util.tcrdist_cpp_available():
        print('conga.tcr_clumping.estimate_background_tcrdist_distributions:: need to compile the C++ tcrdist executables')
        exit(1)

    if tmpfile_prefix is None:
        tmpfile_prefix = Path('./tmp_nbrs{}'.format(random.randrange(1,10000)))
    else:
        tmpfile_prefix = Path(tmpfile_prefix)


    if tcrs_for_background_generation is None:
        # only used when background_alpha_chains and/or background_beta_chains is None
        tcrs_for_background_generation = tcrs

    max_dist = int(0.1+max_dist) ## need an integer

    if background_alpha_chains is None or background_beta_chains is None:
        # parse the V(D)J junction regions of the tcrs to define split-points for shuffling
        junctions_df = tcr_sampler.parse_tcr_junctions(organism, tcrs_for_background_generation)

        # resample shuffled single-chain tcrs
        if background_alpha_chains is None:
            background_alpha_chains = tcr_sampler.resample_shuffled_tcr_chains(
                organism, num_random_samples, 'A', junctions_df)
        if background_beta_chains is None:
            background_beta_chains  = tcr_sampler.resample_shuffled_tcr_chains(
                organism, num_random_samples, 'B', junctions_df)

    # save all tcrs to files
    achains_file = str(tmpfile_prefix) + '_bg_achains.tsv'
    bchains_file = str(tmpfile_prefix) + '_bg_bchains.tsv'
    tcrs_file = str(tmpfile_prefix) + '_tcrs.tsv'

    pd.DataFrame({'va'   :[x[0] for x in background_alpha_chains],
                  'cdr3a':[x[2] for x in background_alpha_chains]}).to_csv(achains_file, sep='\t', index=False)

    pd.DataFrame({'vb'   :[x[0] for x in background_beta_chains ],
                  'cdr3b':[x[2] for x in background_beta_chains ]}).to_csv(bchains_file, sep='\t', index=False)

    pd.DataFrame({'va':[x[0][0] for x in tcrs], 'cdr3a':[x[0][2] for x in tcrs],
                  'vb':[x[1][0] for x in tcrs], 'cdr3b':[x[1][2] for x in tcrs]})\
      .to_csv(tcrs_file, sep='\t', index=False)


    # compute distributions vs background chains
    if os.name == 'posix':
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , 'calc_distributions')
    else:
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , 'calc_distributions.exe')

    outfile = str(tmpfile_prefix) + '_dists.tsv'

    db_filename = Path.joinpath( Path(util.path_to_tcrdist_cpp_db) , 'tcrdist_info_{}.txt'.format( organism))

    cmd = '{} -f {} -m {} -d {} -a {} -b {} -o {}'\
    .format(exe, tcrs_file, max_dist, db_filename, achains_file, bchains_file, outfile)

    util.run_command(cmd, verbose=True)

    if not exists(outfile):
        print('tcr_clumping:: calc_distributions failed: missing', outfile)
        exit(1)

    counts = np.loadtxt(outfile, dtype=int)
    counts = np.cumsum(counts, axis=1)
    assert counts.shape == (len(tcrs), max_dist+1)
    n_bg_pairs = len(background_alpha_chains) * len(background_beta_chains)
    tcrdist_freqs = np.maximum(pseudocount, counts.astype(float))/n_bg_pairs

    for filename in [achains_file, bchains_file, tcrs_file, outfile]:
        os.remove(filename)

    return tcrdist_freqs


def assess_tcr_clumping(
        adata,
        outfile_prefix,
        radii = [24, 48, 72, 96],
        num_random_samples = 50000, # higher numbers are slower but allow more significant pvalues for extreme clumping
        pvalue_threshold = 1.0,
        verbose=True,
        also_find_clumps_within_gex_clusters=False,
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
        print('conga.tcr_clumping.assess_tcr_clumping:: need to compile the C++ tcrdist executables')
        exit(1)

    if also_find_clumps_within_gex_clusters:
        clusters_gex = np.array(adata.obs['clusters_gex'])

    organism = adata.uns['organism']
    num_clones = adata.shape[0]

    radii = [int(x+0.1) for x in radii] #ensure integers

    outprefix = outfile_prefix + '_tcr_clumping'

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    bg_freqs = estimate_background_tcrdist_distributions(
        adata.uns['organism'], tcrs, max(radii), num_random_samples=num_random_samples, tmpfile_prefix=outprefix)

    tcrs_file = outprefix +'_tcrs.tsv'
    adata.obs['va cdr3a vb cdr3b'.split()].to_csv(tcrs_file, sep='\t', index=False)


    # find neighbors in fg tcrs up to max(radii) #######################################

    if os.name == 'posix':
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors')
    else:
        exe = Path.joinpath( Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors.exe')

    agroups, bgroups = preprocess.setup_tcr_groups(adata)
    agroups_filename = outprefix+'_agroups.txt'
    bgroups_filename = outprefix+'_bgroups.txt'
    np.savetxt(agroups_filename, agroups, fmt='%d')
    np.savetxt(bgroups_filename, bgroups, fmt='%d')

    db_filename = Path.joinpath( Path(util.path_to_tcrdist_cpp_db), f'tcrdist_info_{organism}.txt')

    tcrdist_threshold = max(radii)

    cmd = '{} -f {} -t {} -d {} -o {} -a {} -b {}'\
    .format(exe, tcrs_file, tcrdist_threshold, db_filename, outprefix, agroups_filename, bgroups_filename)

    util.run_command(cmd, verbose=True)

    nbr_indices_filename = outprefix + '_nbr{}_indices.txt'.format( tcrdist_threshold)
    nbr_distances_filename = outprefix + '_nbr{}_distances.txt'.format( tcrdist_threshold)

    if not exists(nbr_indices_filename) or not exists(nbr_distances_filename):
        print('find_neighbors failed:', exists(nbr_indices_filename), exists(nbr_distances_filename))
        exit(1)

    all_nbrs = []
    all_distances = []
    for line1, line2 in zip(open(nbr_indices_filename,'r'), open(nbr_distances_filename,'r')):
        l1 = line1.split()
        l2 = line2.split()
        assert len(l1) == len(l2)
        #ii = len(all_nbrs)
        all_nbrs.append([int(x) for x in l1])
        all_distances.append([int(x) for x in l2])
    assert len(all_nbrs) == num_clones

    clone_sizes = adata.obs['clone_sizes']

    # use poisson to find nbrhoods with more tcrs than expected; have to handle agroups/bgroups
    dfl = []

    is_clumped = np.full((num_clones,), False)

    n_bg_pairs = num_random_samples * num_random_samples

    for ii in range(num_clones):
        ii_freqs = bg_freqs[ii]
        ii_dists = all_distances[ii]
        for radius in radii:
            num_nbrs = np.sum(x<=radius for x in ii_dists)
            max_nbrs = np.sum( (agroups!=agroups[ii]) & (bgroups!=bgroups[ii]))
            if num_nbrs:
                # adjust for number of tests
                mu = max_nbrs * ii_freqs[radius]
                pval = len(radii) * num_clones * poisson.sf( num_nbrs-1, mu )
                if pval< pvalue_threshold:
                    is_clumped[ii] = True
                    raw_count = ii_freqs[radius]*n_bg_pairs # if count was 0, will be pseudocount
                    if verbose:
                        print('tcr_nbrs_global: {:2d} {:9.6f} radius: {:2d} pval: {:9.1e} {:9.1f} tcr: {:3d} {} {}'\
                              .format( num_nbrs, mu, radius, pval, raw_count, clone_sizes[ii],
                                       ' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3])))
                    dfl.append( OrderedDict(clump_type='global',
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
                if also_find_clumps_within_gex_clusters:
                    ii_nbrs = all_nbrs[ii]
                    ii_cluster = clusters_gex[ii]
                    ii_cluster_mask = (clusters_gex==ii_cluster)
                    num_nbrs = np.sum( (x<=radius and clusters_gex[y]==ii_cluster) for x,y in zip(ii_dists, ii_nbrs))
                    if num_nbrs:
                        max_nbrs = np.sum( (agroups!=agroups[ii]) & (bgroups!=bgroups[ii]) & ii_cluster_mask)
                        mu = max_nbrs * ii_freqs[radius]
                        pval = len(radii) * num_clones * poisson.sf( num_nbrs-1, mu )
                        if pval< pvalue_threshold:
                            is_clumped[ii] = True
                            raw_count = ii_freqs[radius]*n_bg_pairs # if count was 0, will be pseudocount
                            if verbose:
                                print('tcr_nbrs_intra: {:2d} {:9.6f} radius: {:2d} pval: {:9.1e} {:9.1f} tcr: {:3d} {} {}'\
                                      .format( num_nbrs, mu, radius, pval, raw_count, clone_sizes[ii],
                                               ' '.join(tcrs[ii][0][:3]), ' '.join(tcrs[ii][1][:3])))
                            dfl.append( OrderedDict(clump_type='intra_gex_cluster',
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
        return results_df

    # identify groups of related hits?
    all_clumped_nbrs = {}
    for l in results_df.itertuples():
        ii = l.clone_index
        radius = l.nbr_radius
        clumped_nbrs = set(x for x,y in zip(all_nbrs[ii], all_distances[ii]) if y<= radius and is_clumped[x])
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
            new_nbr = min(nbr, np.min([all_smallest_nbr[x] for x in all_clumped_nbrs[ii]]))
            if nbr != new_nbr:
                all_smallest_nbr[ii] = new_nbr
                updated = True
        if not updated:
            break
    # define clusters, choose cluster centers
    clusters = np.array([0]*num_clones) # 0 if not clumped

    cluster_number=0
    for ii in clumped_inds:
        nbr = all_smallest_nbr[ii]
        if ii==nbr:
            cluster_number += 1
            members = [ x for x,y in all_smallest_nbr.items() if y==nbr]
            clusters[members] = cluster_number

    for ii, nbrs in all_clumped_nbrs.items():
        for nbr in nbrs:
            assert clusters[ii] == clusters[nbr] # confirm single-linkage clusters

    assert not np.any(clusters[is_clumped]==0)
    assert np.all(clusters[~is_clumped]==0)

    results_df['clumping_group'] = [ clusters[x.clone_index] for x in results_df.itertuples()]


    return results_df

