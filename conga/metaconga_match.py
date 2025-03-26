######################################################################################88

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sys import exit
from os.path import exists
import os
from pathlib import Path
from collections import Counter, namedtuple#, OrderedDict
#import time #debugging
#import random
#matplotlib.use('Agg')
#from scipy import stats
#from sklearn.metrics import pairwise_distances
#from sklearn.decomposition import PCA
#from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind, poisson
#import scipy.sparse as sps
#from scipy.sparse import issparse, csr_matrix
#from statsmodels.stats.multitest import multipletests

# package imports
from .tags import *
from . import preprocess
from . import tcr_scoring
from . import util
from . import correlations
from . import plotting
from . import tcr_clumping
#from .tcrdist.all_genes import all_genes
#import sys

# these are annoying
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

## setup some annotation dicts
all_aacluster_tags = {
    'cd4':{i:x for i,x in enumerate(
        '''btreg2 btreg5 ttreg3 btreg1 ttreg2
        thymic ttreg4 cyto5 cyto3 unk1
        btreg3 cyto2 btreg4 ttreg1 cyto1
        unk2 cyto4 small1 small2'''.split())},

    'cd8':{i:x for i,x in enumerate(
        '''helios2 temra3 temra2 helios3 temra1
        hobit2 hobhel unk1 cxcr5 maitnk
        nk2 helios1 nk1 thymic3 hobit1
        thymic1 thymic2 unk2 zen misc
        small'''.split())},
}


aacluster_groups = {
    'cd4': [
        [[0,1,3,10,12], 'blood-treg'],
        [[2,4,6,13], 'tissue-treg'],
        #[[0,1,2,3,4,6,10,12,13], 'treg'],
        [[7,8,11,14,16], 'cytotoxic'],
        [[9], 'unknown'],
        [[5], 'thymic'],
        #[[15,17], 'unk2_small1'],
    ], # drop 18
    'cd8': [
        [[0,3,5,6,11,14], 'hobit-helios'],
        [[1,2,4], 'temra'],
        [[10,12], 'nk12'],
        [[8], 'cxcr5'],
        [[13,15,16], 'thymic'],
        #[[9], 'mait'],
        #[[7], 'cd4'],
        #[[17,18,19], 'small'],
    ], # drop 20
}
aacluster_groups = {
    'cd4': [
        [[0,1,3,10,12], 'blood-treg'],
        [[2,4,6,13], 'tissue-treg'],
        #[[0,1,2,3,4,6,10,12,13], 'treg'],
        [[7,8,11,14,16], 'cytotoxic'],
        [[9], 'unknown'],
        [[5], 'Thymic'],
        #[[15,17], 'unk2_small1'],
    ], # drop 18
    'cd8': [
        [[0,3,5,6,11,14], 'HOBIT-HELIOS'],
        [[1,2,4], 'TEMRA'],
        [[10,12], 'dNKT'],
        [[8], 'CXCR5-IL10'],
        [[13,15,16], 'Thymic'],
        #[[9], 'mait'],
        #[[7], 'cd4'],
        #[[17,18,19], 'small'],
    ], # drop 20
}


aacluster2group = {}
for cd48, groups in aacluster_groups.items():
    for ii,(group,name) in enumerate(groups):
        for cluster in group:
            aacluster2group[(cd48, cluster)] = ii


######################################



path_to_mc_data = util.path_to_data / 'metaconga'

def get_categorical_colors(ncolors):
    cmap = plt.get_cmap('tab10') if ncolors <= 10 else plt.get_cmap('tab20')
    if ncolors <= 20:
        return cmap.colors[:ncolors]
    else:
        vals = np.linspace(0,1,ncolors-20)
        extra_colors = plt.get_cmap('viridis')(vals)[:,:3] # drop alpha
        return list(cmap.colors)+list(extra_colors)



def encode_tcr_seqs(
        tcr_df,
        min_cdr3len,
        max_cdr3len,
        va_genes,
        vb_genes,
        fg_trim=4,
        verbose=False,
):
    ''' This is for computing CDR3aa-bias-cluster enrichment scores

    cdr3a aas
    cdr3b aas
    cdr3a lens
    cdr3b lens
    va genes
    vb genes
    '''
    from .tcrdist.amino_acids import amino_acids


    num_lens = max_cdr3len-min_cdr3len+1
    num_vas = len(va_genes)
    num_vbs = len(vb_genes)

    num_cols = 2*20 + 2*num_lens + num_vas + num_vbs

    vecs = []
    va2ind = {x:i for i,x in enumerate(va_genes)}
    vb2ind = {x:i for i,x in enumerate(vb_genes)}

    for ii, l in enumerate(tcr_df.itertuples()):
        if verbose and ii%1000==0:
            print(ii, '%', 100*ii/tcr_df.shape[0])
        v_aa = []
        v_len = [0.]*(2*num_lens)
        for cdr3, lenshift in zip([l.cdr3a, l.cdr3b], [0, num_lens]):
            fg = cdr3[fg_trim:-fg_trim]
            norm = 1. if not fg else 1.0/len(fg)
            v_aa.extend([norm*fg.count(x) for x in amino_acids])

            ind = min(max_cdr3len, max(min_cdr3len, len(cdr3))) - min_cdr3len + lenshift
            v_len[ind] = 1.

        v_va = [0]*num_vas
        ind = va2ind.get(l.va, -1)
        if ind >= 0:
            v_va[ind] = 1.

        v_vb = [0]*num_vbs
        ind = vb2ind.get(l.vb, -1)
        if ind >= 0:
            v_vb[ind] = 1.

        vecs.append(v_aa + v_len + v_va + v_vb)
        assert len(vecs[-1]) == num_cols

    vecs = np.stack(vecs)
    assert vecs.shape == (tcr_df.shape[0], num_cols)

    return vecs



def get_cdr3aa_bias_tcr_scores(
        tcr_df,
        min_cdr3len=7,
        max_cdr3len=21,
):
    '''
    '''
    from .tcrdist.amino_acids import amino_acids

    # read the aacluster enrichments (V, cdr3len, cdr3aa)
    fname = Path.joinpath(
        path_to_mc_data, f'round8_v5cd4_run84_xribo_200_process_v1_leiden2_cdr3aa_'
        'sig_1e-06_cluster_feature_enrichments.tsv')
    enrichdf_cd4 = pd.read_table(fname)
    fname = Path.joinpath(
        path_to_mc_data, f'round8_v5cd8_run84_xribo_200_process_v1_leiden2_cdr3aa_'
        'sig_1e-06_cluster_feature_enrichments.tsv')
    enrichdf_cd8 = pd.read_table(fname)

    cdr3_cols = [f'cdr3{ab}_{aa}_log_enrich'
                 for ab in 'ab' for aa in amino_acids]

    len_cols = [f'cdr3len_{ab}{x}_log_enrich'
                for ab in 'ab' for x in range(min_cdr3len, max_cdr3len+1)]

    num_lens = max_cdr3len - min_cdr3len+1

    va_cols = sorted(set(
        [x for x in list(enrichdf_cd4.columns) + list(enrichdf_cd8.columns)
         if x.startswith('va') and x.endswith('_log_enrich')]))

    vb_cols = sorted(set(
        [x for x in list(enrichdf_cd4.columns) + list(enrichdf_cd8.columns)
         if x.startswith('vb') and x.endswith('_log_enrich')]))

    for col in va_cols+vb_cols:
        if col not in enrichdf_cd8.columns:
            enrichdf_cd8[col] = 0.0
        if col not in enrichdf_cd4.columns:
            enrichdf_cd4[col] = 0.0

    va_genes = [x[3:-11] for x in va_cols]
    vb_genes = [x[3:-11] for x in vb_cols]


    cols = cdr3_cols + len_cols + va_cols + vb_cols

    enrich_cd4 = enrichdf_cd4[cols].to_numpy().T
    enrich_cd8 = enrichdf_cd8[cols].to_numpy().T

    vecs = encode_tcr_seqs(tcr_df, min_cdr3len, max_cdr3len, va_genes, vb_genes)

    assert vecs.shape == (tcr_df.shape[0], enrich_cd8.shape[0])

    cd4_scores = vecs @ enrich_cd4
    cd8_scores = vecs @ enrich_cd8

    return cd4_scores, cd8_scores



def get_cdr3aa_bias_deg_scores(
        adata,
        max_pval_adj=0.05,
        max_degs=50,
        skip_ribosomal=True,
        topn_for_good_score=5, # ie, the 3rd highest score
):
    '''
    '''

    dfl = []

    for cd48, degs_runtag in zip(['cd4','cd8'], ['run105','run106']):
        degsfile = Path.joinpath(
            path_to_mc_data, f'{degs_runtag}_{cd48}_deg_results.tsv')

        dfdegs = pd.read_table(degsfile)

        dfdegs = dfdegs[dfdegs.pval_adj <= max_pval_adj].copy()

        if skip_ribosomal:
            dfdegs = dfdegs[~dfdegs.is_ribo].copy()

        num_clusters = dfdegs.leiden.max()+1

        adata_genes = set(adata.raw.var_names)
        badmask = ~dfdegs.gene.isin(adata_genes)
        if badmask.sum():
            print(f'dropping degs not present in adata.raw: {badmask.sum()}',
                  dfdegs.shape[0])
            dfdegs = dfdegs[~badmask].copy()

        all_scores = []
        for c in range(num_clusters):
            df = dfdegs[dfdegs.leiden==c].sort_values('neg_abs_score').head(max_degs)
            good_score = df.neg_abs_score.head(topn_for_good_score).median()
            #print('good_score:', c, good_score)

            df['weight'] = df.neg_abs_score * df.log2fold/ good_score

            # print(c)
            # print(df['gene weight'.split()].head())

            genes = list(df.gene)

            inds = [adata.raw.var_names.get_loc(x) for x in genes]

            X = adata.raw.X[:,inds].toarray()

            # subtract mean, std
            Xstd = np.maximum(1e-3, X.std(axis=0, keepdims=True)) # no nans!
            X = (X - X.mean(axis=0, keepdims=True))/Xstd
            clip = 10
            X = np.minimum(clip, np.maximum(-clip, X))

            weights = np.array(df.weight)
            assert X.shape == (adata.shape[0], weights.shape[0])

            scores = (X * weights[None,:]).sum(axis=1)

            all_scores.append(scores[:,None])

        if cd48=='cd4':
            cd4_scores = np.hstack(all_scores)
        else:
            cd8_scores = np.hstack(all_scores)

    return cd4_scores, cd8_scores



## hack: copied this little block from devel.py since I don't want to import devel here
##
_tfs_listfile= util.path_to_data / 'Lambert_et_al_PMID_29425488_TF_names_v_1.01.txt'
human_tf_gene_names = [x.split()[0] for x in open(_tfs_listfile,'r')]


def find_degs_for_subset(
        inds,
        adata,
        pval_threshold=0.05,
):


    inds = set(inds)
    assert 0 <= min(inds) and max(inds) < adata.shape[0]

    # look at differentially expressed genes in the subset vs the rest
    obs_tag = 'tmp_groups'
    adata.obs[obs_tag] = [str(int(x in inds)) for x in range(adata.shape[0])]
    key_added = 'degs_for_subset'
    rank_method = 'wilcoxon'
    group = '1'
    print('run rank genes:', len(inds), 'vs', adata.shape[0]-len(inds))

    # not sure why this is happening, maybe because the adata we downloaded from
    #  the original study had some messed up metadata
    #
    if 'log1p' in adata.uns_keys() and 'base' not in adata.uns['log1p']:
        #print('log1p:', adata.uns['log1p'])
        adata.uns['log1p']['base'] = None

    sc.tl.rank_genes_groups(
        adata,
        groupby=obs_tag,
        method=rank_method,
        groups=[group],
        reference='rest',
        key_added=key_added,
    )

    dfl = []
    for igene,gene in enumerate( adata.uns[key_added]['names'][group] ):
        log2fold = adata.uns[key_added]['logfoldchanges'][group][igene]
        pval_adj = adata.uns[key_added]['pvals_adj'][group][igene]
        score    = adata.uns[key_added]['scores'][group][igene]
        if pval_adj < pval_threshold:
            dfl.append(dict(
                rank=igene,
                score=score,
                gene=gene,
                pval_adj=pval_adj,
                log2fold=log2fold,
                is_tf=gene in human_tf_gene_names,
                is_ribo=gene[:3] in ['RPL','RPS'] and gene[3].isdigit(),
            ))

    results = pd.DataFrame(dfl)
    return results


Matchinfo = namedtuple('Matchinfo', ['pvals','degs','obs'])

def find_aacluster_matches(
        adata_in,
        cd48,
        NO_NBRS = False,
        subset_fracs = [0.01, 0.025, 0.05, 0.1],
        pval_threshold = 1e-3, # originally 1e-2
        pval_threshold_for_rank_genes = 1e-6,
        fdr_threshold=0.25,
        max_gex_cluster_mait_frac = 0.33,
):
    ''' returns  Matchinfo namedtuple

    assumes that adata has already been set up, for example:
    * X_gex_2d in adata.obsm
    * clusters_gex in adata.obs

    will recalculate GEX neighbors since we want a particular nbr_frac
    '''
    from scipy.stats import hypergeom, norm, ttest_ind, mannwhitneyu, linregress
    from statsmodels.stats.multitest import multipletests
    from .preprocess import add_mait_info_to_adata_obs

    assert cd48 in ['cd4','cd8']

    adata = adata_in.copy() # safer, since we might modify??

    obs = adata.obs.copy() # fill this with new info
    obs['original_index'] = np.arange(obs.shape[0])
    obs['cd48_for_aacluster_matching'] = cd48
    
    # store gex umap coords for plotting
    obs['X_gex_2d_0'] = adata.obsm['X_gex_2d'][:,0]
    obs['X_gex_2d_1'] = adata.obsm['X_gex_2d'][:,1]


    ## drop MAIT/iNKT clones, both by GEX cluster and by clone TCR sequence
    if 'is_invariant' not in adata.obs.columns:
        add_mait_info_to_adata_obs(adata)

    mait_fracs = adata.obs.groupby(
        'clusters_gex').is_invariant.mean().sort_values(ascending=False)

    mait_clusters = [x for x,y in mait_fracs.items() if y>max_gex_cluster_mait_frac]
    drop_mask = np.array(adata.obs.is_invariant)

    if mait_clusters:
        drop_mask |= np.array(adata.obs.clusters_gex.isin(mait_clusters))
        print('mait_clusters:', mait_clusters)
        print('mait_fracs:', mait_fracs, sep='\n')

    if drop_mask.sum():
        print('drop MAIT/iNKT by cluster and by clone, num_to_drop:', drop_mask.sum())
        adata = adata[~drop_mask].copy()
        obs = obs[~drop_mask].copy()


    all_gex_scores_cd4, all_gex_scores_cd8 = get_cdr3aa_bias_deg_scores(adata)
    all_tcr_scores_cd4, all_tcr_scores_cd8 = get_cdr3aa_bias_tcr_scores(adata.obs)

    if cd48 == 'cd4':
        all_gex_scores = all_gex_scores_cd4
        all_tcr_scores = all_tcr_scores_cd4
    else:
        all_gex_scores = all_gex_scores_cd8
        all_tcr_scores = all_tcr_scores_cd8

    num_clusters = min(all_gex_scores.shape[1], all_tcr_scores.shape[1])
    gex_cluster_counts = Counter(adata.obs.clusters_gex)

    def get_overlap_cluster_fracs(inds):
        counts = Counter(adata.obs.iloc[list(inds)].clusters_gex).most_common(2)
        s = ','.join(f'{x}:{int(np.round(100*y/len(inds)))}'
                     f':{int(np.round(100*y/gex_cluster_counts[x]))}'
                     for x,y in counts)
        return s

    if not NO_NBRS: # nbr computation
        nbr_frac = 0.025 if adata.shape[0]<30000 else 0.01

        all_nbrs = preprocess.calc_nbrs(
            adata, [nbr_frac], obsm_tag_gex='X_pca_gex', obsm_tag_tcr = None,
        )

        nbrs = all_nbrs[nbr_frac][0] # gex nbrs
        num_nbrs = nbrs.shape[1]


    dfl = []
    for c in range(num_clusters):
        best_pval = 1.0 # for storing the inds, etc
        best_log2enrich = 0

        ctag = all_aacluster_tags[cd48][c]
        gex_scores, raw_tcr_scores = all_gex_scores[:,c], all_tcr_scores[:,c]
        obs[f'gex_bias_score_{c}'] = gex_scores
        obs[f'tcr_bias_score_{c}'] = raw_tcr_scores

        inds_gex = np.argsort(gex_scores)[::-1]
        N = adata.shape[0]

        for subset_frac in subset_fracs:
            topn = int(np.round(subset_frac*N))
            if topn<5:
                continue
            mask_gex = np.zeros((N,), dtype=bool)
            mask_gex[inds_gex[:topn]] = True

            for tcr_topn_type in ['raw','fdr','avg']:
                if NO_NBRS and tcr_topn_type != 'raw':
                    continue

                # normalize
                tcr_scores = ((raw_tcr_scores-np.mean(raw_tcr_scores))/
                              np.std(raw_tcr_scores))

                mask_tcr = np.zeros((N,), dtype=bool)
                if tcr_topn_type == 'raw':
                    inds_tcr = np.argsort(tcr_scores)[::-1]
                    mask_tcr[inds_tcr[:topn]] = True
                else:
                    # nbr-avg
                    tcr_scores = ((np.sum(tcr_scores[nbrs], axis=1)+tcr_scores)/
                                  (num_nbrs+1))
                    tcr_scores *= np.sqrt(num_nbrs+1)
                    obs[f'avg_tcr_bias_score_{c}'] = tcr_scores

                    if tcr_topn_type == 'avg':
                        inds_tcr = np.argsort(tcr_scores)[::-1]
                        mask_tcr[inds_tcr[:topn]] = True
                    else:
                        assert tcr_topn_type == 'fdr'
                        res = multipletests(
                            norm.sf(tcr_scores), fdr_threshold, method='fdr_bh')
                        reject, pvals = res[:2]
                        good_inds_tcr = np.nonzero(reject)[0]
                        mask_tcr[good_inds_tcr] = True

                combo_scores = (np.argsort(np.argsort(gex_scores)) + # bigger is
                                np.argsort(np.argsort(tcr_scores)))/(2*N)# better

                if tcr_topn_type != 'fdr':
                    obs[f'{tcr_topn_type}_combo_bias_score_{c}'] = combo_scores

                overlap = (mask_tcr&mask_gex).sum()
                if not overlap:
                    continue
                topn_gex, topn_tcr = mask_gex.sum(), mask_tcr.sum()
                assert topn_gex == topn and (topn_tcr==topn or tcr_topn_type=='fdr')
                pval = hypergeom.sf(overlap-1, N, topn_gex, topn_tcr)
                expect = topn_gex*topn_tcr/adata.shape[0]
                real_overlap = overlap-expect
                log2enrich = np.log2(overlap/expect)

                mait_fraction = adata.obs.is_invariant[mask_tcr&mask_gex].mean()
                overlap_inds = np.nonzero(mask_tcr&mask_gex)[0]
                overlap_scores = combo_scores[overlap_inds]
                sorted_overlap_inds = overlap_inds[np.argsort(overlap_scores)[::-1]]
                real_overlap_inds= sorted_overlap_inds[:int(np.round(real_overlap))]
                if pval <= pval_threshold:
                    pval_type = f'{tcr_topn_type}_overlap_{subset_frac}'
                    dfl.append(dict(
                        cd48=cd48,
                        aacluster=c,
                        ctag=ctag,
                        pval=pval,
                        pval_type = pval_type,
                        overlap=overlap, log2enrich=log2enrich,
                        real_overlap=real_overlap,
                        overlap_clusters=get_overlap_cluster_fracs(overlap_inds),
                        real_overlap_clusters=get_overlap_cluster_fracs(
                            real_overlap_inds),
                        topn_gex=topn_gex, topn_tcr=topn_tcr,
                        subset_frac=subset_frac,
                        num_clonotypes=N, mait_fraction=mait_fraction,
                    ))

                    if (pval < best_pval or
                        pval == best_pval and log2enrich > best_log2enrich):
                        best_pval = pval
                        best_log2enrich = log2enrich
                        mask = np.zeros((N,), dtype=bool)
                        mask[real_overlap_inds] = True
                        obs[f'real_overlap_mask_{c}'] = mask
                        obs[f'real_overlap_mask_{c}_pval_type'] = pval_type
                        obs[f'real_overlap_mask_{c}_pval'] = pval


                    print(f'{cd48} {c:2d} {ctag:7s} {subset_frac:.3f} '
                          f'{tcr_topn_type} {pval:9.2e} {overlap:3d} '
                          f' {real_overlap:5.1f} {log2enrich:5.2f} '
                          f'mf {mait_fraction:.2f}',
                          get_overlap_cluster_fracs(overlap_inds),
                          get_overlap_cluster_fracs(real_overlap_inds))

    # maybe some deg calcs
    dfl_degs = []
    for c in range(num_clusters):
        ctag = all_aacluster_tags[cd48][c]
        col = f'real_overlap_mask_{c}_pval'
        overlap_pval = obs[col].iloc[0] if col in obs.columns else 10.
        if overlap_pval <= pval_threshold_for_rank_genes:
            mask = np.array(obs[f'real_overlap_mask_{c}'])
            inds = np.nonzero(mask)[0]

            dfdegs = find_degs_for_subset(inds, adata)
            if dfdegs.shape[0] > 0:
                overlap_pval_type = obs[col+'_type'].iloc[0]
                dfdegs['neg_abs_score'] = -1 * np.abs(dfdegs.score)
                dfdegs.sort_values('neg_abs_score', inplace=True)
                dfdegs['rank'] = np.arange(dfdegs.shape[0])
                dfdegs['aacluster'] = c
                dfdegs['ctag'] = ctag
                dfdegs['cd48'] = cd48
                dfdegs['overlap_size'] = len(inds)
                dfdegs['overlap_pval'] = overlap_pval
                dfdegs['overlap_pval_type'] = overlap_pval_type

                dfl_degs.append(dfdegs)

                # drop ribo and under expressed
                dfshow = dfdegs[(~dfdegs.is_ribo) & (dfdegs.log2fold>0)].copy()
                print('degs:', c, ctag, cd48, overlap_pval,
                      overlap_pval_type)
                cols = 'rank gene score pval_adj log2fold'.split()
                print(dfshow[cols].head(10))

    dfpvals = pd.DataFrame(dfl)
    dfdegs = pd.concat(dfl_degs) if dfl_degs else pd.DataFrame()


    return Matchinfo(pvals=dfpvals, degs=dfdegs, obs=obs)



def reduce_to_single_aacluster_match_per_clonotype(
        matches,
        pval_threshold=1e-3,
        DROP_THYMIC_HITS=True,
        DROP_SMALL_HITS=True,
):
    ''' return hits DataFrame which has one row for each matched clonotype,
    otherwise kind of like matches.obs
    '''
    MAX_CLUSTERS = 25 # just big
    assert matches.obs.cd48_for_aacluster_matching.nunique() == 1
    cd48 = matches.obs.cd48_for_aacluster_matching.iloc[0]
    
    ## now reduce matches to best hit per clonotype

    drop_clusters = {'cd8':[], 'cd4':[]}

    if DROP_SMALL_HITS:
        drop_clusters['cd8'].extend([17, 18, 19, 20])
        drop_clusters['cd4'].extend([17, 18])

    if DROP_THYMIC_HITS:
        drop_clusters['cd8'].extend([13,15,16])
        drop_clusters['cd4'].extend([5])


    # assuming that we already dropped invariants at the matching stage

    # sorted...
    cols = ['real_overlap_mask_'+str(x) for x in range(MAX_CLUSTERS)
            if 'real_overlap_mask_'+str(x) in matches.obs.columns and
            matches.obs[f'real_overlap_mask_{x}_pval'].iloc[0] <= pval_threshold]
    if not cols:
        print('no hits:', cd48)
        return pd.DataFrame() # empty dataframe

    all_masks = []
    all_scores = []
    clusters = []

    for col in cols:
        cluster = int(col.split('_')[-1])
        if cluster in drop_clusters[cd48]:
            continue
        mask = np.array(matches.obs[col])
        scores = np.array(matches.obs[f'avg_combo_bias_score_{cluster}'])
        all_masks.append(mask[:,None])
        all_scores.append(scores[:,None])
        clusters.append(cluster)

    if not clusters:
        print('no undropped clusters:', cd48)
        return pd.DataFrame()

    all_masks = np.hstack(all_masks)
    all_scores = np.hstack(all_scores)
    all_scores *= all_masks
    any_mask = all_masks.sum(axis=1)>0.5
    top_cluster = np.argmax(all_scores, axis=1)
    print(f'{cd48} any_frac: {any_mask.sum()/matches.obs.shape[0]:.3f}')


    dfl = []
    for ic, cluster in enumerate(clusters):
        mask = all_masks[:,ic] & (top_cluster==ic)
        res = matches.obs[mask].copy()
        res['cd48'] = cd48
        res['aacluster'] = cluster
        res['aacluster_tag'] = all_aacluster_tags[cd48][cluster]
        dropcols = [x for x in res.columns if 'bias_score' in x or
                    'overlap_mask' in x] # to drop later
        cols = [x for x in res.columns
                if 'bias_score' in x and x.endswith(f'_{cluster}')]
        assert len(cols) == 5
        for col in cols:
            newcol = '_'.join(col.split('_')[:-1])
            res[newcol] = res[col]
        for suf in ['pval_type','pval']:
            col = f'real_overlap_mask_{cluster}_{suf}'
            newcol = f'real_overlap_mask_{suf}'
            res[newcol] = res[col]

        res.drop(columns=dropcols, inplace=True)
        #print(cluster, res.shape[0])
        dfl.append(res)


    hits = pd.concat(dfl)
    return hits


                
def plot_aacluster_matches(
        adata, # so we can stash the pngfiles in adata.uns['conga_results']
        matches,
        outfile_prefix,
        pval_threshold_for_bars = 1e-3,
        pval_threshold_for_overlap_umaps = 1e-6,
        MAX_UMAP_ROWS=10,
        PLOT_GEX_CLUSTERS_UMAP = False, # set True to show umaps colored by GEX clusters
        gexcluster_degs = None, # pass to also show degs from the closest GEX cluster
        DROP_SMALL_HITS = True,
        DROP_THYMIC_HITS = True,
):
    ''' matches is the namedtuple returned by find_aacluster_matches
    has (pvals, degs, obs)
    
    '''
    assert matches.pvals.cd48.nunique() == 1
    cd48 = matches.pvals.cd48.iloc[0]

    obs = matches.obs # convenience

    hits = reduce_to_single_aacluster_match_per_clonotype(
        matches,
        pval_threshold=pval_threshold_for_bars,
        DROP_SMALL_HITS=DROP_SMALL_HITS,
        DROP_THYMIC_HITS=DROP_THYMIC_HITS,
    )

    # save some TSV files
    for df, tag in zip([matches.pvals, matches.degs, matches.obs, hits],
                       ['pvals','degs','obs','hits']):
        outfile = f'{outfile_prefix}_aacluster_match_{tag}.tsv'
        df.to_csv(outfile, sep='\t', index=False)
        print('made:', outfile)

    ###############################
    ## make bar plots of matches, relative to total clones, total cells

    catcolors = get_categorical_colors(21) # should be plenty for cd4 or cd8

    num_cgroups = max(y+1 for x,y in aacluster2group.items() if x[0] == cd48)
    catcolors_cgroups = get_categorical_colors(num_cgroups)
    cgroup_names = [x[1] for x in aacluster_groups[cd48]]

    nrows, ncols, plotno = 2, 3, 0
    plt.figure(figsize=(15, 7.5))

    overlap_pval_map = hits.drop_duplicates('aacluster')\
                           .set_index('aacluster').real_overlap_mask_pval
    
    for USE_CGROUPS in [False, True]:
        for BY_CELLS in [False, True]:


            bars = []
            texts = []
            plotno += 1
            plt.subplot(nrows, ncols, plotno)

            # the keys for counts are aacluster numbers
            if BY_CELLS:
                frac_denom = obs.clone_sizes.sum()
                counts = Counter()
                for l in hits.itertuples():
                    counts[l.aacluster] += l.clone_sizes
            else:
                frac_denom = obs.shape[0] #num_clonotypes
                counts = Counter(hits.aacluster)

            if USE_CGROUPS:
                new_counts = Counter()
                for c, count in counts.items():
                    if (cd48,c) in aacluster2group:
                        new_counts[aacluster2group[(cd48,c)]] += count

                total = 0.
                for cg in range(num_cgroups):
                    frac = new_counts[cg] / frac_denom
                    bars.append( (0, frac, total, catcolors_cgroups[cg]))
                    total += frac
            else:
                total = 0.
                for c, count in counts.items():
                    frac = count / frac_denom
                    bars.append( (0, frac, total, catcolors[c]))
                    texts.append((0.5, total+frac/2, catcolors[c],
                                  f'overlap_pval: {overlap_pval_map[c]:9.2e}'))
                    total += frac

            plt.bar(
                x=[x[0] for x in bars],
                height=[x[1] for x in bars],
                bottom=[x[2] for x in bars],
                color=[x[3] for x in bars],
                width=0.8,
            )
            for x,y,color,text in texts:
                plt.text(x,y,text, ha='left', va='center', fontsize=9, color=color)
            plt.xlim((-.5,1.5))
            ymx = max(x[1]+x[2] for x in bars)
            plt.ylim((0, 1.55*ymx))
            if BY_CELLS:
                plt.ylabel('fraction of cells')
            else:
                plt.ylabel('fraction of clonotypes')

            if USE_CGROUPS:
                for ii, tag in enumerate(cgroup_names):
                    if tag.lower() == 'thymic' and DROP_THYMIC_HITS:
                        continue
                    plt.bar(x=0, height=0, bottom=0, color=catcolors_cgroups[ii],
                            label=tag)
                plt.legend(ncols=2, loc='upper right', title='AAcluster group',
                           fontsize=7)
            else:
                for cluster in counts.keys():
                    plt.bar(x=0, height=0, bottom=0, color=catcolors[cluster],
                            label=all_aacluster_tags[cd48][cluster])
                plt.legend(ncols=2, loc='upper right', title='AAcluster',
                           fontsize=7)
            plt.xticks([],[])

            by = 'cell' if BY_CELLS else 'clonotype'
            plt.title(f'{cd48.upper()} AAcluster matches by {by}')


        ## add a umap plot
        plotno += 1
        plt.subplot(nrows, ncols, plotno)

        plt.scatter(obs.X_gex_2d_0, obs.X_gex_2d_1, c='#EEEEEE', s=5)
        plt.xlabel('GEX UMAP_1')
        plt.ylabel('GEX UMAP_2')
        plt.xticks([],[])
        plt.yticks([],[])

        if USE_CGROUPS:
            seen = set()
            for c, count in hits.aacluster.value_counts().items():
                if (cd48,c) in aacluster2group:
                    grp = aacluster2group[(cd48,c)]
                    color = catcolors_cgroups[grp]
                    res = hits[hits.aacluster==c]
                    if grp in seen:
                        label = None
                    else:
                        label = cgroup_names[grp]
                        seen.add(grp)
                    plt.scatter(res.X_gex_2d_0, res.X_gex_2d_1, color=color, s=5,
                                label=label)
            plt.legend(fontsize=7)
        else:
            total = 0.
            for c, count in hits.aacluster.value_counts().items():
                res = hits[hits.aacluster==c]
                plt.scatter(res.X_gex_2d_0, res.X_gex_2d_1, color=catcolors[c],
                            s=5, label=all_aacluster_tags[cd48][c])
            plt.legend(fontsize=7)

            

    plt.tight_layout()
    pngfile = outfile_prefix+'_aacluster_match_bars.png'
    plt.savefig(pngfile)
    adata.uns.setdefault('conga_results',{})[AACLUSTER_MATCH_BARS] = pngfile
    print('made:', pngfile)


    ################################################################
    ## now a plot showing the UMAP overlaps
    #pval_threshold = 1e-3
    degs_pval_threshold = 0.05

    dfpvals = matches.pvals[matches.pvals.pval<=pval_threshold_for_overlap_umaps].copy()

    dfpvals['sortscore'] = [x.pval if x.pval>0 else -x.log2enrich
                            for x in dfpvals.itertuples()]
    dfpvals = dfpvals.sort_values('sortscore').drop_duplicates('aacluster')

    nrows = min(MAX_UMAP_ROWS, dfpvals.shape[0])
    ncols = 4 + PLOT_GEX_CLUSTERS_UMAP + (gexcluster_degs is not None)
    plotno = 0

    plt.figure(figsize=(ncols*3, nrows*3))

    for row in dfpvals.reset_index().head(nrows).itertuples():
        last_row = (row.Index+1 == nrows)
        title_fontsize = 10
        cluster = row.aacluster
        ctag = all_aacluster_tags[cd48][cluster]
        cctag = f'C{cluster}-{ctag}'

        # read the obs file
        tcr_scores = np.array(obs[f'avg_tcr_bias_score_{cluster}'])
        gex_scores = np.array(obs[f'gex_bias_score_{cluster}'])
        combo_scores = np.array(obs[f'avg_combo_bias_score_{cluster}'])
        overlap_mask = np.array(obs[f'real_overlap_mask_{cluster}'])
        assert overlap_mask.sum() == int(np.round(row.real_overlap))

        # find greatest overlap with better-scoring mask
        cols = [x for x in obs.columns if x.startswith('real_overlap_mask_') and
                x.endswith('_pval')]
        best_jaccard = 0
        for col2 in cols:
            c2 = int(col2.split('_')[-2])
            if c2==cluster:
                continue
            pval2 = obs[col2].iloc[0]
            if pval2 < row.pval:
                mask2 = np.array(obs[col2[:-5]])
                jaccard = (overlap_mask&mask2).sum() / (overlap_mask|mask2).sum()
                if jaccard>best_jaccard:
                    best_jaccard = jaccard
                    best_jaccard_cctag = f'C{c2}-{all_aacluster_tags[cd48][c2]}'
        if best_jaccard:
            print('best_jaccard:', best_jaccard, best_jaccard_cctag)


        # which is the closest gex cluster??
        num_gex_clusters = obs.clusters_gex.max()+1
        best_gc_jaccard = 0
        for gc in range(num_gex_clusters):
            gcmask = np.array(obs.clusters_gex==gc)
            jaccard = (overlap_mask&gcmask).sum()/(overlap_mask|gcmask).sum()
            if jaccard>best_gc_jaccard:
                best_gc_jaccard = jaccard
                best_gc = gc
        #print('best_gc_jaccard:', best_gc_jaccard)
        print(row.pval, row.log2enrich, cd48, cctag, best_jaccard, best_gc_jaccard)

        ##### now start plotting.........
        xy = np.array([[x,y] for x,y in zip(obs.X_gex_2d_0, obs.X_gex_2d_1)])

        if PLOT_GEX_CLUSTERS_UMAP:
            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            num_clusters = np.max(obs.clusters_gex)+1
            catcolors = get_categorical_colors(num_clusters)
            obs['color'] = [catcolors[x] for x in obs.clusters_gex]

            plt.scatter(xy[:,0], xy[:,1], c=obs.color, s=5)
            plt.title('GEX clusters', fontsize=title_fontsize)
            plt.xticks([],[]) ; plt.yticks([],[])
            plt.ylabel('GEX UMAP_2')
            if last_row:
                plt.xlabel('GEX UMAP_1')

        plotno += 1
        plt.subplot(nrows, ncols, plotno)

        mx = np.max(np.abs(tcr_scores))
        reorder = np.argsort(np.abs(tcr_scores))
        plt.scatter(xy[reorder,0], xy[reorder,1], c=tcr_scores[reorder],
                    cmap='coolwarm', vmin=-mx, vmax=mx, s=5)
        plt.xticks([],[]) ; plt.yticks([],[])
        if last_row:
            plt.xlabel('GEX UMAP_1')
        if not PLOT_GEX_CLUSTERS_UMAP:
            plt.ylabel('GEX UMAP_2')
        plt.title(f'{cctag} TCR scores', fontsize=title_fontsize)
        #plt.title(f'{cctag} TCR scores, minmax= {mx:.2f}', fontsize=title_fontsize)

        plotno += 1
        plt.subplot(nrows, ncols, plotno)

        mx = np.max(np.abs(gex_scores))
        reorder = np.argsort(np.abs(gex_scores))
        plt.scatter(xy[reorder,0], xy[reorder,1], c=gex_scores[reorder],
                    cmap='coolwarm', vmin=-mx, vmax=mx, s=5)
        plt.xticks([],[]) ; plt.yticks([],[])
        if last_row:
            plt.xlabel('GEX UMAP_1')
        plt.title(f'{cctag} GEX scores', fontsize=title_fontsize)
        #plt.title(f'{cctag} GEX scores, minmax= {mx:.2f}', fontsize=title_fontsize)

        plotno += 1
        plt.subplot(nrows, ncols, plotno)

        mx = np.max(np.abs(gex_scores))
        plt.scatter(obs.X_gex_2d_0, obs.X_gex_2d_1, c='#DDDDDD', s=5)

        colors = combo_scores[overlap_mask]
        reorder = np.argsort(colors)
        plt.scatter(xy[overlap_mask,0][reorder], xy[overlap_mask,1][reorder],
                    c=colors[reorder], s=5)
        title = (f'overlap: {overlap_mask.sum()} of {obs.shape[0]}, '
                  f'pval= {row.pval:9.2e}')
        if best_jaccard>.3:
            msg = (f'WARNING: Significant\noverlap w/ {best_jaccard_cctag}\n'
                   f'jaccard= {best_jaccard:.3f}')
            plt.text(0.01,0.01,msg, ha='left', va='bottom', fontsize=7,
                     transform=plt.gca().transAxes)
        #if best_jaccard>.3:
        #    title += f'\nbetter pval jaccard: {best_jaccard:.3f} {best_jaccard_cctag}'
        if last_row:
            plt.xlabel('GEX UMAP_1')
        plt.title(title, fontsize=title_fontsize)
        plt.xticks([],[]) ; plt.yticks([],[])


        topn = 25
        max_deg_pval = 0.05
        plotno += 1
        dfdegs = matches.degs
        if dfdegs.shape[0]:

            dfd = dfdegs[(dfdegs.pval_adj <= degs_pval_threshold) &
                         (dfdegs.aacluster==cluster) &
                         (~dfdegs.is_ribo) &
                         (dfdegs.score>0)].sort_values(
                             'score', ascending=False).head(topn).reset_index()

            assert (dfd.overlap_size==overlap_mask.sum()).all()

            if dfd.shape[0]:
                plt.subplot(nrows, ncols, plotno)

                plt.scatter(dfd.index, dfd.score, c='white', zorder=1)
                for l in dfd.itertuples():
                    plt.text(l.Index, l.score, l.gene, rotation='vertical',
                             va='bottom', ha='center', fontsize=7, zorder=2)

                plt.title(f'{cctag} DEGs', fontsize=title_fontsize)
                plt.xticks([],[])
                ymn,ymx = plt.ylim()
                plt.ylim((ymn,ymn+(ymx-ymn)*1.2))
                plt.xlabel('DEG rank')
                
        if gexcluster_degs is not None:
            topn = 25
            # read gex cluster degs
            gexcluster_degs['is_ribo'] = gexcluster_degs.gene.str.contains(
                '^RP[SL]', regex=True) # sometimes we require [0-9] after the [SL]!
            dfd = gexcluster_degs[
                ( gexcluster_degs.leiden == best_gc) &
                ( gexcluster_degs.pval_adj <= degs_pval_threshold) &
                (~gexcluster_degs.is_ribo) &
                ( gexcluster_degs.score > 0)].sort_values(
                    'score', ascending=False).head(topn).reset_index()

            plotno += 1
            plt.subplot(nrows, ncols, plotno)

            plt.scatter(dfd.index, dfd.score, c='white', zorder=1)
            for l in dfd.itertuples():
                plt.text(l.Index, l.score, l.gene, rotation='vertical',
                         va='bottom', ha='center', fontsize=7, zorder=2,
                         color=catcolors[best_gc])

            plt.title(
                f'{cctag} gexleiden {best_gc} DEGs (jaccard={best_gc_jaccard:.2f})',
                color=catcolors[best_gc], fontsize=title_fontsize)
            plt.xticks([],[])

    if nrows:
        plt.tight_layout()
        pngfile = f'{outfile_prefix}_aacluster_match_umaps.png'
        adata.uns.setdefault('conga_results',{})[AACLUSTER_MATCH_UMAPS] = pngfile
        plt.savefig(pngfile, dpi=150)
        print('made:', pngfile)


