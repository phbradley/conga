################################################################################
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import sys
from sys import exit
import pandas as pd
from pathlib import Path
import math
from sklearn.decomposition import PCA
from . import svg_basic
from . import util
from . import preprocess
from . import pmhc_scoring
from . import correlations
from . import tcr_scoring
from .tags import *
from .tcrdist.make_tcr_logo import make_tcr_logo_for_tcrs
from .tcrdist.tcr_distances import TcrDistCalculator
from .tcrdist.util import assign_colors_to_conga_tcrs
from .tcrdist.all_genes import all_genes
from .tcrdist.make_tcr_trees import make_tcr_tree_svg_commands
import os
from os.path import exists
from collections import Counter, OrderedDict
import matplotlib.image as mpimg
from matplotlib.collections import LineCollection
from scipy.cluster import hierarchy
from scipy.spatial import distance
import time
import random


default_logo_genes = {
    'human': ['CD4','CD8A','CD8B','CCR7','SELL',
              'GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1',
              'CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
    'mouse': ['Cd4', 'Cd8a', 'Cd8b1', 'Ccr7', 'Sell',
              'Itgal', 'Prf1', 'Gzma', 'Il2rb', 'Gzmk', 'Ifng',
              'Ccl5', 'Cxcr3', 'Zbtb16', 'Nkg7', 'Klrd1'],
    # should probably specialize these
    'human_gd': ['CD4','CD8A','CD8B','CCR7','SELL',
                 'GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1',
                 'CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
    'mouse_gd': ['Cd4', 'Cd8a', 'Cd8b1', 'Ccr7', 'Sell',
                 'Itgal', 'Prf1', 'Gzma', 'Il2rb', 'Gzmk', 'Ifng',
                 'Ccl5', 'Cxcr3', 'Zbtb16', 'Nkg7', 'Klrd1'],
    # b cells
    'human_ig': ['IL4R','TCL1A','SELL','CRIP1','CD27',
                 'ZFP36','HLA-C','HLA-DRB1','COTL1','JCHAIN','XBP1',
                 'IGHD','IGHM','IGHG1','IGHA1','IGHA2']
}


default_gex_header_genes = {
    'human': ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
    'mouse': ['clone_sizes','Cd4', 'Cd8a', 'Cd8b1', 'Sell', 'Itgal', 'Gzma', 'Ccl5', 'Il2rb', 'Ikzf2', 'Pdcd1', 'Zbtb16'],
    'human_gd': ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
    'mouse_gd': ['clone_sizes','Cd4', 'Cd8a', 'Cd8b1', 'Sell', 'Itgal', 'Gzma', 'Ccl5', 'Il2rb', 'Ikzf2', 'Pdcd1', 'Zbtb16'],
    'human_ig': ['clone_sizes','IL4R','TCL1A','SELL','IGKC','HLA-DRB1','CD27','JCHAIN','XBP1','COTL1','EGR1','IGHM'],
}


# this is a bit of a hack
DONT_SHOW_LOW_MAIT_SCORES = True

def get_integers_color_dict( num_categories, cmap=None):
    C = {}
    if cmap is None:
        if num_categories <= 10:
            cmap = plt.get_cmap('tab10')
        else:
            cmap = plt.get_cmap('tab20')
        for i in range(num_categories):
            C[i] = cmap.colors[i%len(cmap.colors)]
    else:
        # this is not the best, eg if num_categories==1
        for i in range(num_categories):
            C[i] = cmap( float(i)/(num_categories-1))
    return C

def add_categorical_legend(
        ax,
        categories,
        colors,
        legend_fontsize=None,
        ncol=None
):
    for idx, label in enumerate(categories):
        color = colors[idx]
        # use empty scatter to set labels
        ax.scatter([], [], c=[color], label=label)
    if ncol is None:
        ncol=(1 if len(categories) <= 20 # was 14
              else 2 if len(categories) <= 30 else 3)
    ax.legend(
        frameon=False, loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        ncol=ncol,
        fontsize=legend_fontsize,
        borderpad=0.0,
        handletextpad=0.0)
    # Shrink current axis by 10% to fit legend and match
    # size of plots that are not categorical
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.91, box.height])

def add_integers_legend( ax, color_dict, ncol=None ):
    categories = [ str(x) for x in sorted( color_dict.keys() ) ]
    colors = [ color_dict[int(x)] for x in categories ]
    add_categorical_legend( ax, categories, colors, ncol=ncol )

def make_rank_genes_logo_stack( ranks, upper_left, logo_width, max_logo_height,
                                top_pval_for_max_height = 1e-30,
                                min_pval_for_scaling = 1e-300,
                                signcolors=False,
                                num_genes_to_show = 10 ):
    ''' ranks is a list of (gene,l2r,pval)
    '''

    def pval_factor( pval, min_pval ):
        return math.sqrt( max(1e-6, -1 * math.log10( max(min_pval, pval) ) ))
        #return -1 * math.log10( max(min_pval,pval) )

    top_pval = ranks[0][2]
    logo_height = max_logo_height * ( pval_factor(top_pval, top_pval_for_max_height) /
                                      pval_factor(top_pval_for_max_height, top_pval_for_max_height) )

    if logo_height<1e-3:
        return []

    height_factors = [ pval_factor( x[2], min_pval_for_scaling ) for x in ranks[:num_genes_to_show] ]

    height_scale = logo_height / sum(height_factors)

    x0,y_offset = upper_left
    y_offset += ( max_logo_height - logo_height )

    cmds = []
    for gene,l2r,pval in ranks[:num_genes_to_show]:
        height = height_scale * pval_factor(pval, min_pval_for_scaling)
        if height<0.01:
            continue
        x1 = x0 + logo_width
        y0 = y_offset
        y1 = y0 + height
        if signcolors:
            if l2r <0:
                cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='blue' ) )
            else:
                cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='red' ) )
        else:
            if l2r < 1:
                cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='lightgray' ) )
            else:
                cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='black' ) )
        # elif l2r < 1:
        #     cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='blue' ) )
        # elif l2r < 1.5:
        #     cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='orange' ) )
        # else:
        #     cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='red' ) )

        y_offset+= height
    return cmds


def make_single_rank_genes_logo(
        ranks,
        svgfile,
        margin = 10,
        logo_width = 300,
        logo_max_height = 400,
        top_pval_for_max_height = 1e-50,
        min_pval_for_scaling = 1e-300,
        num_genes_to_show = 10,
        signcolors=False,
        create_png = True,
        create_html = False,
):
    cmds = []

    upper_left=  [margin, margin]
    cmds.extend( make_rank_genes_logo_stack(
        ranks, upper_left, logo_width, logo_max_height,
        top_pval_for_max_height=top_pval_for_max_height,
        min_pval_for_scaling=min_pval_for_scaling,
        signcolors=signcolors,
        num_genes_to_show=num_genes_to_show ) )

    svg_basic.create_file(
        cmds, logo_width+2*margin, logo_max_height + 2*margin, svgfile,
        create_png=create_png, create_html=create_html,
        background_color='white' )


def _make_lit_matches_logo(
        clone_indices,
        lit_matches,
        svgfile,
        margin = 10,
        logo_width = 300,
        max_logo_height = 400,
        create_png = True,
        create_html = False,
):
    clone_indices = set(clone_indices)
    mask = [x in clone_indices for x in lit_matches.clone_index]
    if np.sum(mask)==0:
        return False

    lit_matches = lit_matches[mask]

    matched_clone_indices = set()
    match2clone_indices = {}
    for l in lit_matches.itertuples():
        mhc = 'XXX' if pd.isna(l.db_mhc) else l.db_mhc
        epi = 'XXX' if pd.isna(l.db_epitope) else l.db_epitope
        if ':' in mhc:
            mhc = mhc[:mhc.index(':')]
        mhc = mhc.replace('*','')
        match = f'{mhc}-{epi[:3]}'
        match2clone_indices.setdefault(match,set()).add(l.clone_index)
        matched_clone_indices.add(l.clone_index)

    x_mhc, x_epi = 0.1, 0.1 # penalty for missing either info
    countslist = [
        (len(y) - x_mhc*x.startswith('XXX') - x_epi*x.endswith('XXX'), x)
        for x,y in match2clone_indices.items()
    ]
    countslist.sort(reverse=True)
    #print('_make_lit_matches_logo:', len(clone_indices), countslist)

    matched_fraction = len(matched_clone_indices)/len(clone_indices)
    total_height = max_logo_height * matched_fraction

    cmds = []
    y0 = margin + max_logo_height - total_height
    count_sum = sum(x[0] for x in countslist)
    for count, match in countslist:
        height = total_height * count / count_sum
        cmd = svg_basic.text_in_box(
            [margin, y0], [margin+logo_width,y0+height], match,
            color='black')
        cmds.append(cmd)
        y0+=height

    # # for starters: just take the top one
    # top_count, top_match = countslist[0]

    # logo_height = max_logo_height * top_count / len(clone_indices)
    # x0, y0, x1, y1 = (
    #     margin, margin+max_logo_height-logo_height,
    #     margin+logo_width, margin+max_logo_height
    #     )

    # cmd = svg_basic.text_in_box([x0,y0], [x1,y1], top_match, color='black')

    svg_basic.create_file(
        cmds, logo_width+2*margin, max_logo_height + 2*margin, svgfile,
        create_png=create_png, create_html=create_html,
        background_color='white')

    return True


def _parse_bicluster_rank_genes( adata, uns_tag = 'rank_genes_good_biclusters' ):
    ''' returns all_ranks dict mapping from cluspair to list [ (gene, l2r, pval), ... ]
    '''

    names = pd.DataFrame(adata.uns[uns_tag]['names'])
    organism = adata.uns['organism']

    all_ranks = {}
    for clptag in names.columns:
        assert clptag.startswith('clp_') and clptag.count('_') == 2
        clp = ( int(clptag.split('_')[1]), int(clptag.split('_')[2]) ) # gex_cluster, tcr_cluster
        ranks = []
        for igene, gene in enumerate( adata.uns[uns_tag]['names'][clptag] ):
            if pd.isnull(gene):
                continue
            if util.is_vdj_gene(gene, organism, include_constant_regions=True):
                continue  # NEWNEWNEW
            log2fold = adata.uns[uns_tag]['logfoldchanges'][clptag][igene]
            pval_adj = adata.uns[uns_tag]['pvals_adj'][clptag][igene]
            ranks.append( ( gene, log2fold, pval_adj ) )
        all_ranks[clp] = ranks
    return all_ranks



# def old_make_tcr_logo( tcrs, ab, organism, pngfile):
#     clonesfile = pngfile+'_clones.tsv'
#     util.make_clones_file( tcrs, clonesfile )

#     exe = '{}make_tcr_logo.py'.format( util.TCRDIST_REPO )
#     if not exists(exe):
#         print( 'ERROR: unable to locate python script in tcr-dist repo:', exe)
#         exit()

#     #> /dev/null 2> /dev/null'\
#     cmd = '{} {} --organism {} --ABs {} --outfile_prefix {} --clones_file {}'\
#         .format( util.PYTHON2_EXE, exe, organism, ab, pngfile[:-4], clonesfile )
#     util.run_command(cmd)
#     if not exists(pngfile):
#         print('make_tcr_logo:: cmd failed:', cmd)
#         exit()
#     os.remove(clonesfile)



def make_logo_plots(
        adata,
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size,
        logo_pngfile,
        logo_genes=None,
        gene_logo_width=6, # len(logo_genes) == 3*gene_logo_width-2
        clusters_gex= None,
        clusters_tcr= None,
        clusters_gex_names=None,
        clusters_tcr_names=None,
        ignore_tcr_cluster_colors = False,
        show_real_clusters_gex = False, # for the tcr clumping hack
        good_bicluster_tcr_scores=None,
        rank_genes_uns_tag = 'rank_genes_good_biclusters',
        include_alphadist_in_tcr_feature_logos=False,
        max_expn_for_gene_logo = 2.5, # or max over the clps, if larger
        show_pmhc_info_in_logos = False,
        nocleanup = False, # dont delete temporary image files (debugging)
        conga_scores = None,
        conga_scores_name = 'CoNGA', # for labeling the figure panels
        good_score_mask = None,
        make_batch_bars = None,
        batch_keys = None,
        make_cluster_gex_logos = True,
        draw_edges_between_conga_hits = True,
        add_conga_scores_colorbar = False,
        add_gex_logos_colorbar = False,
        pretty = False,

        ## controls for the gene expression thumbnails that come before
        ##  the actual logos:
        gex_header_genes=None,
        make_gex_header=True,
        make_gex_header_raw=True,
        make_gex_header_nbrZ=True,
        gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
        include_full_tcr_cluster_names_in_logo_lines=False,
        lit_matches=None, # show an additional 'logo' with lit-matches

        ## makes pdf version make_graph_vs_graph_logos and make_tcr_clumping_plots, expect large file sizes
        save_pdf = False

):
    ''' need:
    * gex/tcr clusters: (obsm)
    * clone sizes (obs: clone_sizes)
    * tcrs (obs)
    * 2d projections: gex and tcr (obsm: X_gex_2d, X_tcr_2d)
    * conga scores (obs: conga_scores) UNLESS passed in
    * rank genes info for each cluster-pair

    This code is a bit of a mess. OK, actually a huge mess.
    It started out simple but then when through
    multiple versions and complications and variable renamings. Sorry!

    '''

    if not make_gex_header:
        make_gex_header_raw, make_gex_header_nbrZ = False, False

    # show the distribution of the clones among the different batches
    if make_batch_bars is None:
        make_batch_bars = 'batch_keys' in adata.uns_keys()
    if make_batch_bars and batch_keys is None:
        batch_keys = adata.uns['batch_keys']

    ## unpack data from adata arrays ##################################
    clone_sizes = adata.obs['clone_sizes']
    if clusters_gex is None:
        clusters_gex = np.array(adata.obs['clusters_gex'])
        if clusters_gex_names is None:
            clusters_gex_names = adata.uns.get(
                'clusters_gex_names',
                [str(x) for x in range(np.max(clusters_gex)+1)])
    elif clusters_gex_names is None:
        clusters_gex_names = [str(x) for x in range(np.max(clusters_gex)+1)]
    if show_real_clusters_gex:
        real_clusters_gex = np.array(adata.obs['clusters_gex'])
        real_clusters_gex_names = adata.uns.get(
            'clusters_gex_names',
            [str(x) for x in range(np.max(real_clusters_gex)+1)])

    if clusters_tcr is None:
        clusters_tcr = np.array(adata.obs['clusters_tcr'])
        if clusters_tcr_names is None:
            clusters_tcr_names = adata.uns.get(
                'clusters_tcr_names',
                [str(x) for x in range(np.max(clusters_tcr)+1)])
    elif clusters_tcr_names is None:
        clusters_tcr_names = [str(x) for x in range(np.max(clusters_tcr)+1)]

    X_gex_2d = adata.obsm['X_gex_2d']
    X_tcr_2d = adata.obsm['X_tcr_2d']
    if conga_scores is None:
        conga_scores = np.array(adata.obs['conga_scores'])
    if good_score_mask is None:
        good_score_mask = np.array(conga_scores) <= 1.0
    organism = adata.uns['organism']

    if show_pmhc_info_in_logos:
        if 'X_pmhc' not in adata.obsm_keys():
            print('ERROR: include_pmhc_info_in_dendrogram=True but',
                  'no X_pmhc info in adata.obsm_keys()')
            sys.exit()
        pmhc_var_names = adata.uns['pmhc_var_names']
        X_pmhc = adata.obsm['X_pmhc']
        assert X_pmhc.shape == ( adata.shape[0], len(pmhc_var_names))

    if (lit_matches is None and
        'conga_results' in adata.uns_keys() and
        TCR_DB_MATCH in adata.uns['conga_results']):
        lit_matches = adata.uns['conga_results'][TCR_DB_MATCH]

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    ########################################## no more unpacking below here...
    #(actually we also unpack the logo/header genes GEX from adata.raw below)

    help_message = f"""
This figure summarizes the results of a CoNGA analysis that produces
scores ({conga_scores_name}) and clusters. At the top are six
2D UMAP projections of clonotypes in the dataset based on GEX similarity
(top left three panels) and TCR similarity (top right three panels),
colored from left to right by GEX cluster assignment;
{conga_scores_name} score; joint GEX:TCR cluster assignment for
clonotypes with significant {conga_scores_name} scores,
using a bicolored disk whose left half indicates GEX cluster and whose right
half indicates TCR cluster; TCR cluster; {conga_scores_name}; GEX:TCR cluster
assignments for {conga_scores_name} hits, as in the third panel.

Below are two rows of GEX landscape plots colored by (first row, left)
expression of selected marker genes, (second row, left) Z-score normalized and
GEX-neighborhood averaged expression of the same marker genes, and
(both rows, right) TCR sequence features (see CoNGA manuscript Table S3 for
TCR feature descriptions).

GEX and TCR sequence features of {conga_scores_name} hits in clusters with
{min_cluster_size} or more hits are summarized by a series
of logo-style visualizations, from left to right:
differentially expressed genes (DEGs); TCR sequence logos showing the V and
J gene usage and CDR3 sequences for the TCR alpha and beta chains; biased
TCR sequence scores, with red indicating elevated scores and blue indicating
decreased scores relative to the rest of the dataset (see CoNGA manuscript
Table S3 for score definitions); GEX 'logos' for each cluster
consisting of a panel of marker genes shown with red disks colored by
mean expression and sized according to the fraction of cells expressing
the gene (gene names are given above).

DEG and TCRseq sequence logos are scaled
by the adjusted P value of the associations, with full logo height requiring
a top adjusted P value below 10-6. DEGs with fold-change less than 2 are shown
in gray. Each cluster is indicated by a bicolored disk colored according to
GEX cluster (left half) and TCR cluster (right half). The two numbers above
each disk show the number of hits within the cluster (on the left) and
the total number of cells in those clonotypes (on the right). The dendrogram
at the left shows similarity relationships among the clusters based on
connections in the GEX and TCR neighbor graphs.

The choice of which marker genes to use for the GEX umap panels and for the
cluster GEX logos can be configured using run_conga.py command line flags
or arguments to the conga.plotting.make_logo_plots function.
    """

    num_clones = adata.shape[0]
    num_good_clones = np.sum(good_score_mask)

    if lit_matches is not None:
        if lit_matches.shape[0]==0:
            lit_matches = None
        else:
            required_fields = 'clone_index db_mhc db_epitope'.split()
            missing = False
            for f in required_fields:
                if f not in lit_matches.columns:
                    print('make_logo_plots: lit_matches missing required field',
                          f)
                    missing = True
            if missing:
                lit_matches = None

    if logo_genes is None:
        logo_genes = default_logo_genes[organism]

    if gex_header_genes is None:
        header2_genes = default_gex_header_genes[organism]
    else:
        header2_genes = gex_header_genes[:]
    if gex_header_tcr_score_names:
        header2_tcr_scores = tcr_scoring.make_tcr_score_table(
            adata, gex_header_tcr_score_names)
        assert header2_tcr_scores.shape == (adata.shape[0],
                                            len(gex_header_tcr_score_names))
    else:
        header2_tcr_scores = None

    # extract GEX info for the logo and header genes from the raw array
    raw_var_names = list(adata.raw.var_names)
    X_igex_genes = sorted(set(x for x in logo_genes+header2_genes
                              if x in raw_var_names))
    X_igex_indices = [raw_var_names.index(x) for x in X_igex_genes]
    X_igex = adata.raw[:,X_igex_indices].X.toarray()

    if 'clone_sizes' in header2_genes:
        X_igex = np.hstack(
            [X_igex, np.log(np.array(clone_sizes)[:,np.newaxis])])
        X_igex_genes.append('clone_sizes')
        assert X_igex.shape == (num_clones, len(X_igex_genes))
    if 'nndists_gex' in header2_genes:
        if 'nndists_gex' in adata.obs_keys():
            X_igex = np.hstack(
                [X_igex, np.array(adata.obs['nndists_gex'])[:,np.newaxis]])
            X_igex_genes.append('nndists_gex')
            assert X_igex.shape == (num_clones, len(X_igex_genes))
        else:
            print('WARNING nndists_gex not found in adata.obs')
            header2_genes.remove('nndists_gex')

    gene_width = gene_logo_width


    assert len(logo_genes) == 3*gene_width - 2

    # for making the tcr logos
    tcrdist_calculator = TcrDistCalculator(organism)


    # create nbrhood averaged X_igex array
    gex_nbrhood_X_igex = []
    gex_nbrhood_tcr_scores = []
    for ii in range(num_clones):
        # gex nbrhood
        nbrhood_mask = np.full( (num_clones,), False)
        nbrhood_mask[ nbrs_gex[ii] ] = True
        nbrhood_mask[ ii ] = True
        num = np.sum(nbrhood_mask)
        gex_nbrhood_X_igex.append(np.sum(X_igex[nbrhood_mask,:], axis=0)/num)
        if header2_tcr_scores is not None:
            gex_nbrhood_tcr_scores.append(
                np.sum(header2_tcr_scores[nbrhood_mask,:], axis=0)/num)

    gex_nbrhood_X_igex = np.array(gex_nbrhood_X_igex)
    assert gex_nbrhood_X_igex.shape == X_igex.shape
    if header2_tcr_scores is not None:
        gex_nbrhood_tcr_scores = np.array(gex_nbrhood_tcr_scores)
        assert gex_nbrhood_tcr_scores.shape == header2_tcr_scores.shape

    # set max expression levels for gex logo


    ## read info on louvains and tcr
    ## node2cluspair only includes good nodes
    node2cluspair = {
        i:(x,y) for i,(x,y,m) in enumerate(zip(clusters_gex,
                                               clusters_tcr,
                                               good_score_mask)) if m}
    num_clusters_gex = np.max(clusters_gex)+1
    num_clusters_tcr = np.max(clusters_tcr)+1

    cluspair2nodes = {}
    for node, clp in node2cluspair.items():
        cluspair2nodes.setdefault(clp,[]).append(node)


    # read the 2D projections
    all_xy_gex = { i:xy for i,xy in enumerate(X_gex_2d) }
    all_xy_tcr = { i:xy for i,xy in enumerate(X_tcr_2d) }

    # setup dbl nbrs
    all_dbl_nbrs = {}
    for ii in node2cluspair:
        gex_nbrs = frozenset( nbrs_gex[ii])
        all_dbl_nbrs[ii] = [ x for x in nbrs_tcr[ii]
                             if good_score_mask[x] and x in gex_nbrs ]


    ## parse the X_igex matrix
    all_scores = {} # map from clp to [ num_clusters_tcr, fracs, means]


    logo_gene_indices = [X_igex_genes.index(x) if x in X_igex_genes else None
                         for x in logo_genes]

    all_scores = {}
    for clp, nodes in cluspair2nodes.items():
        if len(nodes) < min_cluster_size:
            continue
        X = X_igex[nodes,:]
        means = np.mean(X,axis=0)
        fracs = np.sum(X>1e-3, axis=0) / float(len(nodes))

        fracs = [ 0 if x is None else fracs[x] for x in logo_gene_indices ]
        means = np.array([ 0.0 if x is None else means[x] for x in logo_gene_indices ])
        all_scores[clp] = [ len(nodes), fracs, means ]

    # if any of the means are bigger than max_expn_for_gene_logo, downscale that column
    max_means = np.maximum( np.max( np.array( [x[2] for x in all_scores.values()] ), axis=0 ), max_expn_for_gene_logo )
    means_rescale = max_expn_for_gene_logo / max_means

    for clp in all_scores:
        all_scores[clp][2] *= means_rescale

    # read rank_genes info
    if rank_genes_uns_tag is None:
        all_ranks = None
    else:
        all_ranks = _parse_bicluster_rank_genes( adata, uns_tag=rank_genes_uns_tag)


    clps = [x[1] for x in reversed(sorted([(y[0],x)
                                           for x,y in all_scores.items()]))]

    # aim for 10in width
    dendro_width = 1.0
    single_batch_bar_width = 0.25
    batch_bars_width = 0 if not make_batch_bars else \
                       single_batch_bar_width * len(batch_keys)
    title_logo_width = 0.75
    rg_logo_width = 1.5
    score_logo_width = 0.0 if good_bicluster_tcr_scores is None else 0.5
    lit_logo_width = 0.0 if lit_matches is None else 0.5
    tcr_logo_width = 8
    logo_height = 0.5
    frac_scale = 130
    margin = 0.2 # left and right margins actually
    top_margin = 0.2
    bottom_margin = 0.75 if make_batch_bars else 2.0*margin # was 1.3*margin
    yspacer_above_logos = 0.1
    yspacer_below_header = 0.15 if pretty else 0.0
    xspacer_within_header = 0.15 if pretty else 0.0

    gex_logo_width = logo_height * ( gene_width / (math.sqrt(3)+1) ) \
                     if make_cluster_gex_logos else 0.0

    fig_width = (2*margin
                 + dendro_width
                 + title_logo_width
                 + batch_bars_width
                 + rg_logo_width
                 + tcr_logo_width
                 + score_logo_width
                 + lit_logo_width
                 + gex_logo_width)
    conga_scores_colorbar_width = 0.035*fig_width if add_conga_scores_colorbar \
                                  else 0.0
    header_height = (fig_width
                     - conga_scores_colorbar_width
                     - xspacer_within_header
                     - 2*margin)/6.

    num_header2_rows = int(make_gex_header_nbrZ)+int(make_gex_header_raw)
    gex_logo_key_width = 2 if make_cluster_gex_logos else 0 # in header2 cols

    if make_gex_header:
        num_header2_tcr_score_cols = max(
            gex_logo_key_width,
            (gex_logo_key_width+len(gex_header_tcr_score_names))//num_header2_rows)
        gap_between_gex_and_tcr_panels = logo_height/10. \
                                         if num_header2_tcr_score_cols else 0
        header2_height = (fig_width-2*margin-gap_between_gex_and_tcr_panels)/\
                         (len(header2_genes) + num_header2_tcr_score_cols)

    else:
        header2_height = 0.

    fig_height = (
        top_margin +
        bottom_margin +
        header_height +
        num_header2_rows * header2_height +
        yspacer_below_header +
        (len(clps)>0)*yspacer_above_logos +
        len(clps)*logo_height
    )

    # frac_scale is in units of pts^2, and sqrt(frac_scale) is the MAX width in points of the marker
    # there are 72 points per inch
    # we want the markers to have a width of 1 data unit
    # so how many points in a data unit? the gex_logo is gene_width data units wide
    #
    if 2:
        fudge_factor = 0.85 # dunno why.
        marker_width_inches = gex_logo_width / gene_width
        marker_width_points = marker_width_inches*72.
        frac_scale = fudge_factor * (marker_width_points)**2
        #exit()

    plt.figure( figsize = (fig_width, fig_height ) )


    gex_colors = get_integers_color_dict( num_clusters_gex )
    tcr_colors = get_integers_color_dict( num_clusters_tcr )

    # set markersize for some small panel scatter plots
    if num_clones < 1000:
        small_markersize=5
    elif num_clones < 5000:
        small_markersize=4
    elif num_clones < 10000:
        small_markersize=3
    elif num_clones < 20000:
        small_markersize=2
    else:
        small_markersize=1.5

    # store distances between biclusters
    D = np.full( (len(clps),len(clps)), 1.0 ) # dist is 1.0 if no edges between bicluster
    for i in range(len(clps)):
        D[i,i] = 0.0

    ##############################################
    # make a header: 6 square scatter plots across the top
    # left to right:
    # 1 GEX UMAP, colored by GEX clusters
    # 2 GEX UMAP, colored by conga score
    # 3 GEX UMAP, only good cells, colored by cluster-pairs
    # 4 TCR UMAP, colored by TCR clusters
    # 5 TCR UMAP, colored by conga score
    # 6 TCR UMAP, only good cells, colored by cluster-pairs

    for icol, all_xy in enumerate( [all_xy_gex, all_xy_tcr] ):
        proj_tag = 'GEX TCR'.split()[icol]

        #############################################
        ## first a plot colored by the clusters
        if icol==0 or icol==1 and ignore_tcr_cluster_colors:
            if show_real_clusters_gex:
                clusters = real_clusters_gex
                cluster_names = real_clusters_gex_names
            else:
                clusters = clusters_gex
                cluster_names = clusters_gex_names
        else:
            clusters = clusters_tcr
            cluster_names = clusters_tcr_names
        bottom = (fig_height-top_margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+3*icol*header_height+icol*xspacer_within_header)/fig_width # make these plots squares
        width = header_height/fig_width
        plt.axes( [left,bottom,width,height] )
        num_clusters = np.max(clusters)+1
        C = get_integers_color_dict( num_clusters)
        ks = sorted( all_xy.keys() )
        xy = np.array( [ all_xy[x] for x in ks ] )
        colors = [ C[x] for x in clusters]
        plt.scatter( xy[:,0], xy[:,1], s=small_markersize, c=colors )
        plt.xticks([],[])
        plt.yticks([],[])
        clusters_tag = proj_tag if not ignore_tcr_cluster_colors else 'combo'
        if pretty:
            plt.title('{} clusters'.format(clusters_tag), fontsize=9, pad=1) #pad=8)
            plt.xlabel(f'{proj_tag} UMAP1', fontsize=7, labelpad=1)
            plt.ylabel(f'{proj_tag} UMAP2', fontsize=7, labelpad=0.0)
        else:
            plt.title('{} clusters ({}2D)'.format(clusters_tag, proj_tag), fontsize=9, pad=1) #pad=8)

        if icol==0:
            plt.text(0.01,0.01,'{} clonotypes'.format(len(all_xy)), ha='left', va='bottom', fontsize=8,
                     transform=plt.gca().transAxes)
            plt.text(0.99,0.01,'{} cells'.format(sum(clone_sizes)), ha='right', va='bottom', fontsize=8,
                     transform=plt.gca().transAxes)

        ## add cluster labels to the plot, try drawing at centroid location
        for c in set(clusters):
            cmask = (clusters==c)
            if (icol == 1 and show_real_clusters_gex and
                (c==0 or np.sum(cmask)<min_cluster_size)):
                # tcr_clumping hack
                continue
            centroid =np.mean(xy[cmask], axis=0)
            fontsize = 10 if len(cluster_names[c])<=2 else 8
            plt.text( centroid[0], centroid[1], cluster_names[c], color='k', fontsize=fontsize,
                      bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=1))

        #########################################3
        ## now a plot colored by conga score
        bottom = (fig_height-top_margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+(3*icol+1)*header_height+icol*xspacer_within_header)/fig_width # make these plots square
        width = header_height/fig_width
        plt.axes( [left,bottom,width,height] )

        ## make a plot colored by pval
        colors = np.array([max(0, -1*np.log10(conga_scores[x])) for x in ks])
        sort_order = np.argsort( colors )

        vmin = math.sqrt( -1*math.log10( 1.0 ) ) # was 0.05 then 0.2
        vmax = math.sqrt( -1*math.log10( 1e-5 ) )
        plt.scatter( xy[sort_order,:][:,0], xy[sort_order,:][:,1], s=2*small_markersize,
                     c=np.sqrt( colors[ sort_order ] ), vmin=vmin, vmax=vmax )
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.xlabel('{} UMAP1'.format(TAG))
        #plt.ylabel('{} UMAP2'.format(TAG), fontsize=12,labelpad=1)
        #plt.title('kNN-overlap E-values (Eval=Pval*num_clones)')
        # if add_conga_scores_colorbar:
        #     cbar = plt.colorbar(ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
        #     cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'])
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        if pretty:
            plt.title(f'{conga_scores_name} scores',fontsize=9,pad=1)
            plt.xlabel(f'{proj_tag} UMAP1', fontsize=7, labelpad=1)
        else:
            plt.title(f'{conga_scores_name} scores ({proj_tag}2D)',
                      fontsize=9, pad=1)
        #plt.title('kNN overlap E-values {}-UMAP'.format('GEX' if icol==0 else 'TCR'),fontsize=12,pad=1)


        ##############################################
        ## now just showing the top scoring clones
        bottom = (fig_height-top_margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+(3*icol+2)*header_height+icol*xspacer_within_header)/fig_width # make these plots square
        width = header_height/fig_width
        plt.axes( [left,bottom,width,height] )


        # first draw edges
        intra_lines = []
        inter_lines = []
        seen_line = set()
        for ii in node2cluspair:
            for jj in all_dbl_nbrs[ii]:
                if (jj,ii) not in seen_line:
                    seen_line.add((ii,jj))
                    if node2cluspair[ii] == node2cluspair[jj]:
                        intra_lines.append( [ [ all_xy[ii][0], all_xy[ii][1]], [ all_xy[jj][0], all_xy[jj][1] ] ] )
                    else:
                        inter_lines.append( [ [ all_xy[ii][0], all_xy[ii][1]], [ all_xy[jj][0], all_xy[jj][1] ] ] )

        if draw_edges_between_conga_hits:
            # the inter-cluspair lines
            light_gray = '#D3D3D3'
            line_segments = LineCollection( inter_lines, linewidths=0.25, zorder=1, colors=light_gray, alpha=0.25)
            plt.gca().add_collection(line_segments)

            # the intra-cluspair lines
            dark_gray = '#888888'
            line_segments = LineCollection( intra_lines, linewidths=0.25, zorder=1, colors='k')#dark_gray, alpha=0.5)
            plt.gca().add_collection(line_segments)


        for clp, nodes in cluspair2nodes.items():
            xvals = np.array([ all_xy[x][0] for x in nodes ])
            yvals = np.array([ all_xy[x][1] for x in nodes ])
            cscores = np.array([ conga_scores[x] for x in nodes ])

            if ignore_tcr_cluster_colors:
                plt.plot( xvals, yvals,
                          zorder=2,
                          marker='o',
                          linestyle='None',
                          color=gex_colors[clp[0]],
                          markeredgecolor='none',
                          markeredgewidth=0 )
            else:
                for score_threshold, markersize in [[1., 6],
                                                    [1e-1, 7],
                                                    [1e-2, 8]]:
                    mask = cscores <= score_threshold
                    if np.sum(mask)==0:
                        continue # nothing to plot
                    # most significant last:
                    reorder = np.argsort(cscores[mask])[::-1]
                    if show_real_clusters_gex:
                        left_color = 'none'
                    else:
                        left_color = gex_colors[clp[0]]
                    plt.plot(xvals[mask][reorder], yvals[mask][reorder],
                             zorder=2,
                             marker='o',
                             linestyle='None',
                             markersize=markersize,
                             color=left_color,
                             markerfacecoloralt=tcr_colors[clp[1]],
                             fillstyle='left',
                             markeredgecolor='none',
                             markeredgewidth=0 )
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.ylabel('{} UMAP2'.format(TAG), fontsize=12,labelpad=1)
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        if pretty:
            plt.title(f'{conga_scores_name} hits', fontsize=9, pad=1) #pad=8)
            plt.xlabel(f'{proj_tag} UMAP1', fontsize=7, labelpad=1)
        else:
            plt.title(f'{conga_scores_name} hits ({proj_tag}2D)',
                      fontsize=9, pad=1)

        if add_conga_scores_colorbar and icol==1:
            ax_width = conga_scores_colorbar_width/8.
            ax_height = 0.72*header_height
            bottom = (fig_height-top_margin-header_height+0.08*header_height)/fig_height
            height = ax_height/fig_height
            left = (margin+6*header_height+xspacer_within_header+0.5*ax_width)/fig_width
            width = ax_width/fig_width
            plt.axes( [left,bottom,width,height] )
            sm = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax),
                cmap='viridis')
            sm.set_array(np.sqrt(np.array([max(0, -1*np.log10(conga_scores[x]))
                                           for x in ks ])))

            #plt.scatter([],[],c=[],cmap='viridis', vmin=vmin, vmax=vmax)
            #plt.gca().set_visible(False)
            cbar = plt.colorbar(mappable=sm, cax=plt.gca(),
                                ticks=[math.sqrt(-1*math.log10(x))
                                       for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5]])
            # cbar = plt.colorbar(
            #     cax=plt.gca(), ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
            cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'], fontsize=7)
            #cbar.ax.set_ylabel('CoNGA score')
            cbar.ax.set_title(f'{conga_scores_name}\nscore', loc='left',
                              ha='left', fontsize=9)



        # for making the dendrogram
        cluspair_edge_counts = Counter()
        for ii, dbl_nbrs in all_dbl_nbrs.items():
            for jj in dbl_nbrs:
                cluspair_edge_counts[ (node2cluspair[ii], node2cluspair[jj]) ] += 1
                cluspair_edge_counts[ (node2cluspair[jj], node2cluspair[ii]) ] += 1 # since dbl_nbrs not nec sym

        for (clp1,clp2),count in cluspair_edge_counts.items():
            if clp2<=clp1:
                assert count == cluspair_edge_counts[(clp2,clp1)]
                continue
            sz1 = len(cluspair2nodes[clp1])
            sz2 = len(cluspair2nodes[clp2])
            if sz1<min_cluster_size or sz2<min_cluster_size:
                continue
            # how thick to make the edge?
            frac = float(count)/(sz1*sz2*2) # since added both ways...
            dist = 1.0 - frac
            D[ clps.index(clp1), clps.index(clp2) ] = dist
            D[ clps.index(clp2), clps.index(clp1) ] = dist




    if make_gex_header:
        # make a line of GEX plots, right now for the logo_genes
        for ig, gene in enumerate(header2_genes):
            if gene not in X_igex_genes:
                print('missing header2_gene:', gene,'from X_igex_genes')
                continue
            index = X_igex_genes.index(gene)

            if make_gex_header_raw:
                bottom = (fig_height-top_margin-header_height-header2_height-yspacer_below_header)/fig_height
                height = header2_height/fig_height
                left = (margin+ig*header2_height)/fig_width # make these plots square, so height=width
                width = header2_height/fig_width
                plt.axes( [left,bottom,width,height] )
                colors = X_igex[:,index]
                reorder = np.argsort(colors)
                vmax = max(1.0, np.max(colors))
                plt.scatter( X_gex_2d[:,0][reorder], X_gex_2d[:,1][reorder], s=small_markersize, c=colors[reorder],
                             cmap='Reds', vmin=0, vmax=vmax )
                plt.xticks([],[])
                plt.yticks([],[])
                xmn,xmx = plt.xlim()
                _,ymx = plt.ylim()
                pad = 0.02*(xmx-xmn)
                if ig==0:
                    plt.ylabel('GEX UMAP2', fontsize=6, labelpad=1)
                    if not make_gex_header_nbrZ:
                        plt.xlabel('GEX UMAP1', fontsize=6, labelpad=1)
                if make_gex_header_nbrZ:
                    plt.text( xmn+pad, ymx, 'raw', va='top', ha='left', fontsize=8)
                    gene_text = gene # no space for both 'clonotype_size' and 'raw'
                else:
                    gene_text = gene.replace('clone_sizes', 'clonotype_size')
                plt.text( xmx, ymx, gene_text, va='top', ha='right', fontsize=8)

            if make_gex_header_nbrZ: ## make another plot with the gex-nbrhood averaged scores
                bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height-yspacer_below_header)/fig_height
                height = header2_height/fig_height
                left = (margin+ig*header2_height)/fig_width # make these plots square, so height=width
                width = header2_height/fig_width
                plt.axes( [left,bottom,width,height] )
                vals = X_igex[:,index]
                mean, std = np.mean(vals), np.std(vals)
                if std<0.01:
                    print('low_std:',gene,std)
                    std = 0.01
                colors = (gex_nbrhood_X_igex[:,index]-mean)/std

                mx = max( 0.2, np.max(np.abs(colors)))
                plt.scatter( X_gex_2d[:,0], X_gex_2d[:,1], s=small_markersize, c=colors,cmap='coolwarm',
                             vmin=-mx, vmax=mx )
                plt.xticks([],[])
                plt.yticks([],[])
                xmn,xmx = plt.xlim()
                _,ymx = plt.ylim()
                pad = 0.02*(xmx-xmn)
                if ig==0:
                    plt.xlabel('GEX UMAP1', fontsize=6, labelpad=1)
                    plt.ylabel('GEX UMAP2', fontsize=6, labelpad=1)
                plt.text( xmn+pad, ymx, 'nbr', va='top', ha='left', fontsize=8)
                plt.text( xmx, ymx, gene, va='top', ha='right', fontsize=8)

        if gex_header_tcr_score_names:
            for ii, scoretag in enumerate(gex_header_tcr_score_names):
                irow=ii//num_header2_tcr_score_cols
                icol=ii%num_header2_tcr_score_cols
                ## make another plot with the gex-nbrhood averaged TCR scores
                bottom = (fig_height-top_margin-header_height-(irow+1)*header2_height-yspacer_below_header)/fig_height
                height = header2_height/fig_height
                left = (fig_width-margin-(num_header2_tcr_score_cols-icol)*header2_height)/fig_width
                width = header2_height/fig_width
                plt.axes( [left,bottom,width,height] )
                vals = header2_tcr_scores[:,ii]
                mean, std = np.mean(vals), np.std(vals)
                if std<0.01:
                    print('low_std:', scoretag, std)
                    std = 0.01
                colors = (gex_nbrhood_tcr_scores[:,ii]-mean)/std

                mx = max( 0.2, np.max(np.abs(colors)))
                plt.scatter( X_gex_2d[:,0], X_gex_2d[:,1], s=small_markersize, c=colors,cmap='coolwarm',
                             vmin=-mx, vmax=mx )
                plt.xticks([],[])
                plt.yticks([],[])
                xmn,xmx = plt.xlim()
                _,ymx = plt.ylim()
                pad = 0.02*(xmx-xmn)
                plt.text( xmx-pad, ymx,
                          scoretag.replace('alphadist','alphaD').replace('nndists_tcr','NNd').replace('cdr3len','len'),
                          va='top', ha='right', fontsize=8)
                plt.text( xmn+pad, ymx, 'tcr', va='top', ha='left', fontsize=8)



    if header2_height > 0.001 and len(clps)>0 and make_cluster_gex_logos: ## make a key for the gex logo
        downscale = 0.84 if add_gex_logos_colorbar else 1.0
        bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height-yspacer_below_header)/fig_height
        height = downscale * header2_height/fig_height
        left = (fig_width-margin-2*header2_height)/fig_width # make these plots square, so height=width
        width = downscale * 2*header2_height/fig_width
        plt.axes( [left,bottom,width,height] )
        yshift = math.sqrt(3)/2 # 30 60 90 triangle
        xy = []
        for r in range(3):
            xoff = 0 if r==1 else 0.5
            yoff = yshift * ( 1-r )
            ng = gene_width if r==1 else gene_width-1
            for ii in range(ng):
                xy.append( [xoff+ii, yoff] )
        xy = np.array(xy)
        color = plt.get_cmap('Reds')(0.2)
        plt.scatter(xy[:,0], xy[:,1], c=[color], s=310*downscale, zorder=1) # was 350
        for gene,(x,y) in zip( logo_genes, xy ):
            plt.text(x,y,gene,fontsize=7,ha='center',va='center',rotation=30, zorder=2) # was 8
        plt.xlim( [-0.5, gene_width-0.5] )
        plt.ylim( [-yshift-0.5, yshift+0.5 ] )
        plt.xticks([],[])
        plt.yticks([],[])
        if add_gex_logos_colorbar:
            plt.title('GEX logo genes', fontsize=8, va='bottom', pad=1)

            ax_width = 0.06*width
            ax_height = 0.75*height
            ax_bottom = bottom + 0.05*height
            #bottom = (fig_height-top_margin-header_height+0.08*header_height)/fig_height
            #height = ax_height/fig_height
            ax_left = left+width+0.05*width
            plt.axes( [ax_left,bottom,ax_width,ax_height] )
            sm = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=0, vmax=1), cmap='Reds')
            sm.set_array(np.linspace(vmin, vmax, 10))

            #plt.scatter([],[],c=[],cmap='viridis', vmin=vmin, vmax=vmax)
            #plt.gca().set_visible(False)
            cbar = plt.colorbar(mappable=sm, cax=plt.gca(), ticks = [0,1])
            #ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
            # cbar = plt.colorbar(
            #     cax=plt.gca(), ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
            cbar.ax.set_yticklabels(['0', 'max'], fontsize=7)
            cbar.ax.set_title('Mean\nexpn', loc='left', ha='left', fontsize=8)


    if len(clps)>1:

        # make a dendrogram along the LHS
        bottom = bottom_margin/fig_height # there's an extra margin-width at the bottom
        height = (fig_height-top_margin-bottom_margin-header_height-num_header2_rows*header2_height)/fig_height
        left = (margin)/fig_width
        width = dendro_width/fig_width
        ax = plt.axes( [left,bottom,width,height] )
        tree_method = 'average'
        tree_optimal_ordering = False

        Z = hierarchy.linkage( distance.squareform(D,force='tovector'),
                               method=tree_method, optimal_ordering=tree_optimal_ordering )

        R = hierarchy.dendrogram( Z, orientation='left', ax=ax, link_color_func= lambda x:'k' )
        def get_leafnum_from_y(y):
            " 5.0 15.0 25.0 --> 0 1 2 "
            return int((y-4.99)/10)

        def get_thickness_for_clp(clp):
            lw_perc_scale = 8.0
            max_lw = 25
            top_clp_size = max(len(x) for x in cluspair2nodes.values())
            top_lw = lw_perc_scale * 100 * top_clp_size / num_clones
            if top_lw > max_lw:
                lw_perc_scale *= max_lw/top_lw
                #print('NOTE: big clusters, decrease lw_perc_scale from 5.0 to',lw_perc_scale)
            clp_size = len(cluspair2nodes[clp])
            perc = 100*clp_size/num_clones # percent of all clones in this clust-pair
            # a big perc is like 5%
            lw = max( 1, lw_perc_scale * perc)
            #print(lw, perc, clp)
            return lw


        for xs,ys in zip(R['dcoord'], R['icoord']):
            for i0,i1 in [[0,1],[3,2]]:
                if xs[i0] == 0:# a terminal edge
                    y = ys[i0]
                    num = get_leafnum_from_y(y) # 0 at the bottom to numleaves-1 at the top
                    clp = clps[ R['leaves'][num] ] # R['leaves'] is indices into clps, numbered from bottom
                    lw = get_thickness_for_clp(clp)
                    #lw = 1
                    plt.plot( [0,xs[i1]], [y, y], '-k', linewidth=lw, solid_capstyle='butt')

        plt.xlim((1.03,0.0))
        plt.axis('off')
        plt.text(0.5, 0.25/len(clps),
                 'Clusters\n(size>{:d})'.format(int(min_cluster_size)-1),
                 ha='center', va='top', transform=plt.gca().transAxes)
        # plt.text(0.0,0.0, 'Clusters\n(size>{:d})'.format(int(min_cluster_size)-1),
        #          ha='left', va='top', transform=plt.gca().transAxes)
        leaves = R['leaves'][:] #list( hierarchy.leaves_list( Z ) )
        leaves.reverse() # since we are drawing them downward, but the leaf-order increases upward
    elif clps:
        leaves = [0]
    else:
        leaves = []

    tmpfiles = []
    for irow, iclp in enumerate(leaves):
        first_row = ( irow == 0)
        last_row = ( irow == len(leaves)-1 )
        print('making cluster logos:',irow,len(leaves),logo_pngfile)
        clp = clps[iclp]
        bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height
                  -yspacer_below_header-yspacer_above_logos - (irow+1)*logo_height)/fig_height
        height = logo_height/fig_height

        num_nodes, fracs, means = all_scores[clp]
        nodes = cluspair2nodes[clp]
        assert len(nodes) == num_nodes
        num_cells = sum( clone_sizes[x] for x in nodes )

        # make the title logo
        left = (margin+dendro_width)/fig_width
        width = title_logo_width/fig_width
        plt.axes( [left,bottom,width,height] )
        if ignore_tcr_cluster_colors:
            plt.plot([0.0], [0.0], marker='o', linestyle='None',
                     color=gex_colors[clp[0]], markeredgecolor='none',
                     markeredgewidth=0, markersize=30)
        else:
            if show_real_clusters_gex:
                # for the tcr clumping hack
                plt.plot([0.0], [0.0], marker='o', linestyle='None',
                         color='white',
                         markerfacecoloralt=tcr_colors[clp[1]],
                         fillstyle='left', markeredgecolor='none',
                         markeredgewidth=0, markersize=30, zorder=1)

                # make a vertical bar plot showing the gex clusters
                clp_cluscounts = Counter(real_clusters_gex[x] for x in nodes)
                numclus = np.max(real_clusters_gex)+1
                C = get_integers_color_dict(numclus)
                colors = [C[x] for x in range(numclus)]
                fractions = np.array([clp_cluscounts[x]/len(nodes)
                                      for x in range(numclus)])
                stack_width = 0.54
                plt.bar([-0.03]*numclus, width=-stack_width,
                        height=1.8*fractions,
                        bottom=-0.9+1.8*(np.cumsum(fractions)-fractions),
                        color=colors, zorder=2,
                        align='edge') # w/ neg width --> right alignment
                # label explaining that these bars show GEX clusters
                if last_row:
                    plt.text(-0.25, -0.95, 'GEX clust', rotation=45,
                             va='top', ha='right', fontsize=7)

            else:
                plt.plot([0.0], [0.0], marker='o', linestyle='None',
                         color=gex_colors[clp[0]],
                         markerfacecoloralt=tcr_colors[clp[1]],
                         fillstyle='left', markeredgecolor='none',
                         markeredgewidth=0, markersize=30)
        # if gex_cluster_names is None:
        #     plt.text(-0.1,0,str(clp[0]),va='center',ha='right')
        # else:
        short_name_fontsize, long_name_fontsize = 11, 6
        if not show_real_clusters_gex:
            name = clusters_gex_names[clp[0]]
            if len(name) <= 2:
                plt.text(-0.1,0,name,va='center',ha='right',
                         fontsize=short_name_fontsize)
            else:
                plt.text(-0.1,0,name,va='center',ha='right',
                         fontsize=long_name_fontsize,
                         bbox=dict(facecolor='white', alpha=0.5,
                                   edgecolor='white', pad=1))

        if not ignore_tcr_cluster_colors:
            if include_full_tcr_cluster_names_in_logo_lines:
                name = clusters_tcr_names[clp[1]]
                if len(name)<=2:
                    plt.text(0.1, 0, name, va='center', ha='left',
                             fontsize=short_name_fontsize)
                else:
                    plt.text(0.1, 0, name, va='center', ha='left',
                             fontsize=long_name_fontsize,
                             bbox=dict(facecolor='white', alpha=0.5,
                                       edgecolor='white', pad=1))
            else:
                plt.text( 0.1,0,str(clp[1]),va='center',ha='left',
                          fontsize=short_name_fontsize)
        plt.text(-1,1,'{}'.format(num_nodes),va='top',ha='left',fontsize=8)
        plt.text( 1,1,'{}'.format(num_cells),va='top',ha='right',fontsize=8)
        if first_row:
            plt.text( 0.0, 1, 'clonotypes', va='bottom', ha='right', fontsize=6)
            plt.text( 1.0, 1, 'cells', va='bottom', ha='right', fontsize=6)

        if show_pmhc_info_in_logos:
            X_pmhc_clp = np.sum(X_pmhc[nodes,:],axis=0)/len(nodes)
            sorted_pmhc_indices = np.argsort(X_pmhc_clp)[::-1]
            pmhc0, avglog1pcount0 = pmhc_var_names[sorted_pmhc_indices[0]], X_pmhc_clp[sorted_pmhc_indices[0]]
            pmhc1, avglog1pcount1 = pmhc_var_names[sorted_pmhc_indices[1]], X_pmhc_clp[sorted_pmhc_indices[1]]
            plt.text(-1,0,'{}\n{}\n{:.2f}'.format(pmhc0[:3], pmhc0[4:7], avglog1pcount0),
                     va='center',ha='left',fontsize=6)
            plt.text( 1,0,'{}\n{}\n{:.2f}'.format(pmhc1[:3], pmhc1[4:7], avglog1pcount1),
                      va='center',ha='right',fontsize=6)

        plt.xlim((-1,1))
        plt.ylim((-1,1))
        plt.axis('off')

        # show the batch distribution
        if make_batch_bars:
            for ib, k in enumerate(batch_keys):
                batch_counts_for_clp = np.array(adata.obsm[k])[nodes,:]
                left = (margin+dendro_width+title_logo_width+ib*single_batch_bar_width)/fig_width
                width = single_batch_bar_width/fig_width
                plt.axes( [left,bottom,width,height] )
                # bar plot
                counts = np.sum(batch_counts_for_clp, axis=0)
                num_batch_key_choices = counts.shape[0]
                assert num_batch_key_choices >1 # we add a fake one in preprocess.py if necessary
                fractions = counts.astype(float)/np.sum(counts)
                cmap = plt.get_cmap('tab10') if num_batch_key_choices<=10 else \
                       plt.get_cmap('tab20')
                colors = [cmap.colors[x%20]
                          for x in range(num_batch_key_choices)]
                plt.bar([0]*num_batch_key_choices, height=fractions, width=0.8,
                        bottom=np.cumsum(fractions)-fractions, color=colors, align='center')
                # overall counts/fractions
                batch_counts_for_clp = np.array(adata.obsm[k])[:,:]
                counts = np.sum(batch_counts_for_clp, axis=0)
                fractions = counts.astype(float)/np.sum(counts)
                plt.bar([0.6]*num_batch_key_choices, height=fractions, width=0.4,
                        bottom=np.cumsum(fractions)-fractions, color=colors, align='center')
                plt.axis('off')
                plt.ylim((0,1.05))

                # label with the batch_key
                if last_row:
                    plt.text(0.6, -0.05, k, rotation=45, va='top', ha='right',
                             fontsize=7)



        # make the rank genes logo
        if all_ranks is not None:
            clp_rank_genes = all_ranks[clp]
            for r in range(3):
                pngfile = '{}.tmp_{}_{}_{}.png'.format(logo_pngfile, clp[0], clp[1], r)
                start_rank, stop_rank = 3*r, 3*r+3
                left = (margin+dendro_width+title_logo_width+batch_bars_width+r*0.333*rg_logo_width)/fig_width
                width = 0.333*rg_logo_width/fig_width
                plt.axes( [left,bottom,width,height] )
                if len(clp_rank_genes) > start_rank and clp_rank_genes[start_rank][2] < 0.99:
                    # if none of the pvals are better than 1 then we get an all white image--> black box (?)
                    make_single_rank_genes_logo( clp_rank_genes[start_rank:stop_rank], pngfile[:-3]+'svg',
                                                 logo_width=800, logo_max_height=1000, top_pval_for_max_height=1e-6)

                    image = mpimg.imread(pngfile)
                    tmpfiles.append(pngfile)
                    plt.imshow(image, cmap='Greys_r')
                    xmn,xmx = plt.xlim()
                    ymn,ymx = plt.ylim()
                    #plt.plot([xmn,xmx],[ymn,ymn],'-k')
                    plt.plot([xmn,xmx],[ymx,ymx],'-k')
                    plt.xlim(xmn,xmx)
                    plt.ylim(ymn,ymx)
                else:
                    plt.scatter([],[])
                plt.axis('off')
                if last_row and r==1:
                    plt.text(0.5,-0.05, 'Top 9 DE genes logo', ha='center', va='top', transform=plt.gca().transAxes)

            #plt.xticks([],[])
            #plt.yticks([],[])

        # make a tcr logo
        for iab,ab in enumerate('AB'):
            left = (margin+dendro_width+title_logo_width+batch_bars_width+rg_logo_width+iab*0.5*tcr_logo_width)/fig_width
            width = 0.5*tcr_logo_width/fig_width

            plt.axes( [left,bottom,width,height] )

            pngfile = '{}_tmp_{}_{}_{}.png'.format(logo_pngfile, clp[0], clp[1], ab)
            tmpfiles.append(pngfile)
            # if 1: # old way
            #     old_make_tcr_logo( [ tcrs[x] for x in nodes ], ab, organism, pngfile )
            # else: # new way
            make_tcr_logo_for_tcrs( [ tcrs[x] for x in nodes ], ab, organism, pngfile,
                                    tcrdist_calculator=tcrdist_calculator )
            image = mpimg.imread(pngfile)
            plt.imshow(image)
            plt.axis('off')
            if last_row:
                plt.text(0.5,-0.05, 'TCR{} logo'.format(ab), ha='center', va='top', transform=plt.gca().transAxes)


        # make the tcr_scores logo
        if good_bicluster_tcr_scores is not None:
            clp_rank_scores = good_bicluster_tcr_scores[clp]
            if DONT_SHOW_LOW_MAIT_SCORES:
                clp_rank_scores = [ x for x in clp_rank_scores if x[0] != 'mait' or x[1]>0 ]
            if not include_alphadist_in_tcr_feature_logos:
                clp_rank_scores = [ x for x in clp_rank_scores if x[0] != 'alphadist' ]
            pngfile = '{}.tmpsc_{}_{}.png'.format(logo_pngfile, clp[0], clp[1])
            left = (margin+dendro_width+title_logo_width+batch_bars_width+rg_logo_width+tcr_logo_width)/fig_width
            width = score_logo_width/fig_width
            plt.axes( [left,bottom,width,height] )
            if len(clp_rank_scores):
                # if none of the pvals are better than 1 then we get an all white image--> black box (?)
                make_single_rank_genes_logo( clp_rank_scores[:3], pngfile[:-3]+'svg',
                                             logo_width=900, logo_max_height=1000, top_pval_for_max_height=1e-6,
                                             signcolors=True)

                #image = mpimg.imread(pngfile)
                image = plt.imread(pngfile)
                tmpfiles.append(pngfile)
                plt.imshow(image)
                xmn,xmx = plt.xlim()
                ymn,ymx = plt.ylim()
                #plt.plot([xmn,xmx],[ymn,ymn],'-k')
                plt.plot([xmn,xmx],[ymx,ymx],'-k')
                plt.xlim(xmn,xmx)
                plt.ylim(ymn,ymx)
            else:
                plt.scatter([],[])
            plt.axis('off')
            plt.xticks([],[])
            plt.yticks([],[])
            if last_row:
                plt.text(0.9, -0.05, 'TCRseq\nfeatures', ha='right', va='top', fontsize=8,
                         transform=plt.gca().transAxes)

        if lit_matches is not None:
            # maybe make a lit-matches logo for this cluster
            left = (margin+dendro_width+title_logo_width+batch_bars_width
                    +rg_logo_width+tcr_logo_width+score_logo_width)/fig_width
            width = lit_logo_width/fig_width

            plt.axes( [left,bottom,width,height] )
            pngfile = '{}.tmp_{}_{}_lit.png'\
                      .format(logo_pngfile, clp[0], clp[1])

            success = _make_lit_matches_logo(
                nodes, lit_matches, pngfile[:-3]+'svg',
                logo_width=1900*lit_logo_width,
                max_logo_height=2000*logo_height)

            if success: # there were lit matches for this cluster
                image = mpimg.imread(pngfile)
                tmpfiles.append(pngfile)
                plt.imshow(image, cmap='Greys_r')
            plt.axis('off')
            if last_row:
                plt.text(0.5,-0.05, 'Lit.DB\nmatch', ha='center', va='top',
                         transform=plt.gca().transAxes, fontsize=8)



        if make_cluster_gex_logos: # make the gex logo
            left = (margin+dendro_width+title_logo_width+batch_bars_width
                    +rg_logo_width+tcr_logo_width+lit_logo_width
                    +score_logo_width)/fig_width
            width = gex_logo_width/fig_width

            plt.axes( [left,bottom,width,height] )

            yshift = math.sqrt(3)/2 # 30 60 90 triangle

            xy = []
            for r in range(3):
                xoff = 0 if r==1 else 0.5
                yoff = yshift * ( 1-r )
                ng = gene_width if r==1 else gene_width-1
                for ii in range(ng):
                    xy.append( [xoff+ii, yoff] )

            xy = np.array(xy)
            sizes = [ x*frac_scale for x in fracs ]
            colors = means

            plt.scatter( xy[:,0], xy[:,1], s=sizes, c=colors, cmap = plt.get_cmap('Reds'), vmin=0,
                         vmax=max_expn_for_gene_logo )
            plt.xlim( [-0.5, gene_width-0.5] )
            plt.ylim( [-yshift-0.5, yshift+0.5 ] )
            plt.xticks([],[])
            plt.yticks([],[])
            if last_row:
                plt.text(0.5,-0.05, 'GEX logo', ha='center', va='top', transform=plt.gca().transAxes, fontsize=9)
                plt.text(0.5,-0.32, '(radius ~ % pos. cells)', ha='center', va='top', fontsize=7,
                         transform=plt.gca().transAxes)
                # msg = 'GEX logo genes: row1= {} row2= {} row3= {}'\
                #     .format(','.join(logo_genes[:gene_width-1]),
                #             ','.join(logo_genes[gene_width-1:2*gene_width-1]),
                #             ','.join(logo_genes[2*gene_width-1:]) )
                # plt.text(1.0,-0.5, msg, ha='right', va='top', transform=plt.gca().transAxes, fontsize=7)



    print('making:', logo_pngfile)
    plt.savefig(logo_pngfile, dpi=300)

    if save_pdf:
        logo_pdffile = logo_pngfile.replace(".png", ".pdf")
        print('making:', logo_pdffile)
        plt.savefig(logo_pdffile, dpi=300)

    if not nocleanup:
        for tmpfile in tmpfiles:
            if exists(tmpfile):
                os.remove(tmpfile)
            svgfile = tmpfile[:-3]+'svg'
            if exists(svgfile):
                os.remove(svgfile)


    return help_message


def make_n_pseudopoints( n, xy, radius_in_sdevs=0.25 ):
    ''' Find the centroid of the 2D xy array, pick n points around there
    '''
    assert xy.shape[1] == 2
    center = np.mean(xy, axis=0)
    sdev = np.sqrt( np.sum( np.square(xy-center[np.newaxis,:]))/xy.shape[0])
    radius = radius_in_sdevs * sdev
    #print('make_n_pseudopoints:', n, sdev, radius, center)
    rots = np.linspace(1, 1+2*np.pi, n+1)[:-1] # stagger things a bit
    vecs = np.vstack( [ radius*np.cos(rots), radius*np.sin(rots) ] ).transpose()
    # vecs = np.zeros((1,2)) # at center
    # if n>1:
    #     sdev = np.sqrt( np.sum( np.square(xy-center[np.newaxis,:]))/xy.shape[0])
    #     radius = radius_in_sdevs * sdev
    #     print('make_n_pseudopoints:', n, sdev, radius, center)
    #     rots = np.linspace(0, 2*np.pi, n)[:-1]
    #     new_vecs = np.vstack( [ radius*np.cos(rots), radius*np.sin(rots) ] ).transpose()
    #     vecs = np.vstack( [ vecs, new_vecs] )
    assert vecs.shape == (n,2)
    points = vecs + center[np.newaxis,:]
    return points


def plot_ranked_strings_on_cells(
        adata,
        results_df,
        xy_tag, # for retrieving from obsm
        clone_index_column, # for mapping to xy array
        ranking_column, # lower is better
        ranking_threshold, # larger than this, we won't show the text
        string_column,
        pngfile=None,
        exclude_strings = [],
        direction_column=None, # for diverging plots
        color_text=False,
        unique_by_cluster=True, ## requires cluster information
        gex_cluster_column='gex_cluster',
        tcr_cluster_column='tcr_cluster',
        strings_per_cluster=3,
        figsize=(6,6),
        ax=None,
        title=None,
):
    ''' Visualize the location of neighborhoods with skewed feature scores
    on a UMAP plot, labeled by the feature name

    returns a help message
    '''
    assert ax or pngfile

    if results_df.shape[0] == 0: # no results to plot
        print('plot_ranked_strings_on_cells: empty results')
        return


    # new wrinkle: clone_index column may be -1, if we included the
    # cluster-level results... figure out if these are gex or tcr clusters
    if 'gex' in xy_tag:
        clusters = np.array(adata.obs['clusters_gex'])
        clp_index_for_minus1 = 0
    else:
        assert 'tcr' in xy_tag
        clusters = np.array(adata.obs['clusters_tcr'])
        clp_index_for_minus1 = 1


    # unpack data
    xy = adata.obsm[xy_tag]
    cols = [ranking_column, clone_index_column, string_column, string_column]
    if direction_column is not None:
        cols[3] = direction_column
    if unique_by_cluster:
        cols.extend([gex_cluster_column,tcr_cluster_column])

    df = results_df[ cols ].sort_values(ranking_column, ascending=True)

    # not sure why but sometimes these integers get converted to floats...
    df[clone_index_column] = df[clone_index_column].astype(int)
    if unique_by_cluster:
        for col in [gex_cluster_column, tcr_cluster_column]:
            df[col] = df[col].astype(int)


    #first restrict to a single text per cell
    mask = []
    seen = set()
    if unique_by_cluster:
        clp_count=Counter()
        # max one text per cell
        # max strings_per_cluster texts per cluster
        # no string repeats for cluster
        for row in df.itertuples():
            pval = row[1]
            ii = row[2]
            gene = row[3]
            clp = (row[5], row[6])
            if ii==-1: # for these just use a single integer as clp
                clp = clp[clp_index_for_minus1]
            if (DONT_SHOW_LOW_MAIT_SCORES and gene=='mait' and 'gex' in xy_tag
                and row[4]<0):
                mask.append(False)
            elif (pval > ranking_threshold or gene in exclude_strings or
                  (ii!=-1 and ii in seen) or (gene,clp) in seen or
                  clp_count[clp] >= strings_per_cluster):
                mask.append(False)
            else:
                seen.add(ii)
                seen.add((gene,clp))
                clp_count[clp] += 1
                mask.append(True)
    else:
        for row in df.itertuples():
            pval = row[1]
            ii = row[2]
            gene = row[3]
            if pval>ranking_threshold or gene in exclude_strings or ii in seen:
                mask.append(False)
            else:
                seen.add(ii)
                mask.append(True)
    df = df[mask]

    # now sort so that the best stuff is at the end
    df.sort_values(ranking_column,ascending=False,inplace=True)


    if pngfile:
        plt.figure(figsize=figsize)
    else:
        plt.sca(ax)

    plt.scatter( xy[:,0], xy[:,1], c='gray')

    # map pvalues to sqrt(-1*log(pval/ranking_threshold))
    max_color_pval = 1e-10
    max_color_value = np.sqrt(-1*np.log10(max_color_pval/ranking_threshold))

    cmap = plt.get_cmap('bwr')
    #cmap = plt.get_cmap('seismic')
    #cmap = plt.get_cmap('coolwarm')

    fake_xys = {} # dictionary holding fake xy locations for results of cluster-based comparisons (not nbrhood-based)
    #                     clone_index will be -1 for these so there's no way to get xy
    for row in df.itertuples():
        _, pval, ii, gene = row[:4]
        if direction_column is None:
            dirn = 1
        else:
            dirn = 1 if row[4]>0 else -1
        assert gene not in exclude_strings
        pval = max(pval, 1e-299) # don't want problems in log10
        raw_color_value = np.sqrt( max(0, -1*np.log10( pval/ranking_threshold) ) )
        scaled_color_value = min(1.0, raw_color_value/max_color_value) # 0-->1
        signed_color_value = 0.5 + dirn * 0.5 * scaled_color_value
        text_color = cmap( signed_color_value)
        #print('plotting:', pval, ii, gene, dirn, raw_color_value, scaled_color_value, text_color, pval, row )

        if ii==-1:
            clp = row[5] if clp_index_for_minus1==0 else row[6]
            if clp not in fake_xys:
                fake_xys[clp] = make_n_pseudopoints(strings_per_cluster, xy[clusters==clp])
                assert fake_xys[clp].shape == (strings_per_cluster, 2)
            ii_xy = fake_xys[clp][0]
            fake_xys[clp] = fake_xys[clp][1:,:] # remove the point we just used
            gene_suffix='*'
        else:
            ii_xy = xy[ii]
            gene_suffix = ''

        if color_text:
            plt.text( ii_xy[0], ii_xy[1], gene+gene_suffix, fontweight='black',
                      bbox=dict(facecolor='white', alpha=0.5), color=text_color )
        else:
            plt.text( ii_xy[0], ii_xy[1], gene+gene_suffix, fontsize=8, fontweight=1+998*scaled_color_value,
                      bbox=dict(facecolor=text_color, alpha=0.5), color='black' )

    #plt.text( 0, 0 , 'hi there', fontsize=8, fontweight=1,
    #          bbox=dict(facecolor='white', alpha=0.5), color='black' )
    plt.xticks([],[])
    plt.yticks([],[])
    plt.xlabel('{}[0]'.format(xy_tag))# if xlabel is None else xlabel )
    plt.ylabel('{}[1]'.format(xy_tag))# if ylabel is None else ylabel )

    if title is not None:
        plt.title(title)

    if pngfile:
        print('making:', pngfile)
        plt.savefig(pngfile)

    help_message = f"""This plot summarizes the results of a graph
versus features analysis by labeling the clonotypes at the center of
each biased neighborhood with the name of the feature biased in that
neighborhood. The feature names are drawn in colored boxes whose
color is determined by the strength and direction of the feature score bias
(from bright red for features that are strongly elevated to bright blue
for features that are strongly decreased in the corresponding neighborhoods,
relative to the rest of the dataset).

At most one feature (the top scoring) is shown for each clonotype
(ie, neighborhood). The UMAP xy coordinates for this plot are
stored in adata.obsm['{xy_tag}']. The score used for ranking correlations
is '{ranking_column}'. The threshold score for displaying a feature is
{ranking_threshold}. The feature column is '{string_column}'. Since
we also run graph-vs-features using "neighbor" graphs that are defined
by clusters, ie where each clonotype is connected to all the other
clonotypes in the same cluster, some biased features may be associated with
a cluster rather than a specific clonotype. Those features are labeled with
a '*' at the end and shown near the centroid of the clonotypes belonging
to that cluster.
"""

    return help_message

def make_summary_figure(
        adata,
        outfile_prefix,
        markersize=20,
        pval_threshold_for_tcr_genes_results=0.05, # but bonferroni is harder on tcrs since there are so many lame genes
        pval_threshold_for_gex_scores_results=0.05,
):
    ''' Make a multi-panel summary figure showing clusters, graph-vs-graph
    conga scores, and graph-vs-feature hits
    '''

    figure_tag = GRAPH_VS_SUMMARY
    pngfile = f'{outfile_prefix}_{figure_tag}.png'

    util.setup_uns_dicts(adata) # shouldn't be necessary

    tcr_genes_results = adata.uns['conga_results'].get(
        TCR_GRAPH_VS_GEX_FEATURES, None)

    gex_scores_results = adata.uns['conga_results'].get(
        GEX_GRAPH_VS_TCR_FEATURES, None)

    if tcr_genes_results is None or gex_scores_results is None:
        print('ERROR make_summary_figure: adata is missing graph-vs-features',
              'results')

    conga_scores = np.array(adata.obs['conga_scores'])

    # make a nice combo fig
    nrows, ncols = 2,3
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*6,nrows*6))

    cbar_axs = []
    for irow, xy_tag in enumerate(['gex','tcr']):
        XY_TAG = xy_tag.upper()
        xy = adata.obsm['X_{}_2d'.format(xy_tag)]
        clusters = adata.obs['clusters_{}'.format(xy_tag)]

        #############################33
        ## first a plot colored by the clusters
        plt.sca(axs[irow,0])

        num_clusters = np.max(clusters)+1
        cluster_names_tag = 'clusters_{}_names'.format(xy_tag)
        if cluster_names_tag in adata.uns:
            cluster_names = adata.uns[cluster_names_tag]
        else:
            cluster_names = [ str(x) for x in range(num_clusters)]
        C = get_integers_color_dict( num_clusters)
        colors = [ C[x] for x in clusters]
        plt.scatter( xy[:,0], xy[:,1], s=markersize, c=colors )
        plt.xticks([],[])
        plt.yticks([],[])
        plt.title('{} clusters'.format(XY_TAG))
        plt.xlabel('{} UMAP1'.format(XY_TAG))
        plt.ylabel('{} UMAP2'.format(XY_TAG))
        #add_integers_legend( plt.gca(), C)
        add_categorical_legend(plt.gca(), cluster_names,
                               [C[x] for x in range(num_clusters)] )
        if irow==0:
            plt.text(0.01,0.01,'{} clonotypes'.format(adata.shape[0]),
                     ha='left', va='bottom', fontsize=8,
                     transform=plt.gca().transAxes)


        ############################
        ## now a plot colored by pval
        plt.sca(axs[irow,1])
        colors = np.maximum( -1*np.log10(conga_scores), 0.0 )
        sort_order = np.argsort( colors )

        vmin = math.sqrt( -1*math.log10( 1.0 ) ) # was 0.05 then 0.2
        vmax = math.sqrt( -1*math.log10( 1e-5 ) )

        plt.scatter( xy[sort_order,:][:,0], xy[sort_order,:][:,1], s=markersize,
                     c=np.sqrt( colors[ sort_order ] ), vmin=vmin, vmax=vmax )
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel('{} UMAP1'.format(XY_TAG))
        #plt.ylabel('{} UMAP2'.format(XY_TAG))
        #plt.title('kNN-overlap E-values (Eval=Pval*num_clones)')
        cbar = plt.colorbar(ticks=[math.sqrt(-1*math.log10(x))
                                   for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
        cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'])
        cbar_axs.append(cbar.ax)
        #xmin,xmax = plt.xlim()
        #ymin,ymax = plt.ylim()
        plt.title('CoNGA scores')



        ## now a plot of the nbrhood gene/score enrichments
        if xy_tag == 'tcr':
            plot_ranked_strings_on_cells(
                adata, tcr_genes_results, 'X_tcr_2d', 'clone_index',
                'mwu_pvalue_adj', pval_threshold_for_tcr_genes_results,
                'feature', ax=axs[irow,2])
            plt.sca(axs[irow,2])
            plt.title("Differential gene expression in TCR nbrhoods")
            plt.xlabel('{} UMAP1'.format(XY_TAG))
            plt.ylabel('')
        else:
            plot_ranked_strings_on_cells(
                adata, gex_scores_results, 'X_gex_2d', 'clone_index',
                'mwu_pvalue_adj', pval_threshold_for_gex_scores_results,
                'feature', direction_column='ttest_stat', ax=axs[irow,2])
            plt.title('TCR sequence feature bias in GEX nbrhoods')
            plt.xlabel('{} UMAP1'.format(XY_TAG))
            plt.ylabel('')

    #plt.tight_layout()
    plt.subplots_adjust(left=0.05, hspace=0.16, wspace=0.19, bottom=0.05,
                        right = 0.97, top=0.97)
    # fiddle with sizes
    for irow in range(2):
        ax = axs[irow,1]
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 1.17, box.height])
        ax = cbar_axs[irow]
        box = ax.get_position()
        ax.set_position([box.x0+0.04, box.y0, box.width, box.height])
    print('making:', pngfile)
    plt.savefig(pngfile)

    # store results and help message
    adata.uns['conga_results'][figure_tag] = pngfile
    help_message = """Summary figure for the graph-vs-graph and
    graph-vs-features analyses."""
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)


def make_tcr_db_match_plot(
        adata,
        outfile_prefix,
        pval_adj_threshold=1.,
):
    ''' Plot the location of the lit matches in the GEX and TCR UMAP landscapes
    should be called after calling conga.tcr_clumping.match_adata_tcrs_to_db_tcrs
    '''

    table_tag = TCR_DB_MATCH
    figure_tag = TCR_DB_MATCH_PLOT

    util.setup_uns_dicts(adata) # shouldn't be necessary

    results = adata.uns['conga_results'].get(table_tag, None)
    if table_tag is None:
        print('WARNING: conga.plotting.make_db_matches_figure::',
              "lit match results not found in adata.uns['conga_results']",
              'do you need to call conga.tcr_clumping.match_adata_tcrs_to_db_tcrs?')
        return

    if results.shape[0]==0:
        print('conga.plotting.make_db_matches_figure:: no significant hits')
        return

    results = results.sort_values('pvalue_adj').drop_duplicates('clone_index')

    missing_epitope = results.db_epitope==''
    results.loc[missing_epitope,'db_epitope'] = results.db_epitope_gene[missing_epitope]

    missing_epitope = results.db_epitope==''
    results.loc[missing_epitope, 'db_epitope'] = 'UNKNOWN'

    results['pmhc_label'] = (results.db_epitope+'_'+results.db_mhc_trim.astype(str))
    counts = results.pmhc_label.value_counts()

    pngfile = f'{outfile_prefix}_{figure_tag}.png'
    plt.figure(figsize=(12,6))

    for col, xytag in enumerate('gex tcr'.split()):
        plt.subplot(1,2,col+1)
        xy = adata.obsm[f'X_{xytag}_2d']

        # show all the points in light gray
        plt.scatter(xy[:,0], xy[:,1], c='#dddddd', s=5)

        # now show the lit matches colored
        if counts.shape[0] <=10:
            cmap = plt.get_cmap('tab10')
        else:
            cmap = plt.get_cmap('tab20')
        for ii, pmhc_label in enumerate(counts.index):
            inds = np.array(results[results.pmhc_label==pmhc_label].clone_index)
            print(ii, pmhc_label, len(inds))
            plt.scatter(xy[inds,0], xy[inds,1], c=[cmap.colors[ii%len(cmap.colors)]],
                        s=5, label=pmhc_label)
        plt.xticks([],[])
        plt.yticks([],[])
        plt.xlabel(f'{xytag.upper()} UMAP1')
        plt.ylabel(f'{xytag.upper()} UMAP2')
        plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(pngfile)
    print('made:', pngfile)

    # store results and help message
    adata.uns['conga_results'][figure_tag] = pngfile
    help_message = """GEX and TCR UMAPs showing the location of significant
    matches to the literature TCR database"""
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)


def make_clone_gex_umap_plots(
        adata,
        outfile_prefix,
        max_clones = 16,
        dpi=200,
        color_by_clusters=False,
):

    ''' This is called before we've condensed to a single cell per clone
    '''

    clusters_gex = np.array(adata.obs['clusters_gex'])

    tcrs = preprocess.retrieve_tcrs_from_adata(adata, include_subject_id_if_present=True)

    unique_tcrs = sorted(set(tcrs))
    # the clone ids are integer indices into the unique_tcrs list

    tcr2clone_id = { x:i for i,x in enumerate(unique_tcrs)}
    clone_ids = np.array( [ tcr2clone_id[x] for x in tcrs ])

    clone_counts = Counter(clone_ids)

    nrows = int(np.sqrt(max_clones))
    ncols = (max_clones-1)//nrows + 1

    plt.figure(figsize=(3*ncols,3*nrows))
    plotno=0

    xy = adata.obsm['X_gex_2d']
    for (clone_id, clone_size) in clone_counts.most_common(max_clones):
        plotno += 1
        plt.subplot(nrows, ncols, plotno)
        plt.scatter(xy[:,0], xy[:,1], s=16, c='gray',alpha=0.2)
        mask = clone_ids==clone_id
        colors = clusters_gex[mask] if color_by_clusters else 'blue'
        plt.scatter(xy[mask,0], xy[mask,1], s=16, c=colors, alpha=0.5)
        plt.xlabel('GEX UMAP1')
        plt.ylabel('GEX UMAP2')
        plt.xticks([], [])
        plt.yticks([], [])
        plt.text(0, 0, '{} cells'.format(clone_size), ha='left', va='bottom',
                 transform=plt.gca().transAxes)
    plt.tight_layout()
    pngfile = outfile_prefix+'_clone_gex_umaps.png'
    print('making:', pngfile)
    plt.savefig(pngfile, dpi=dpi)
    # now also save pdf
    plt.savefig(pngfile[:-3]+'pdf')



def make_clone_batch_clustermaps(
        adata,
        outfile_prefix,
        batch_keys = None, # if None will look in adata.uns
        max_clones = 75,
        min_clone_size = 5,
        dpi=200,
        conga_scores = None,             # np.array of pvals
        tcr_clumping_pvalues = None,      # --""--
        batch_bias_results = None, # tuple of dataframes: nbrhood_results, hotspot_results
        cmap_for_row_scores = 'viridis',
        show_mait_and_inkt_clones = True,
):

    '''
    '''
    try:
        import seaborn as sns
    except:
        print('ERROR seaborn is not installed')
        return

    if batch_keys is None:
        if 'batch_keys' not in adata.uns_keys():
            print('make_clone_batch_clustermaps: no batch_keys in adata.uns')
            return
        batch_keys = adata.uns['batch_keys']

    cmap_for_row_scores = plt.get_cmap(cmap_for_row_scores)

    num_clones = adata.shape[0]

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    organism = adata.uns['organism']

    # sort the clonotypes by clone size
    clone_sizes = np.array(adata.obs['clone_sizes'])
    sortl = sorted( [(x,i) for i,x in enumerate(clone_sizes)])
    sortl.reverse() # now in decreasing order of clone size

    top_clone_indices = [x[1] for x in sortl[:max_clones] if x[0] >= min_clone_size]

    if len(top_clone_indices)<2:
        print('make_clone_batch_clustermaps:: fewer than 2 clones meet size requirements')
        return

    for batch_key in batch_keys:
        batch_freqs = adata.obsm[batch_key].astype(float)/clone_sizes[:,np.newaxis]
        assert batch_freqs.shape[0] == num_clones
        num_choices = batch_freqs.shape[1]

        colorbar_row_colors = []

        def transform_adjusted_pvalues_into_unit_range(pvalues,
                                                       max_pval=1.0, min_pval=1e-5):
            vals = np.sqrt(np.maximum(0., -1*np.log10(np.maximum(
                1e-299, np.array(pvalues)/max_pval))))
            vmax = np.sqrt(-1*np.log10(min_pval/max_pval))
            retvals = np.maximum(0.001, np.minimum(0.999, vals/vmax))
            #assert False
            return retvals

        def get_row_colors_from_pvalues(pvalues, cmap=cmap_for_row_scores):
            vals =  transform_adjusted_pvalues_into_unit_range(pvalues)
            retvals = [ cmap(x) for x in vals]
            #assert False
            return retvals


        if conga_scores is not None: # these are adjusted pvalues
            colorbar_row_colors.append(get_row_colors_from_pvalues(conga_scores[top_clone_indices]))

        if tcr_clumping_pvalues is not None:
            colorbar_row_colors.append(get_row_colors_from_pvalues(tcr_clumping_pvalues[top_clone_indices]))

        if batch_bias_results is not None:
            nbrhood_results, hotspot_results = batch_bias_results
            tcr_nbrhood_pvals = np.full((num_clones,), num_clones).astype(float)
            for l in nbrhood_results.itertuples():
                if l.nbrs_tag == 'tcr' and l.batch_key == batch_key:
                    tcr_nbrhood_pvals[l.clone_index] = min(
                        l.pvalue_adj, tcr_nbrhood_pvals[l.clone_index])
            colorbar_row_colors.append(get_row_colors_from_pvalues(tcr_nbrhood_pvals[top_clone_indices]))

            hotspot_pval = 1.
            for l in hotspot_results.itertuples():
                if l.nbrs_tag == 'tcr' and l.batch_key == batch_key:
                    hotspot_pval = min(l.pvalue_adj, hotspot_pval)

            title = 'batch_key: {} num_choices: {} best_hotspot_pval: {:.1e}'\
                    .format( batch_key, num_choices, hotspot_pval)
        else:
            title = 'batch_key: {} num_choices: {}'\
                    .format( batch_key, num_choices)


        if show_mait_and_inkt_clones and organism in ['human','mouse']:
            colors = []
            for ii in top_clone_indices:
                tcr = tcrs[ii]
                if organism == 'human':
                    celltype = 2 * tcr_scoring.is_human_mait_alpha_chain(tcr[0])+\
                               tcr_scoring.is_human_inkt_tcr(tcr)
                elif organism == 'mouse':
                    celltype = 2 * tcr_scoring.is_mouse_mait_alpha_chain(tcr[0])+\
                               tcr_scoring.is_mouse_inkt_alpha_chain(tcr[0])
                colors.append(['white','C0','C1'][celltype])
            colorbar_row_colors.append(colors)

        # show clustermap of clone frequencies for the top clones
        A = []
        clone_labels = []
        def trim_allele(s):
            return s[:s.index('*')]
        def get_gene_strings(tcr, organism):
            if organism.endswith('_ig'):
                va_prefix = tcr[0][0][2]
            else:
                va_prefix = ''
            (va,ja,*_), (vb,jb,*_) = tcr
            return (va_prefix+trim_allele(va[4:]), trim_allele(ja[4:]),
                    trim_allele(vb[4:]), trim_allele(jb[4:]))
        gene_strings = [get_gene_strings(tcrs[x], organism) for x in top_clone_indices]
        max_valen = max(len(x[0]) for x in gene_strings)
        max_jalen = max(len(x[1]) for x in gene_strings)
        max_vblen = max(len(x[2]) for x in gene_strings)
        max_jblen = max(len(x[3]) for x in gene_strings)
        for ii, genes in zip(top_clone_indices, gene_strings):
            A.append(batch_freqs[ii,:])
            va,ja,vb,jb = genes
            cdr3a, cdr3b = tcrs[ii][0][2], tcrs[ii][1][2]

            clone_labels.append('{:3d}  {:{}s} {:{}s} {:18s} {:{}s} {:{}s} {:18s}'\
                                .format(clone_sizes[ii],
                                        va, max_valen, ja, max_jalen, cdr3a,
                                        vb, max_vblen, jb, max_jblen, cdr3b))

        A = np.array(A)

        ## things we might want to add row colors for:
        ##
        ## - tcr nbr batch bias for this batch_key
        ## - conga score
        ## - tcr clumping score
        ##

        fig_width, fig_height = (10,14)
        if not colorbar_row_colors:
            colorbar_row_colors = None
        cm = sns.clustermap(A, metric='euclidean', col_cluster=False, vmin=0, vmax=1,
                            figsize=(fig_width,fig_height), cmap='Reds',
                            row_colors=colorbar_row_colors)

        row_indices = cm.dendrogram_row.reordered_ind
        row_labels = [clone_labels[x] for x in row_indices]
        cm.ax_heatmap.set_yticks( 0.5+np.arange(A.shape[0]))
        cm.ax_heatmap.set_yticklabels( row_labels, rotation='horizontal', fontdict={'fontsize':8, 'fontfamily':'monospace'} )
        cm.ax_heatmap.set_xticks([])
        cm.ax_heatmap.set_xticklabels([])
        #cm.fig.suptitle(title)
        cm.ax_heatmap.set_title(title)


        pngfile = f'{outfile_prefix}_clone_clustermap_batch_{batch_key}.png'
        print('making:', pngfile)
        cm.savefig(pngfile, dpi=200)
        #assert False




def make_feature_panel_plots(
        adata,
        xy_tag,
        all_nbrs,
        results_df,
        pngfile,
        feature_type,
        max_panels_per_bicluster=3,
        max_pvalue=0.05,
        nrows=8, #6,
        ncols=5, #4,
        panel_width_inches=2.5,
        use_nbr_frac=None,
        title=None,
        cmap='viridis',
        sort_order=True, # plot points sorted by feature value
        point_size=10,
):
    ''' Assumes results_df has column mwu_pvalue_adj for sorting

    returns a help message string.
    '''
    assert feature_type in ['gex','tcr']

    if results_df.empty:
        return

    xy = adata.obsm['X_{}_2d'.format(xy_tag)]

    # first let's figure out how many panels we could have
    # sort the results by pvalue
    df = results_df.sort_values('mwu_pvalue_adj')# makes a copy

    nbr_frac_for_cluster_results = (max(all_nbrs.keys()) if use_nbr_frac is None
                                    else use_nbr_frac)

    clp_counts = Counter()
    seen = set()
    inds=[]
    for row in df.itertuples():
        if row.mwu_pvalue_adj > max_pvalue:
            break
        feature = row.feature
        if feature in seen:
            continue
        clp = (row.gex_cluster, row.tcr_cluster)
        if clp_counts[clp] >= max_panels_per_bicluster:
            continue
        seen.add(feature)
        clp_counts[clp] += 1
        inds.append(row.Index)
    if not inds:
        print('no results to plot!')
        return

    if len(inds) < nrows*ncols:
        # make fewer panels
        nrows = max(1, int(np.sqrt(len(inds))))
        ncols = (len(inds)-1)//nrows + 1

    df = df.loc[inds[:nrows*ncols],:]
    assert df.shape[0] <= nrows*ncols

    # try to get the raw values
    feature_to_raw_values = {}

    var_names = list(adata.raw.var_names)
    features = set(df['feature'])
    for f in features:
        if feature_type=='gex' and f in var_names:
            # careful since TR/IG gene names are present in var_names but its
            # better to use the VDJ information
            feature_to_raw_values[f] = adata.raw.X[:, var_names.index(f)].toarray()[:,0]
        elif f in adata.obs:
            feature_to_raw_values[f] = np.array(adata.obs[f])
            if f=='clone_sizes':
                feature_to_raw_values[f] = np.log1p(feature_to_raw_values[f])
        elif f=='nndists_gex_rank':
            if 'nndists_gex' in adata.obs_keys():
                nndists_gex = np.array(adata.obs['nndists_gex'])
            else:
                print('WARNING nndists_gex not in adata.obs!')
                nndists_gex = np.zeros(num_clones)
            feature_to_raw_values[f] = np.log1p(np.argsort(-1*nndists_gex))
        else:
            assert feature_type=='tcr'
            feature_score_table = tcr_scoring.make_tcr_score_table(adata, [f])
            feature_to_raw_values[f] = feature_score_table[:,0]


    figsize= (ncols*panel_width_inches, nrows*panel_width_inches)
    plt.figure(figsize=figsize)
    plotno=0
    for row in df.itertuples():
        plotno+=1
        plt.subplot(nrows, ncols, plotno)

        feature = row.feature

        scores = np.array(feature_to_raw_values[feature])

        if sort_order:
            reorder = np.argsort(scores)
        else:
            reorder = np.arange(len(scores))

        row_nbr_frac = use_nbr_frac if use_nbr_frac is not None else \
                       nbr_frac_for_cluster_results if row.nbr_frac==0.0 else \
                       row.nbr_frac
        nbrs = all_nbrs[row_nbr_frac][0] if xy_tag=='gex' else \
               all_nbrs[row_nbr_frac][1]
        assert nbrs.shape[0] == adata.shape[0]
         # this will not work for ragged nbr arrays right now
        num_neighbors = nbrs.shape[1]

        nbr_averaged_scores = (
            scores + scores[nbrs].sum(axis=1))/(num_neighbors+1)

        plt.scatter(xy[reorder,0], xy[reorder,1],
                    c=nbr_averaged_scores[reorder], cmap=cmap,
                    s=point_size)
        # plt.scatter(xy[:,0], xy[:,1], c=nbr_averaged_scores, cmap='coolwarm',
        #             s=15)
        plt.title('{} ({:d},{:d}) {:.1e}'.format(feature,
                                                 int(row.gex_cluster+.1),
                                                 int(row.tcr_cluster+.1),
                                                 row.mwu_pvalue_adj))
        plt.xticks([],[])
        plt.yticks([],[])
        plt.text(0.99, 0.01, f'K={num_neighbors} {xy_tag.upper()} nbr-avged',
                 va='bottom', ha='right', fontsize=9,
                 transform=plt.gca().transAxes)
        plt.text(0.01, 0.99, f'{np.min(nbr_averaged_scores):.2f}',
                 va='top', ha='left', fontsize=8,
                 transform=plt.gca().transAxes)
        plt.text(0.99, 0.99, f'{np.max(nbr_averaged_scores):.2f}',
                 va='top', ha='right', fontsize=8,
                 transform=plt.gca().transAxes)
        if (plotno-1)//ncols == nrows-1:
            plt.xlabel('{} UMAP1'.format(xy_tag.upper()))
        if plotno%ncols == 1:
            plt.ylabel('{} UMAP2'.format(xy_tag.upper()))

    tl_rect = [0,0,1,1]
    if title is not None:
        title_height_inches = 0.25
        plt.suptitle(title)
        tl_rect = [0, 0, 1, 1.0-title_height_inches/figsize[1]]
    plt.tight_layout(rect=tl_rect)
    print('making:', pngfile)
    plt.savefig(pngfile, dpi = 300)

    help_message = f"""Graph-versus-feature analysis was used to identify
a set of {feature_type.upper()} features that showed biased distributions
in {xy_tag.upper()} neighborhoods. This plot shows the distribution of the
top-scoring {feature_type.upper()} features on the {xy_tag.upper()}
UMAP 2D landscape. The features are ranked by 'mwu_pvalue_adj' ie
Mann-Whitney-Wilcoxon adjusted P value (raw P value * number of comparisons).
At most {max_panels_per_bicluster} features from clonotype neighbhorhoods
in each (GEX,TCR) cluster pair are shown. The raw scores for each feature
are averaged over the K nearest neighbors (K is indicated in the lower
right corner of each panel) for each clonotype. The min and max
nbr-averaged scores are shown in the upper corners of each panel.
"""
    if sort_order:
        help_message += "Points are plotted in order of increasing feature score.\n"

    return help_message


def get_raw_feature_scores( feature, adata, feature_type):
    if feature.startswith('gex_cluster'):
        return (np.array(adata.obs['clusters_gex'])==int(feature[11:])).astype(float)
    elif feature.startswith('tcr_cluster'):
        return (np.array(adata.obs['clusters_tcr'])==int(feature[11:])).astype(float)
    elif feature == 'nndists_gex_rank':
        if 'nndists_gex' in adata.obs_keys():
            nndists_gex = np.array(adata.obs['nndists_gex'])
        else:
            print('WARNING nndists_gex not in adata.obs!')
            nndists_gex = np.zeros(adata.shape[0])
        return np.log1p(np.argsort(-1*nndists_gex))
    elif feature == 'clone_sizes': # note that we are taking log1p for compatibility with the code in correlations.py
        return np.log1p(np.array(adata.obs['clone_sizes']))
    elif (feature_type == 'gex' or
          (feature_type==None and feature in adata.raw.var_names)):
        # will this be slow? creating list every time...
        #  probably OK for plotting routines
        return adata.raw.X[:, list(adata.raw.var_names).index(feature)].toarray()[:,0]
    else:
        assert feature_type in [None, 'tcr']
        return tcr_scoring.make_tcr_score_table(
            adata,[feature])[:,0].astype(float)

def make_raw_feature_scores_table(
        features,
        feature_types,
        adata,
        verbose=False,
):
    ''' Make the whole table at once, should be faster for lots of features

    returns numpy array of shape (num_clones, len(features))

    '''
    assert all(x in ['gex','tcr'] for x in feature_types)

    scoretable = np.zeros((adata.shape[0], len(features)))

    # group the features: gene, tcr, other aka 'special'
    is_special_feature = np.array(
        [('_cluster' in x or x=='clone_sizes' or x=='nndists_gex_rank')
         for x in features])
    var_name2index = {x:i for i,x in enumerate(adata.raw.var_names)}
    is_gene_feature = np.array(
        [x in var_name2index and y=='gex'
         for x,y in zip(features, feature_types)])
    is_tcr_feature = ~(is_special_feature | is_gene_feature)

    # get the genes
    # list of (index in features, index in adata.raw.var_names)
    gene_index_pairs = [(i,var_name2index[x])
                        for i,(x,y) in enumerate(zip(features, is_gene_feature))
                        if y]
    inds0 = [x[0] for x in gene_index_pairs]
    inds1 = [x[1] for x in gene_index_pairs]
    gene_scores = adata.raw.X[:,inds1].toarray()
    scoretable[:, inds0] = gene_scores

    if np.sum(is_tcr_feature): # get the tcr features
        ind_feature_pairs = [
            (i,x) for i,(x,y) in enumerate(zip(features, is_tcr_feature)) if y]
        inds = [x[0] for x in ind_feature_pairs]
        tcr_features = [x[1] for x in ind_feature_pairs]
        tcr_scores = tcr_scoring.make_tcr_score_table(
            adata, tcr_features, verbose=verbose)
        #tcr_scores = np.zeros((adata.shape[0], len(tcr_features)))
        scoretable[:, inds] = tcr_scores

    # fill in the 'special' features
    for ind,(feature,is_special) in enumerate(zip(features,is_special_feature)):
        if is_special:
            if feature == 'nndists_gex_rank':
                if 'nndists_gex' in adata.obs_keys():
                    nndists_gex = np.array(adata.obs['nndists_gex'])
                else:
                    print('WARNING nndists_gex not in adata.obs!')
                    nndists_gex = np.zeros(adata.shape[0])
                scoretable[:,ind] = np.log1p(np.argsort(-1*nndists_gex))
            elif feature.startswith('gex_cluster'):
                num = int(feature[11:]) # cluster number
                scoretable[:,ind] = (
                    np.array(adata.obs['clusters_gex']==num)).astype(float)
            elif feature.startswith('tcr_cluster'):
                num = int(feature[11:]) # cluster number
                scoretable[:,ind] = (
                    np.array(adata.obs['clusters_tcr']==num)).astype(float)
            else:
                assert feature == 'clone_sizes'
                scoretable[:,ind] = np.log1p(np.array(adata.obs['clone_sizes']))

    return scoretable

def plot_hotspot_umap(
        adata,
        xy_tag,
        results_df,
        pngfile,
        nbrs = None,
        compute_nbr_averages=True,
        nrows=8,
        ncols=5,
        panel_width_inches=2.5,
        title=None,
        cmap='viridis',
        sort_order=True, #plot points ordered by feature score
        point_size=10,
        max_feature_correlation=0.9, # for filtering redundancy
):
    """
    xy_tag: use 'gex' or 'tcr' to set umap space used for plotting the hotspot features
    results_df : pandas df of hotspot features to plot. Expected columns are pvalue_adj feature feature_type
    where feature_type is either 'gex' or 'tcr'. We need that since some feature strings can be both. Output
    from correlations.find_hotspots can be fed in directly

    returns a help_message string
    """
    if results_df.shape[0]==0:
        print('no results to plot:', pngfile)
        return

    df = results_df.sort_values('Z', ascending=False)\
                   .sort_values('pvalue_adj') # ensure sorted
    #df = df.iloc[:nrows*ncols,:]

    xy = adata.obsm['X_{}_2d'.format(xy_tag)]
    var_names = list(adata.raw.var_names)
    if df.shape[0] < nrows*ncols:
        # make fewer panels
        nrows = max(1, int(np.sqrt(df.shape[0])))
        ncols = (df.shape[0]-1)//nrows + 1

    figsize=(ncols*panel_width_inches, nrows*panel_width_inches)
    plt.figure(figsize=figsize)
    plotno=0
    other_scores = []
    for row in df.itertuples():

        scores = get_raw_feature_scores(row.feature, adata, row.feature_type)
        if compute_nbr_averages:
            assert nbrs.shape[0] == adata.shape[0]
            # this will not work for ragged nbr arrays
            #  (but we could change it to work)
            num_neighbors = nbrs.shape[1]
            scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

        # check if we are too close to a feature we already plotted
        too_close = False
        for ii, oscores in enumerate(other_scores):
            corr = 1-distance.correlation(scores,oscores)
            if corr > max_feature_correlation:
                #print(corr, row.feature, df.iloc[ii].feature)
                too_close = True
                break
        if too_close:
            #print('skipping, too close')
            continue

        #
        plotno+=1
        if plotno>nrows*ncols:
            break
        plt.subplot(nrows, ncols, plotno)

        other_scores.append(scores)

        if sort_order:
            reorder = np.argsort(scores)
        else:
            reorder = np.arange(len(scores))
        plt.scatter(xy[reorder,0], xy[reorder,1], c=scores[reorder],
                    cmap=cmap, s=point_size)
        plt.title('{} {:.1e}'.format(row.feature, row.pvalue_adj))
        plt.xticks([],[])
        plt.yticks([],[])
        if compute_nbr_averages:
            plt.text(0.99, 0.01,
                     f'K={num_neighbors} {xy_tag.upper()} nbr-avged',
                     va='bottom', ha='right', fontsize=9,
                     transform=plt.gca().transAxes)
        plt.text(0.01, 0.99, f'{np.min(scores):.2f}',
                 va='top', ha='left', fontsize=8,
                 transform=plt.gca().transAxes)
        plt.text(0.99, 0.99, f'{np.max(scores):.2f}',
                 va='top', ha='right', fontsize=8,
                 transform=plt.gca().transAxes)
        if (plotno-1)//ncols == nrows-1:
            plt.xlabel('{} UMAP1'.format(xy_tag.upper()))
        if plotno%ncols == 1:
            plt.ylabel('{} UMAP2'.format(xy_tag.upper()))

    tl_rect = [0,0,1,1]
    if title is not None:
        title_height_inches = 0.25
        plt.suptitle(title)
        tl_rect = [0, 0, 1, 1.0-title_height_inches/figsize[1]]

    plt.tight_layout(rect=tl_rect)
    print('making:', pngfile)
    plt.savefig(pngfile)

    help_message = f"""HotSpot analysis (Nir Yosef lab, PMID: 33951459)
was used to identify a set of GEX (TCR) features that showed biased
distributions in TCR (GEX) space. This plot shows the distribution of the
top-scoring HotSpot features on the {xy_tag.upper()}
UMAP 2D landscape. The features are ranked by adjusted P value
(raw P value * number of comparisons). The raw scores for each feature
are averaged over the K nearest neighbors (K is indicated in the lower
right corner of each panel) for each clonotype. The min and max
nbr-averaged scores are shown in the upper corners of each panel.

Features are filtered based on correlation coefficient to reduce
redundancy: if a feature has a correlation of >= {max_feature_correlation}
(the max_feature_correlation argument to conga.plotting.plot_hotspot_umap)
to a previously plotted feature, that feature is skipped.
"""
    if sort_order:
        help_message += "Points are plotted in order of increasing feature score\n"
    return help_message


def filter_sorted_features_by_correlation(
        features,
        #score_matrix,
        correlation_matrix,
        max_clusters,
):
    ''' return the new features, and a dictionary mapping from each
    representative feature to any other cluster members

    assume that features are sorted in decreasing order of "interest"
    will choose

    '''

    num_features = len(features)
    assert correlation_matrix.shape == (num_features, num_features)
    assert np.min(correlation_matrix)>=-1.01
    assert np.max(correlation_matrix)<= 1.01

    if num_features <= max_clusters:
        return list(features), {}

    Y = distance.squareform(1-correlation_matrix, force='tovector')
    Z = hierarchy.linkage(Y, method='ward')
    clusters = hierarchy.fcluster(Z, max_clusters, criterion = 'maxclust')
    assert len(clusters) == num_features
    if len(set(clusters)) > max_clusters:
        print('WHOAH ERROR hierarchy.fcluster failed: desired=', max_clusters,
              '< actual=', len(set(clusters)))

    reps = []
    duplicates = {}
    for c in set(clusters):
        members = sorted(np.nonzero(clusters==c)[0])
        rep = members[0]
        f = features[rep]
        reps.append(rep)
        if len(members)>1:
            duplicates[f] = [features[x] for x in members[1:]]

    # preserve the sort by importance:
    reps.sort()
    new_features = [features[x] for x in reps]
    return new_features, duplicates






def plot_interesting_features_vs_clustermap(
        adata,
        features, # assume that features are sorted by decreasing "interest"
        feature_types, # listlike of either 'gex','tcr'
        pngfile,
        dist_tag, # gex or tcr
        nbrs = None,
        compute_nbr_averages=False,
        feature_labels = None,
        feature_scores = None,
        feature_scores_cmap = 'viridis',
        max_type_features = 50, # in the clustermap
        show_gex_cluster_colorbar = None,
        show_tcr_cluster_colorbar = None,
        show_VJ_gene_segment_colorbars = None,
        show_batch_keys_colorbars = True,
        min_vmax = 0.45,
        max_vmax = 0.45,
        # max_redundant_features = None,
        # redundancy_threshold = 0.9, # correlation
        extra_feature_colors = None,
        rescale_factor_for_self_features = 0.33,
        use_1d_landscape_for_cell_order = False,
        title=None,
        verbose=False,
        dpi=200, # of the final image, which is ~12 inches wide
):
    ''' Makes a seaborn clustermap: cols are cells, sorted by hierarchical clustering wrt X_pca_{dist_tag}
    rows are the genes, ordered by clustermap correlation
    '''
    try:
        import seaborn as sns
    except:
        print('ERROR seaborn is not installed')
        return

    try:
        import fastcluster
    except:
        print('fastcluster is not available. Consider installing for faster performance.')

    # convert to numpy arrays
    features = np.array(features)
    feature_types = np.array(feature_types)

    # for re-indexing feature_types, feature_scores, etc after we filter
    #original_features = list(features)

    assert dist_tag in ['gex','tcr']
    assert all(x in ['gex','tcr'] for x in feature_types)

    if show_gex_cluster_colorbar is None:
        show_gex_cluster_colorbar = dist_tag=='gex'

    if show_tcr_cluster_colorbar is None:
        show_tcr_cluster_colorbar = dist_tag=='tcr'

    if show_VJ_gene_segment_colorbars is None:
        show_VJ_gene_segment_colorbars = dist_tag=='tcr'

    if show_batch_keys_colorbars and 'batch_keys' not in adata.uns:
        show_batch_keys_colorbars = False

    if feature_labels is None:
        feature_labels = features[:]

    var_names = list(adata.raw.var_names)
    #assert len(gene_labels) == len(genes)
    #gene_labels = [ y for x,y in zip(genes,gene_labels) if x in var_names]
    #genes = [ x for x in genes if x in var_names]

    if len(features)<2:
        print('too few features for clustermap:', len(features))
        return

    nrows = len(features)
    ncols = adata.shape[0]

    if use_1d_landscape_for_cell_order:
        landscape_tag = f'X_{dist_tag}_1d'
        if landscape_tag not in adata.obsm_keys():
            print('ERROR plot_interesting_features_vs_clustermap:',
                  'use_1d_landscape_for_cell_order=True but',landscape_tag,
                  'missing from adata.obsm_keys()')
            return
        X_1d = adata.obsm[landscape_tag][:,0]
        cells_order = np.argsort(X_1d)
        cells_linkage = None
        col_cluster = False
    else:
        # compute linkage of cells based on tcr
        X = adata.obsm['X_pca_'+dist_tag] # the kernal pca components
        print(f'computing pairwise X_pca_{dist_tag} distances, {X.shape}')
        Y = distance.pdist(X, metric='euclidean')

        print(f'computing linkage matrix from pairwise X_pca_{dist_tag}',
              'distances')
        cells_order = None
        cells_linkage = hierarchy.linkage(Y, method='ward')
        col_cluster = True

    A = make_raw_feature_scores_table(
        features, feature_types, adata, verbose=verbose).T

    A_mn = np.mean(A, axis=1)
    A_std = np.maximum(np.std(A, axis=1), 1e-9) # no divide by 0
    A = (A-A_mn[:,None])/A_std[:,None]

    if compute_nbr_averages:
        num_neighbors = nbrs.shape[1]
        if verbose:
            print('nbr-avging A', num_neighbors, len(features))
        last_time = time.time()
        Asum = A.copy()
        for i,ii_nbrs in enumerate(nbrs):
            if verbose and i and i%5000==0:
                print(i, time.time()-last_time)
                last_time = time.time()
            for j in ii_nbrs:
                Asum[:,i] += A[:,j]
        A = Asum / (num_neighbors+1)

        # lighten the colors of the features of the same type as the
        # landscape and nbrs, since they will be correlated across neighborhoods
        # and not comparable to the features of the other type (gex or tcr)
        is_dist_feature_type = np.array([x==dist_tag for x in feature_types])
        if verbose:
            print('rescale scores for:', np.sum(is_dist_feature_type),
                  'features by', rescale_factor_for_self_features)
        A[is_dist_feature_type,:] *= rescale_factor_for_self_features



    tiny_lines = None
    feature_mask = np.full((len(features),), True)

    feature_nbrs = {}
    for ftype in ['gex','tcr']:
        # are there too many features of this type?
        ftype_mask = feature_types==ftype
        if np.sum(ftype_mask) > max_type_features:
            # have to eliminate some
            A_ftype = A[ftype_mask,:].copy()
            C_ftype = 1-distance.squareform(
                distance.pdist(A_ftype, metric='correlation'), force='tomatrix')
            new_features, duplicates = filter_sorted_features_by_correlation(
                features[ftype_mask], C_ftype, max_type_features)
            feature_nbrs.update(duplicates)
            # update feature_nbrs based on duplicates
            # for rep, dups in duplicates.items():
            #     newdups = list(dups)
            #     for dup in dups:
            #         newdups += feature_nbrs.get(dup,[])
            #     feature_nbrs[rep] = feature_nbrs.get(rep,[])+newdups
            new_features = set(new_features) # fast membership checking
            for i,(f,m) in enumerate(zip(features, ftype_mask)):
                if m and (f not in new_features):
                    feature_mask[i] = False # excluded

    if np.sum(feature_mask)<len(features): # had to exclude some
        # write out the feature neighbors for inspection
        dfl = []
        for f, fnbrs in feature_nbrs.items():
            dfl.append(OrderedDict(
                feature=f,
                redundant_neighbors=','.join(fnbrs),
                num_redundant_neighbors=len(fnbrs)))
        df = pd.DataFrame(dfl)
        df.sort_values('num_redundant_neighbors', inplace=True,
                       ascending=False)
        df.to_csv(pngfile+'_filtered_feature_info.tsv', sep='\t',
                  index=False)
        tiny_lines = []
        for l in df.itertuples():
            tiny_lines.extend(
                [x+'\n' for x in ([l.feature+':']+
                                  l.redundant_neighbors.split(',')+
                                  [''])])

        #
        A = A[feature_mask,:]
        old_features = list(features)
        features = np.array(features)[feature_mask]
        feature_types = np.array(feature_types)[feature_mask]
        new_feature_labels = []
        for feature, old_label, m in zip(old_features,
                                         feature_labels,
                                         feature_mask):
            if m:
                nbr_count = len(feature_nbrs.get(feature,[]))
                if nbr_count>0:
                    new_feature_labels.append(
                        '{} [+{}]'.format(old_label, nbr_count))
                else:
                    new_feature_labels.append(old_label)
        feature_labels = new_feature_labels[:]
        nrows = len(features)
        if feature_scores is not None:
            feature_scores = np.array(feature_scores)[feature_mask]
        if extra_feature_colors is not None:
            extra_feature_colors = np.array(extra_feature_colors)[feature_mask]

    ## add some cell colors
    organism = adata.uns['organism']
    colorbar_names, colorbar_colors, colorbar_sorted_tuples = [], [], []

    if show_batch_keys_colorbars:
        for batch_key in adata.uns['batch_keys']:
            batch_counts = adata.obsm[batch_key]
            num_batches = batch_counts.shape[1]
            batch_color_dict = get_integers_color_dict(num_batches)
            top_batch = np.argmax(batch_counts, axis=1)
            colorbar_names.append(batch_key)
            colorbar_colors.append( [batch_color_dict[x] for x in top_batch])
            colorbar_sorted_tuples.append(
                [(str(x), batch_color_dict[x])
                 for x in range(num_batches)])



    if show_gex_cluster_colorbar:
        clusters_gex = np.array(adata.obs['clusters_gex'])
        num_clusters_gex = np.max(clusters_gex)+1
        cluster_color_dict = get_integers_color_dict(num_clusters_gex)
        colorbar_names.append( 'gex_cluster')
        colorbar_colors.append( [cluster_color_dict[x] for x in clusters_gex] )
        colorbar_sorted_tuples.append(
            [ ('gexclus'+str(x), cluster_color_dict[x])
              for x in range(num_clusters_gex)])

    if show_tcr_cluster_colorbar:
        clusters_tcr = np.array(adata.obs['clusters_tcr'])
        num_clusters_tcr = np.max(clusters_tcr)+1
        cluster_color_dict = get_integers_color_dict(num_clusters_tcr)
        colorbar_names.append( 'tcr_cluster')
        colorbar_colors.append( [cluster_color_dict[x] for x in clusters_tcr] )
        colorbar_sorted_tuples.append(
            [ ('tcrclus'+str(x), cluster_color_dict[x])
              for x in range(num_clusters_tcr)])

    if show_VJ_gene_segment_colorbars:
        tcrs = preprocess.retrieve_tcrs_from_adata(adata)
        gene_colors, sorted_gene_colors = assign_colors_to_conga_tcrs(
            tcrs, organism, return_sorted_color_tuples=True)
        colorbar_names.extend('VA JA VB JB'.split())
        colorbar_colors.extend(gene_colors)
        num_to_show = 10 # in the text written at the top
        colorbar_sorted_tuples.extend(
            [ x[:num_to_show] for x in sorted_gene_colors ] )

    # add some row colors
    colorbar_row_colors = []
    if extra_feature_colors is not None:
        colorbar_row_colors.append(extra_feature_colors)

    if feature_scores is not None:
        cm = plt.get_cmap(feature_scores_cmap)
        vmin, vmax = min(feature_scores), max(feature_scores)
        if vmax==vmin: vmax +=.001
        normed_scores = [ (x-vmin)/(vmax-vmin) for x in feature_scores]
        colorbar_row_colors.append( [ cm(x) for x in normed_scores])

    type2color = {'tcr':'C0', 'gex':'C1'}
    colorbar_row_colors.append(
        [type2color.get(x,'black') for x in feature_types] )

    # try to only use other-type features for setting max
    mask = [ x!=dist_tag for x in feature_types]
    if np.sum(mask):
        mx = np.max(np.abs(A[mask,:]))
    else:
        mx = np.max(np.abs(A))

    mx = min(max_vmax, max(min_vmax, mx))


    # now reorder using cells_order
    # will do nothing if use_1d_landscape_for_cell_order==False
    #
    if cells_order is not None:
        A = A[:, cells_order]
        colorbar_colors = [[x[i] for i in cells_order] for x in colorbar_colors]

    xmargin = 0.1
    ymargin_top = 0.4
    ymargin_bottom = 0.1
    col_dendro_height = max(1.5 if col_cluster else 0.1,
                            len(colorbar_colors)*0.1)
    col_colors_height = len(colorbar_colors)*0.2
    heatmap_height = len(features)*0.15
    row_dendro_width = 2.5 #inches
    row_colors_width = (0.1 if colorbar_row_colors is None else
                        len(colorbar_row_colors)*0.3)
    heatmap_width = 7.5
    fig_width = 2*xmargin + row_dendro_width + row_colors_width + heatmap_width
    fig_height = (ymargin_top + col_dendro_height + col_colors_height +
                  heatmap_height + ymargin_bottom)

    print('making clustermap; num_features=', len(features))

    x0 = xmargin/fig_width
    x1 = (xmargin+row_dendro_width)/fig_width
    x2 = (xmargin+row_dendro_width+row_colors_width)/fig_width
    x3 = (fig_width-xmargin)/fig_width

    y0 = ymargin_bottom/fig_height
    y1 = (ymargin_bottom+heatmap_height)/fig_height
    y2 = (ymargin_bottom+heatmap_height+col_colors_height)/fig_height
    y3 = (fig_height-ymargin_top)/fig_height

    # not present in older version of seaborn:
    #dendrogram_inches = 2.
    #dendrogram_ratio = (dendrogram_inches/fig_height, dendrogram_inches/fig_width)
    try:
        cm = sns.clustermap(
            A, col_linkage=cells_linkage, metric='correlation', cmap='coolwarm',
            col_cluster=col_cluster, vmin=-mx, vmax=mx,
            col_colors=colorbar_colors, row_colors=colorbar_row_colors,
            figsize=(fig_width,fig_height))#, dendrogram_ratio=dendrogram_ratio)
    except:
        tmpfile = pngfile+'_bad_matrix.txt'
        print('ERROR plot_interesting_features_vs_clustermap',
              'seaborn clustermap failed, writing matrix to file',
              tmpfile)
        for ii, f in enumerate(features):
            print(ii, f, np.mean(A[ii,:]), np.std(A[ii,:]))
        np.savetxt(tmpfile, A)
        return

    cm.ax_row_dendrogram.set_position([x0,y0,x1-x0,y1-y0])
    if cm.ax_row_colors is not None:
        cm.ax_row_colors.set_position([x1,y0,x2-x1,y1-y0])
    cm.ax_heatmap.set_position([x2,y0,x3-x2,y1-y0])
    if cm.ax_col_colors is not None:
        cm.ax_col_colors.set_position([x2,y1,x3-x2,y2-y1])
    if cm.ax_col_dendrogram:
        cm.ax_col_dendrogram.set_position([x2,y2,x3-x2,y3-y2])
    cax_height = ymargin_top + col_dendro_height + col_colors_height - 0.5
    cm.cax.set_position([x0,(fig_height-0.1-cax_height)/fig_height,
                         0.3/fig_width,cax_height/fig_height])

    row_indices = cm.dendrogram_row.reordered_ind
    row_labels = [feature_labels[x] for x in row_indices]
    cm.ax_heatmap.set_yticks( 0.5+np.arange(nrows))
    cm.ax_heatmap.set_yticklabels( row_labels, rotation='horizontal', fontdict={'fontsize':8} )
    cm.ax_heatmap.set_xticks([])
    cm.ax_heatmap.set_xticklabels([])

    # show the top couple genes and their colors
    ax = cm.ax_col_dendrogram
    col_dendro_width= fig_width*(ax.get_position().x1 - ax.get_position().x0 )
    col_dendro_height= fig_height*(ax.get_position().y1 - ax.get_position().y0)
    #offset = 0.05 # inches
    if title is not None:
        #ax.text(0.5, 1.0-offset/col_dendro_height, title,
        ax.text(0.5, 1.02, title,
                va='bottom', ha='center', fontsize=9, color='k',
                transform=ax.transAxes)
        #offset += 0.2

    line_height_inches = 0.09
    fontsize=7
    for ii, sorted_tuples in enumerate(colorbar_sorted_tuples):
        #y = 1.0 - ii/len(colorbar_sorted_tuples)
        y = 1.0 - (ii*line_height_inches)/col_dendro_height
        ax.text(-0.01, y,
                colorbar_names[ii], va='top', ha='right', fontsize=fontsize,
                color='k', transform=ax.transAxes)
        for jj in range(len(sorted_tuples)):
            gene, color = sorted_tuples[jj]
            ax.text(float(jj)/len(sorted_tuples),
                    y,#1.0 - (ii*line_height_inches+offset)/col_dendro_height,
                    gene, va='top', ha='left', fontsize=fontsize, color=color,
                    transform=ax.transAxes)

    # label the column color bars
    if colorbar_colors and cm.ax_col_colors:
        ax = cm.ax_col_colors
        for ii, name in enumerate(colorbar_names):
            ax.text(-0.01, 1.0 - (ii+0.5)/len(colorbar_names), name,
                    va='center', ha='right', fontsize=fontsize,
                    color='k', transform=ax.transAxes)

    # show the tiny_lines
    if tiny_lines is not None:
        ax = cm.ax_row_dendrogram
        row_dendro_width = fig_width*(ax.get_position().x1-ax.get_position().x0)
        row_dendro_height=fig_height*(ax.get_position().y1-ax.get_position().y0)

        max_lines = int(row_dendro_height / 0.08)
        col_width_frac = 0.33 / row_dendro_width
        x_offset = 0
        tiny_fontsize=5
        ax.text(0.0, 1.0, 'Correlated\nfeatures\nomitted:', va='bottom',
                ha='left', fontsize=tiny_fontsize,
                transform=ax.transAxes, zorder=-1, color='black')
        while tiny_lines:
            ax.text(x_offset, 1.0, ''.join(tiny_lines[:max_lines]), va='top',
                    ha='left', fontsize=tiny_fontsize,
                    transform=ax.transAxes, zorder=-1, color='blue')
            x_offset += col_width_frac
            tiny_lines = tiny_lines[max_lines:]
            if tiny_lines and x_offset+col_width_frac/2>1.0:
                ax.text(x_offset, 1.0, '...', va='top', ha='left',
                        fontsize=tiny_fontsize,transform=ax.transAxes)
                break

    print('making:', pngfile)
    cm.savefig(pngfile, dpi=dpi)

    num_neighbors = nbrs.shape[1] if compute_nbr_averages else 0

    help_message = f"""This plot shows the distribution of significant
    features from graph-vs-features or HotSpot analysis plotted across the
    {dist_tag.upper()} landscape. Rows are features and columns are
    individual clonotypes. Columns are ordered by hierarchical clustering
    (if a dendrogram is present above the heatmap) or by a 1D UMAP projection
    (used for very large datasets or if 'X_pca_{dist_tag}' is not present in
    adata.obsm_keys()). Rows are ordered by hierarchical clustering with
    a correlation metric.

    The row colors to the left of the heatmap show the feature type
    (blue=TCR, orange=GEX). The row colors to the left of those
    indicate the strength of the graph-vs-feature correlation
    (also included in the feature labels to the right of the heatmap;
    keep in mind that highly significant P values for some features may shift
    the colorscale so everything else looks dark blue).

    The column colors above the heatmap are {dist_tag.upper()} clusters
    (and TCR V/J genes if plotting against the TCR landscape). The text
    above the column colors provides more info.

    Feature scores are Z-score normalized and then averaged over the
    K={num_neighbors} nearest neighbors (0 means no nbr-averaging).

    The 'coolwarm' colormap is centered at Z=0.

    Since features of the same type (GEX or TCR) as the landscape and
    neighbor graph (ie {dist_tag.upper()} features) are more highly
    correlated over graph neighborhoods, their neighbor-averaged scores
    will show more extreme variation. For this reason, the nbr-averaged
    scores for these features from the same modality as the landscape
    itself are downscaled by a factor of
    rescale_factor_for_self_features={rescale_factor_for_self_features:.2f}.

    The colormap in the top left is for the Z-score normalized,
    neighbor-averaged scores (multiply by {1.0/rescale_factor_for_self_features:.2f}
    to get the color scores for the {dist_tag.upper()} features).

    """

    return help_message

def plot_cluster_gene_compositions(
        adata,
        pngfile
):
    ''' Show for each GEX cluster the top VJ genes, and the TCR cluster composition
    Show GEX cluster composition for each TCR cluster
    '''
    organism = adata.uns['organism']
    num_clones = adata.shape[0]
    organism_genes = all_genes[organism]

    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    gene_colors, sorted_gene_colors = assign_colors_to_conga_tcrs(tcrs, organism, return_sorted_color_tuples=True)

    # gex cluster compositions
    num_clusters_gex = np.max(clusters_gex)+1
    num_clusters_tcr = np.max(clusters_tcr)+1
    #ncols = 1+num_clusters_gex
    #nrows = 5 # tcr cluster, va, ja, vb, jb

    plt.figure(figsize=(10,6))

    for col in range(5):
        plt.subplot(1,5,col+1)
        #
        if col==0:
            features = clusters_tcr
            color_dict = get_integers_color_dict(num_clusters_tcr) # map from feature to color
            feature_name = 'tcr cluster'

        else:
            ii_ab = (col-1)//2
            ii_vj = (col-1)%2
            features = np.array([ organism_genes[x[ii_ab][ii_vj]].count_rep for x in tcrs])
            color_dict = dict(sorted_gene_colors[col-1])
            feature_name = 'VA JA VB JB'.split()[col-1]


        centers, bottoms, heights, colors = [], [], [], []
        for ii in range(-1, num_clusters_gex):
            if ii==-1:
                mask = np.full((num_clones,), True)
            else:
                mask = clusters_gex==ii

            counts = Counter( features[mask])
            total = np.sum(mask)

            if ii==-1:
                feature_order = [x[0] for x in counts.most_common()]

            bottom=1.0
            for feature in feature_order:
                count = counts[feature]
                height = float(count)/total
                bottom -= height
                centers.append(ii)
                bottoms.append(bottom)
                heights.append(height)
                colors.append(color_dict[feature])

        plt.bar(centers, height=heights, width=0.8, bottom=bottoms, color=colors)
        plt.title(feature_name)
        if col:
            plt.yticks([],[])
        plt.xticks(range(-1, num_clusters_gex),
                   ['A']+[str(x) for x in range(num_clusters_gex)])
        plt.xlabel('GEX cluster (or _All)', fontsize=6, labelpad=0)

        # make a legend
        MAX_TO_SHOW = 25
        for feature in feature_order[:MAX_TO_SHOW]:
            plt.scatter([], [], c=[color_dict[feature]], label=feature)

        plt.legend(loc = 'upper center', bbox_to_anchor=(0.5, -0.1),
                   fontsize=6, ncol=2)

    plt.subplots_adjust(bottom=0.35, right=0.96, left= 0.05, wspace=0.15)
    print('making:', pngfile)
    plt.savefig(pngfile)


def make_cluster_logo_plots_figure(
        adata,
        scores,
        max_good_score,
        clusters_gex,
        clusters_tcr,
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size,
        pngfile,
        figure_tag = None,
        **kwargs # passed to make_logo_plots
):
    util.setup_uns_dicts(adata) # convenience

    cluster_tcr_score_names = tcr_scoring.all_tcr_scorenames[:]

    good_mask = (scores <= max_good_score)

    # bic is short for 'bicluster'
    bic_counts = Counter(
        (x,y) for x,y,m in zip(clusters_gex, clusters_tcr, good_mask) if m)

    num_good_biclusters = sum(1 for x,y in bic_counts.items()
                              if y>=min_cluster_size )

    if num_good_biclusters:
        # calc tcr sequence features of good cluster pairs
        good_bicluster_tcr_scores = correlations.calc_good_cluster_tcr_features(
            adata, good_mask, clusters_gex, clusters_tcr,
            cluster_tcr_score_names, min_count=min_cluster_size)

        # run rank_genes on most common bics
        if 'rank_genes_uns_tag' in kwargs:
            rank_genes_uns_tag = kwargs['rank_genes_uns_tag']
            del kwargs['rank_genes_uns_tag']# since we pass it explicitly below
        else:
            rank_genes_uns_tag = (
                f'rg_{num_good_biclusters}_{np.sum(good_mask)}_biclusters')

        correlations.run_rank_genes_on_good_biclusters(
            adata, good_mask, clusters_gex, clusters_tcr,
            min_count=min_cluster_size, key_added= rank_genes_uns_tag)

        help_message = make_logo_plots(
            adata, nbrs_gex, nbrs_tcr, min_cluster_size, pngfile,
            conga_scores = scores,
            clusters_gex = clusters_gex,
            clusters_tcr = clusters_tcr,
            good_score_mask = good_mask,
            good_bicluster_tcr_scores = good_bicluster_tcr_scores,
            rank_genes_uns_tag = rank_genes_uns_tag,
            **kwargs)

        if figure_tag is not None:
            adata.uns['conga_results'][figure_tag] = pngfile
            adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
            util.make_figure_helpfile(figure_tag, adata)





def make_tcr_clumping_plots(
        adata,
        nbrs_gex,
        nbrs_tcr,
        outfile_prefix,
        min_cluster_size_for_logos=3,
        pvalue_threshold_for_logos=1.0, #pvalues are crude bonferroni corrected
        max_color_pvalue = 1e-16, # in the UMAP figure; no brighter beyond there
        **logo_plot_args,
):
    num_clones = adata.shape[0]
    fake_clusters_gex = np.zeros((num_clones,)).astype(int)
    fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
    clumping_pvals = np.full( (num_clones,), num_clones).astype(float)

    if ('conga_results' not in adata.uns_keys() or
        TCR_CLUMPING not in adata.uns['conga_results']):
        print('make_tcr_clumping_plots:: no results in adata',
              'did you call assess_tcr_clumping first?')
        return

    results_df = adata.uns['conga_results'][TCR_CLUMPING]

    for l in results_df.itertuples():
        clumping_pvals[ l.clone_index] = min(l.pvalue_adj,
                                             clumping_pvals[l.clone_index])
        fake_clusters_tcr[l.clone_index] = l.clumping_group

    figure_tag = TCR_CLUMPING_LOGOS
    pngfile = f'{outfile_prefix}_{figure_tag}.png'

    if 'rank_genes_uns_tag' not in logo_plot_args:
        logo_plot_args['rank_genes_uns_tag'] = 'rg_tcr_clumping_biclusters'

    if 'conga_scores_name' not in logo_plot_args:
        logo_plot_args['conga_scores_name'] = 'TCR clumping'

    if 'show_real_clusters_gex' not in logo_plot_args:
        logo_plot_args['show_real_clusters_gex'] = True

    make_cluster_logo_plots_figure(
        adata, clumping_pvals, pvalue_threshold_for_logos,
        fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile, figure_tag=figure_tag,
        **logo_plot_args)

    # make umaps colored by clumping pvals
    plt.figure(figsize=(12,6))
    for icol, tag in enumerate(['gex','tcr']):
        plt.subplot(1,2,icol+1)
        colors = np.sqrt( np.maximum(
            0.0, -1*np.log10(np.maximum(1e-100, clumping_pvals))))

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
        plt.title('TCR clumping pvalues')
    plt.tight_layout()
    pngfile = '{}_tcr_clumping.png'.format(outfile_prefix)
    print('making:', pngfile)
    plt.savefig(pngfile)




def make_tcrdist_trees(
        adata,
        output_prefix,
        group_by = None,
        max_tcrs_per_cluster_for_tcrdists=5000, # if more, will subsample
):
    """
    generate TCRdist trees by gene expression cluster
    group_by: use 'clusters_gex' or 'clusters_tcr' to generate trees
       by respective cluster assignments.
    """
    if group_by is None or group_by == 'clusters_gex':
        group_by = 'clusters_gex'
        tag = 'GEX'
        figure_tag = GEX_CLUSTERS_TCRDIST_TREES
    elif group_by == 'clusters_tcr':
        group_by = 'clusters_tcr'
        tag = 'TCR'
        figure_tag = TCR_CLUSTERS_TCRDIST_TREES

    width = 800
    height = 1000
    xpad = 25
    organism = adata.uns['organism']
    clusters = np.array(adata.obs[group_by])

    num_clusters = np.max(clusters)+1
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    num_clones = adata.shape[0]
    if 'conga_scores' in adata.obs_keys():
        conga_scores = np.maximum(1e-100, np.array(adata.obs['conga_scores']))
        scores = np.sqrt(
            np.maximum(0.0, -1*np.log10(100*conga_scores/num_clones)))
    else:
        scores = np.zeros((adata.shape[0],))

    tcrdister = TcrDistCalculator(organism)

    x_offset = 0
    all_cmds = []

    #color_score_range = [-1*np.log(10), -1*np.log(1e-5)]
    color_score_range = [0, 3.0]

    for clust in range(num_clusters):
        cmask = (clusters==clust)
        csize = np.sum(cmask)

        if csize>max_tcrs_per_cluster_for_tcrdists:
            print('make_tcrdist_trees: too many TCRs in cluster, subsampling',
                  csize, '>', max_tcrs_per_cluster_for_tcrdists)
            cinds = list(np.nonzero(cmask)[0])
            random.shuffle(cinds)
            cmask = np.full((num_clones,), False)
            cmask[cinds[:max_tcrs_per_cluster_for_tcrdists]] = True
            csize = np.sum(cmask)
            assert csize == max_tcrs_per_cluster_for_tcrdists

        ctcrs   = [x for x,y in zip(  tcrs, cmask) if y]
        cscores = [x for x,y in zip(scores, cmask) if y]

        print('computing tcrdist distances:', clust, csize)
        if csize>1000 and util.tcrdist_cpp_available():
            cdists = preprocess.calc_tcrdist_matrix_cpp(
                ctcrs, adata.uns['organism'])
        else:
            cdists = np.array([tcrdister(x,y) for x in ctcrs for y in ctcrs])\
                       .reshape(csize,csize)

        cmds = make_tcr_tree_svg_commands(
            ctcrs, organism, [x_offset,0], [width,height], cdists,
            max_tcrs_for_trees=400, tcrdist_calculator=tcrdister,
            color_scores=cscores, color_score_range = color_score_range,
            title=f'{tag} cluster {clust}')

        x_offset += width + xpad

        all_cmds.extend(cmds)

    pngfile = f'{output_prefix}_{figure_tag}.png'
    print('making:', pngfile)
    svg_basic.create_file(all_cmds, x_offset-xpad, height, pngfile[:-3]+'svg',
                          create_png= True )

    def conga_score_from_color_score(x):
        return 10**(-1*(x**2)) * num_clones / 100.
    loscore = conga_score_from_color_score(color_score_range[0])
    hiscore = conga_score_from_color_score(color_score_range[1])

    help_message = f"""These are TCRdist hierarchical clustering trees
for the {tag} clusters (cluster assignments stored in
adata.obs['{group_by}']). The trees are colored by CoNGA score
with a color score range of {loscore:.2e} (blue) to {hiscore:.2e} (red).
For coloring, CoNGA scores are log-transformed, negated, and square-rooted
(with an offset in there, too, roughly speaking).
"""

    adata.uns.setdefault('conga_results',{})[figure_tag] = pngfile
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)



def make_tcrdist_tree_for_conga_score_threshold(
        adata,
        threshold,
        outfile_prefix,
        min_clonotypes=10,
        max_tcrs_for_tcrdists=5000, # if more, will subsample
):
    ''' Make a tcrdist tree for clonotypes with conga_scores below a threshold
    '''

    if 'conga_scores' not in adata.obs_keys():
        print('ERROR make_tcrdist_tree_for_conga_hits:',
              'conga_scores is not in adata.obs')
        return

    organism = adata.uns['organism']
    tcrdister = TcrDistCalculator(organism)

    conga_scores = np.array(adata.obs['conga_scores'])
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    # recalibrate the scores to give useful colors.
    scores = np.sqrt(np.maximum(0.0, -1*np.log10(conga_scores/threshold)))
    color_score_range = [0, 3.0] #max(3.0, np.max(scores))]
    cmask = (conga_scores<=threshold)
    csize = np.sum(cmask)
    if csize < min_clonotypes:
        print('make_tcrdist_tree_for_conga_hits: too few hits:', csize,
              'min_clonotypes=', min_clonotypes)
        return
    elif csize>max_tcrs_for_tcrdists:
        print('make_tcrdist_tree_for_conga_hits: too many TCRs for tcrdists,'
              'subsampling:', csize, '>', max_tcrs_for_tcrdists)
        cinds = list(np.nonzero(cmask)[0])
        random.shuffle(cinds)
        cmask = np.full((adata.shape[0],), False)
        cmask[cinds[:max_tcrs_for_tcrdists]] = True
        csize = np.sum(cmask)
        assert csize == max_tcrs_for_tcrdists

    ctcrs   = [x for x,y in zip(  tcrs, cmask) if y]
    cscores = [x for x,y in zip(scores, cmask) if y]

    if csize>1000 and util.tcrdist_cpp_available():
        cdists = preprocess.calc_tcrdist_matrix_cpp(
            ctcrs, adata.uns['organism'])
    else:
        print('computing tcrdist distances:', csize)
        cdists = np.array([tcrdister(x,y) for x in ctcrs for y in ctcrs])\
                   .reshape(csize,csize)

    width = 800
    height = 1000
    cmds = make_tcr_tree_svg_commands(
        ctcrs, organism, [0,0], [width,height], cdists,
        max_tcrs_for_trees=400, tcrdist_calculator=tcrdister,
        color_scores=cscores, color_score_range = color_score_range,
        title=f'conga_score_threshold {threshold:.1f}')

    figure_tag = CONGA_THRESHOLD_TCRDIST_TREE
    pngfile = f'{outfile_prefix}_{figure_tag}.png'
    print('making:', pngfile)
    svg_basic.create_file(cmds, width, height, pngfile[:-3]+'svg',
                          create_png=True )

    def conga_score_from_color_score(x):
        return 10**(-1*(x**2)) * threshold
    loscore = conga_score_from_color_score(color_score_range[0])
    hiscore = conga_score_from_color_score(color_score_range[1])

    help_message = f"""This is a TCRdist hierarchical clustering tree
for the clonotypes with CoNGA score less than {threshold}.
The tree is colored by CoNGA score
with a color score range of {loscore:.2e} (blue) to {hiscore:.2e} (red).
For coloring, CoNGA scores are log-transformed, negated, and square-rooted
(with an offset in there, too, roughly speaking).
"""
    adata.uns.setdefault('conga_results',{})[figure_tag] = pngfile
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)



def make_batch_colored_umaps(
        adata,
        outfile_prefix,
):
    ''' Note that this just plots the batch assignment of the center
    (representative) clone for each clonotype.
    '''
    if 'batch_keys' not in adata.uns:
        print("WARNING: make_batch_colored_umaps called,",
              "but no 'batch_keys' entry in adata.uns")
        return

    batch_keys = adata.uns['batch_keys']

    figure_tag = BATCH_UMAPS
    pngfile = f'{outfile_prefix}_{figure_tag}.png'

    nrows, ncols, plotno = 2, len(batch_keys), 0
    plt.figure(figsize=(ncols*4, nrows*4))
    for xy_tag in ['gex','tcr']:
        xy = adata.obsm[f'X_{xy_tag}_2d']
        for batch_key in batch_keys:
            batches = adata.obs[batch_key]
            num_batch_vals = np.max(batches)+1
            if num_batch_vals <= 10:
                cmap_colors = plt.get_cmap('tab10').colors
            else:
                cmap_colors = plt.get_cmap('tab20').colors

            # plot clonotypes in random order
            reorder = np.random.permutation(adata.shape[0])

            plotno += 1
            plt.subplot(nrows, ncols, plotno)
            colors = np.array([cmap_colors[x%len(cmap_colors)]
                               for x in batches])
            plt.scatter(xy[reorder,0], xy[reorder,1], c=colors[reorder],
                        s=5)
            unique_batches = sorted(set(batches))
            add_categorical_legend(
                plt.gca(), [str(x) for x in unique_batches],
                [cmap_colors[x%len(cmap_colors)] for x in unique_batches],
                #legend_fontsize=6,
            )

            plt.title(f'{batch_key}')
            plt.xlabel(f'{xy_tag.upper()} UMAP1')
            plt.ylabel(f'{xy_tag.upper()} UMAP2')
            plt.xticks([],[])
            plt.yticks([],[])
    #plt.tight_layout()
    print('making:', pngfile)
    plt.savefig(pngfile)
    help_message = """This figure shows the TCR and GEX 2D UMAP landscapes
    colored by the batch assignment of the representative cell for each
    clonotype. The batches are the ones present in adata.uns['batch_keys']
    and are given in the titles of the individual panels.
    """
    if 'conga_results' not in adata.uns:
        util.setup_uns_dicts(data)
    adata.uns['conga_results'][figure_tag] = pngfile
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)


## this function is really similar to make_cluster_logo_plots_figure
##  should condense at some point to avoid duplication...
##
def make_graph_vs_graph_logos(
        adata,
        outfile_prefix,
        min_cluster_size,
        nbrs_gex,
        nbrs_tcr,
        gex_nbrhood_tcr_score_names = tcr_scoring.all_tcr_scorenames,
        gex_header_tcr_score_names = None,
        **make_logo_kwargs,
):
    util.setup_uns_dicts(adata) # makes life easier

    # the conga scores
    conga_scores = np.array(adata.obs['conga_scores'])
    good_mask = (conga_scores <= 1.0)

    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    bic_counts = Counter((x,y) for x,y,m in zip(clusters_gex, clusters_tcr,
                                                good_mask) if m )

    num_good_biclusters = sum(1 for x,y in bic_counts.items()
                              if y>=min_cluster_size )

    adata.uns['conga_stats']['num_gvg_hit_clonotypes'] = np.sum(good_mask)
    adata.uns['conga_stats']['num_gvg_hit_biclusters'] = num_good_biclusters


    if not num_good_biclusters:
        return ################################# EARLY RETURN !!!!!!!!!!!!!

    # calc tcr sequence features of good cluster pairs
    good_bicluster_tcr_scores = correlations.calc_good_cluster_tcr_features(
        adata, good_mask, clusters_gex, clusters_tcr,
        gex_nbrhood_tcr_score_names, min_count=min_cluster_size)

    # calc differentially expressed genes for good cluster pairs
    rank_genes_uns_tag = 'rank_genes_good_biclusters'
    correlations.run_rank_genes_on_good_biclusters(
        adata, good_mask, clusters_gex, clusters_tcr,
        min_count=min_cluster_size, key_added= rank_genes_uns_tag)

    if gex_header_tcr_score_names is None:
        if '_ig' in adata.uns['organism']:
            gex_header_tcr_score_names = ['af2', 'cdr3len', 'volume',
                                          'nndists_tcr']
        else:
            gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8',
                                          'nndists_tcr']

    lit_matches = adata.uns['conga_results'].get(TCR_DB_MATCH, None)

    figure_tag = GRAPH_VS_GRAPH_LOGOS
    pngfile = f'{outfile_prefix}_{figure_tag}.png'

    help_message = make_logo_plots(
        adata,
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size,
        pngfile,
        rank_genes_uns_tag = rank_genes_uns_tag,
        good_bicluster_tcr_scores = good_bicluster_tcr_scores,
        gex_header_tcr_score_names = gex_header_tcr_score_names,
        lit_matches=lit_matches,
        **make_logo_kwargs,
    )

    # note that there is an early return if no good clusters
    adata.uns['conga_results'][figure_tag] = pngfile
    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
    util.make_figure_helpfile(figure_tag, adata)


def make_graph_vs_features_plots(
        adata,
        all_nbrs,
        outfile_prefix,
        max_clones_for_dendrograms = 10000, # for memory
        verbose=False,
        clustermap_pvalue_threshold=0.05, # adjusted pvalues
        clustermap_pvalue_threshold_big=1e-3, # for really big feature sets
        clustermap_pvalue_threshold_big_size=500, # ie, more than 500 fts
        clustermap_max_type_features=50,
        min_nbrs_for_averaging=30,
):
    util.setup_uns_dicts(adata) # should already be done


    ## tcr graph vs gex features #####################3
    tcr_graph_results = adata.uns['conga_results'].get(
        TCR_GRAPH_VS_GEX_FEATURES, pd.DataFrame([]))
    tcr_genes_results = adata.uns['conga_results'].get(
        TCR_GENES_VS_GEX_FEATURES, pd.DataFrame([]))
    gex_graph_results = adata.uns['conga_results'].get(
        GEX_GRAPH_VS_TCR_FEATURES, pd.DataFrame([]))

    tcr_graph_results['feature_type'] = 'gex'
    tcr_genes_results['feature_type'] = 'gex'
    gex_graph_results['feature_type'] = 'tcr'

    combo_results = pd.concat(
        [tcr_graph_results, tcr_genes_results, gex_graph_results])

    if combo_results.shape[0]>1: # make combo clustermap
        # choose a smallish nbr_frac for nbr-averaging in the plots
        nbr_fracs = sorted(all_nbrs.keys())
        for nbr_frac in nbr_fracs:
            if nbr_frac*adata.shape[0] >= min_nbrs_for_averaging:
                plot_nbrs_gex, plot_nbrs_tcr = all_nbrs[nbr_frac]
                break
        else:
            plot_nbrs_gex, plot_nbrs_tcr = all_nbrs[max(nbr_fracs)]

        # sort, drop duplicates (due to multiple nbr_fracs)
        combo_results = combo_results.sort_values('mwu_pvalue_adj')\
                                     .drop_duplicates(['feature',
                                                       'feature_type'])
        # threshold on mwu_pvalue_adj
        combo_results = combo_results[combo_results.mwu_pvalue_adj <=
                                      clustermap_pvalue_threshold]

        if combo_results.shape[0] > clustermap_pvalue_threshold_big_size:
            print('using stringent pval threshold for big clustermap:',
                  combo_results.shape[0], clustermap_pvalue_threshold_big)
            combo_results = combo_results[combo_results.mwu_pvalue_adj <=
                                          clustermap_pvalue_threshold_big]

        for dist_tag in ['gex','tcr']:
            # dist_tag is the similarity measure being used for column ordering

            use_dendrogram_for_clustermap_cell_order = (
                f'X_pca_{dist_tag}' in adata.obsm_keys() and
                adata.shape[0] <= max_clones_for_dendrograms)

            use_1d_landscape_for_clustermap_cell_order = (
                f'X_{dist_tag}_1d' in adata.obsm_keys() and
                not use_dendrogram_for_clustermap_cell_order)

            if not (use_dendrogram_for_clustermap_cell_order or
                    use_1d_landscape_for_clustermap_cell_order):
                # cant make clustermap
                continue


            features = list(combo_results.feature)
            feature_types = list(combo_results.feature_type)

            min_pval = 1e-299 # dont want log10 of 0.0
            feature_pvalues = np.maximum(
                min_pval, np.array(combo_results.mwu_pvalue_adj))
            feature_scores = np.sqrt(-1*np.log10(feature_pvalues))

            feature_labels = [
                '{:9.1e} {} {}'.format(x,y,z)
                for x,y,z in zip(combo_results.mwu_pvalue_adj,
                                 combo_results.feature_type,
                                 combo_results.feature)]

            figure_tag = (GRAPH_VS_FEATURES_GEX_CLUSTERMAP
                          if dist_tag=='gex' else
                          GRAPH_VS_FEATURES_TCR_CLUSTERMAP)

            pngfile = f'{outfile_prefix}_{figure_tag}.png'
            nbrs = plot_nbrs_gex if dist_tag == 'gex' else plot_nbrs_tcr

            help_message = plot_interesting_features_vs_clustermap(
                adata, features, feature_types, pngfile, dist_tag,
                nbrs=nbrs,
                compute_nbr_averages=True,
                feature_labels=feature_labels,
                feature_scores=feature_scores,
                max_type_features=clustermap_max_type_features,
                use_1d_landscape_for_cell_order=
                    use_1d_landscape_for_clustermap_cell_order,
                title=figure_tag,
                verbose=verbose,
            )
            adata.uns['conga_results'][figure_tag] = pngfile
            adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
            util.make_figure_helpfile(figure_tag, adata)

    if not tcr_graph_results.empty:
        # labeled umap
        figure_tag = TCR_GRAPH_VS_GEX_FEATURES_PLOT
        pngfile = f'{outfile_prefix}_{figure_tag}.png'
        help_message = plot_ranked_strings_on_cells(
            adata, tcr_graph_results, 'X_tcr_2d', 'clone_index',
            'mwu_pvalue_adj', 1.0, 'feature', pngfile,
            title=figure_tag)
        adata.uns['conga_results'][figure_tag] = pngfile
        adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
        util.make_figure_helpfile(figure_tag, adata)

        # panels
        figure_tag = TCR_GRAPH_VS_GEX_FEATURES_PANELS
        pngfile = f'{outfile_prefix}_{figure_tag}.png'
        xy_tag, feature_type = 'tcr','gex'
        help_message = make_feature_panel_plots(
            adata, xy_tag, all_nbrs, tcr_graph_results, pngfile, feature_type,
            title=figure_tag)
        adata.uns['conga_results'][figure_tag] = pngfile
        adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
        util.make_figure_helpfile(figure_tag, adata)

    else:
        print('conga.plotting.make_graph_vs_features_plots:: missing results',
              'dont forget to call conga.correlations.run_graph_vs_features')

    ## tcr genes vs gex features ##############3
    if not tcr_genes_results.empty:
        # panels
        figure_tag = TCR_GENES_VS_GEX_FEATURES_PANELS
        pngfile = f'{outfile_prefix}_{figure_tag}.png'
        xy_tag, feature_type = 'tcr','gex'
        help_message = make_feature_panel_plots(
            adata, xy_tag, all_nbrs, tcr_genes_results, pngfile, feature_type,
            use_nbr_frac = max(all_nbrs.keys()),
            title=figure_tag)
        adata.uns['conga_results'][figure_tag] = pngfile
        adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
        util.make_figure_helpfile(figure_tag, adata)
    else:
        print('conga.plotting.make_graph_vs_features_plots:: missing results',
              'dont forget to call conga.correlations.run_graph_vs_features')


    ## gex graph vs tcr features #####################
    if not gex_graph_results.empty:
        # labeled umap
        figure_tag = GEX_GRAPH_VS_TCR_FEATURES_PLOT
        pngfile = f'{outfile_prefix}_{figure_tag}.png'
        help_message = plot_ranked_strings_on_cells(
            adata, gex_graph_results, 'X_gex_2d', 'clone_index',
            'mwu_pvalue_adj', 1.0, 'feature', pngfile,
            direction_column='ttest_stat',
            title=figure_tag)
        adata.uns['conga_results'][figure_tag] = pngfile
        adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
        util.make_figure_helpfile(figure_tag, adata)

        # panels
        figure_tag = GEX_GRAPH_VS_TCR_FEATURES_PANELS
        pngfile = f'{outfile_prefix}_{figure_tag}.png'
        xy_tag, feature_type = 'gex', 'tcr'
        help_message = make_feature_panel_plots(
            adata, xy_tag, all_nbrs, gex_graph_results, pngfile, feature_type,
            title=figure_tag)
        adata.uns['conga_results'][figure_tag] = pngfile
        adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = help_message
        util.make_figure_helpfile(figure_tag, adata)
    else:
        print('conga.plotting.make_graph_vs_features_plots:: missing results',
              'dont forget to call conga.correlations.run_graph_vs_features')


def make_hotspot_plots(
        adata,
        all_nbrs,
        outfile_prefix,
        make_raw_feature_plots=False,
        max_clones_for_dendrograms=20000, # control memory
        min_nbrs_for_averaging=30,
        clustermap_max_type_features=50,
):
    ''' Call conga.correlations.find_hotspots_wrapper first!
    '''
    util.setup_uns_dicts(adata) # probably unnecessary

    if HOTSPOT_FEATURES not in adata.uns['conga_results']:
        print('plotting.make_hotspot_plots:: no hotspot results,',
              'did you call conga.correlations.find_hotspots_wrapper first?')
        return

    all_results = adata.uns['conga_results'][HOTSPOT_FEATURES]

    nbr_fracs = sorted(all_nbrs.keys())

    # choose a smallish nbr_frac for nbr-averaging in the plots
    for nbr_frac in nbr_fracs:
        if nbr_frac*adata.shape[0] >= min_nbrs_for_averaging:
            plot_nbrs_gex, plot_nbrs_tcr = all_nbrs[nbr_frac]
            break
    else:
        plot_nbrs_gex, plot_nbrs_tcr = all_nbrs[max(nbr_fracs)]

    for nbr_frac in nbr_fracs:
        # unpack the results for this nbr_frac:
        combo_results = all_results[all_results.nbr_frac == nbr_frac]
        gex_results = combo_results[combo_results.feature_type=='gex']
        tcr_results = combo_results[combo_results.feature_type=='tcr']

        #nbrs_gex, nbrs_tcr = all_nbrs[nbr_frac]

        for tag, results in [ ['gex', gex_results],
                              ['tcr', tcr_results],
                              ['combo', combo_results] ]:
            if results.shape[0]<1:
                continue

            for plot_tag, plot_nbrs in [['gex',plot_nbrs_gex],
                                        ['tcr',plot_nbrs_tcr]]:
                if tag == plot_tag:
                    continue
                pngfile_prefix = '{}_hotspot_{}_features_{:.3f}_nbrs_{}_plot'\
                                 .format(outfile_prefix, tag, nbr_frac,
                                         plot_tag)
                # 2D UMAPs colored by nbr-averaged feature values
                pngfile = f'{pngfile_prefix}_umap_nbr_avg.png'
                print('making:', pngfile)
                help_message = plot_hotspot_umap(
                    adata, plot_tag, results, pngfile, nbrs=plot_nbrs,
                    compute_nbr_averages=True,
                    title=f'HotSpot {tag} features vs {plot_tag} landscape '
                    '(NBR-avged)')

                if nbr_frac == max(nbr_fracs) and tag == 'combo':
                    figure_tag = (HOTSPOT_GEX_UMAP if plot_tag=='gex' else
                                  HOTSPOT_TCR_UMAP)
                    adata.uns['conga_results'][figure_tag] = pngfile
                    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = (
                        help_message)
                    util.make_figure_helpfile(figure_tag, adata)


                if make_raw_feature_plots:
                    pngfile = f'{pngfile_prefix}_umap_raw.png'
                    print('making:', pngfile)
                    plot_hotspot_umap(
                        adata, plot_tag, results, pngfile, nbrs=plot_nbrs,
                        compute_nbr_averages=False,
                        title=f'HotSpot {tag} features vs {plot_tag} landscape '
                        '(raw)')

                if results.shape[0]<2:
                    continue # clustermap not interesting...

                use_dendrogram_for_clustermap_cell_order = (
                    f'X_pca_{plot_tag}' in adata.obsm_keys() and
                    adata.shape[0] <= max_clones_for_dendrograms)

                use_1d_landscape_for_clustermap_cell_order = (
                    f'X_{plot_tag}_1d' in adata.obsm_keys() and
                    not use_dendrogram_for_clustermap_cell_order)

                if not (use_1d_landscape_for_clustermap_cell_order or
                        use_dendrogram_for_clustermap_cell_order):
                    # can't make a clustermap: too big, or no data for dists
                    continue

                ## clustermap of features versus cells
                features = list(results.feature)
                feature_types = list(results.feature_type)
                feature_labels = ['{:9.1e} {} {}'.format(x,y,z)
                                  for x,y,z in zip(results.pvalue_adj,
                                                   results.feature_type,
                                                   results.feature)]
                min_pval = 1e-299 # dont want log10 of 0.0
                feature_pvalues = np.maximum(
                    min_pval, np.array(results.pvalue_adj))
                feature_scores = np.sqrt(-1*np.log10(feature_pvalues))

                # nbr-averaged version
                pngfile = f'{pngfile_prefix}_clustermap_nbr_avg.png'

                help_message = plot_interesting_features_vs_clustermap(
                    adata, features, feature_types, pngfile, plot_tag,
                    nbrs=plot_nbrs,
                    compute_nbr_averages=True,
                    feature_labels=feature_labels,
                    max_type_features=clustermap_max_type_features,
                    feature_scores=feature_scores,
                    use_1d_landscape_for_cell_order=
                        use_1d_landscape_for_clustermap_cell_order,
                    title=f'HotSpot {tag} features vs {plot_tag} landscape '
                    '(NBR-avged)')

                if nbr_frac == max(nbr_fracs) and tag == 'combo':
                    figure_tag = (HOTSPOT_GEX_CLUSTERMAP if plot_tag=='gex' else
                                  HOTSPOT_TCR_CLUSTERMAP)
                    adata.uns['conga_results'][figure_tag] = pngfile
                    adata.uns['conga_results'][figure_tag+HELP_SUFFIX] = (
                        help_message)
                    util.make_figure_helpfile(figure_tag, adata)

                if make_raw_feature_plots:
                    # RAW scores (non-nbr-averaged) version
                    pngfile = f'{pngfile_prefix}_clustermap_raw.png'
                    help_message = plot_interesting_features_vs_clustermap(
                        adata, features, feature_types,
                        pngfile, plot_tag,
                        compute_nbr_averages=False,
                        feature_labels=feature_labels,
                        max_type_features=clustermap_max_type_features,
                        feature_scores=feature_scores,
                        use_1d_landscape_for_cell_order=
                            use_1d_landscape_for_clustermap_cell_order,
                        title=f'HotSpot {tag} features vs {plot_tag} landscape '
                        '(raw)')




def _add_html_newlines(msg):
    if msg is None:
        return msg
    else:
        return msg.replace('\n\n','\n<br><br>\n')


# the order of the elements in the html summary
default_content_order = [
    GRAPH_VS_GRAPH_STATS,
    GRAPH_VS_GRAPH,
    GRAPH_VS_GRAPH_LOGOS,
    TCR_CLUMPING,
    TCR_CLUMPING_LOGOS,
    TCR_DB_MATCH,
    TCR_DB_MATCH_PLOT,
    TCR_GRAPH_VS_GEX_FEATURES,
    TCR_GRAPH_VS_GEX_FEATURES_PLOT,
    TCR_GRAPH_VS_GEX_FEATURES_PANELS,
    TCR_GENES_VS_GEX_FEATURES,
    TCR_GENES_VS_GEX_FEATURES_PANELS,
    GEX_GRAPH_VS_TCR_FEATURES,
    GEX_GRAPH_VS_TCR_FEATURES_PLOT,
    GEX_GRAPH_VS_TCR_FEATURES_PANELS,
    GRAPH_VS_FEATURES_GEX_CLUSTERMAP,
    GRAPH_VS_FEATURES_TCR_CLUSTERMAP,
    GRAPH_VS_SUMMARY,
    GEX_CLUSTERS_TCRDIST_TREES,
    CONGA_THRESHOLD_TCRDIST_TREE,
    HOTSPOT_FEATURES,
    HOTSPOT_GEX_UMAP,
    HOTSPOT_GEX_CLUSTERMAP,
    HOTSPOT_TCR_UMAP,
    HOTSPOT_TCR_CLUSTERMAP,
    BATCH_UMAPS,
]

# figure tags:

def make_html_summary(
        adata,
        html_filename,
        command_string = None,
        title = None,
        content_order = default_content_order,
        max_dataframe_lines = 50,
        show_table_indices = False,
):
    '''
    '''
    if 'conga_results' not in adata.uns.keys():
        print('ERROR conga.plotting.make_html_summary:',
              'conga_results is missing from adata.uns')
        return


    if command_string is None:
        command_section = ''
    else:
        command_section = f"""
<h1>Command:</h1>
{command_string}
<br>
"""
    if title is None:
        title = 'CoNGA results'

    table_style = '' # configure later

    print('making:', html_filename)
    out = open(html_filename, 'w')
    out.write(f"""<!doctype html>
    <html>
    <head>
    <title>{title}</title>
    {table_style}
    </head>
    <body>

    <h1>Overview:</h1>

    This page contains the results of CoNGA analyses.
    Results in tables may have been filtered to reduce redundancy,
    focus on the most important columns, and
    limit length; full tables should exist as OUTFILE_PREFIX*.tsv files.
    <br>
    <br>

    {command_section}
    """)

    if 'conga_stats' in adata.uns_keys():
        #write out the stats
        out.write('<h1>Stats</h1>\n')
        for tag,val in adata.uns['conga_stats'].items():
            out.write(f'{tag}: {val} <br>\n')

    for content_tag in content_order:
        if content_tag not in adata.uns['conga_results']:
            print('make_html_summary:: no results for content_tag', content_tag)
            continue

        results = adata.uns['conga_results'][content_tag]

        if type(results) is str:
            # filename of an image
            # assume we don't need the path/folder stuff for this html file
            pngfile = Path(results).name

            # hacking: limit the width/height on specific plots
            extra_tag = ''
            if content_tag.endswith('_logos'): # default is too big
                extra_tag += ' width="1250" '
            if content_tag.endswith('_clustermap'): # default is too big
                extra_tag += ' width="1500" '
            if '_tcrdist_tree' in content_tag: # full size is 1000, too big?
                extra_tag += ' height="750" '
            out.write(f'<h1>{content_tag}</h1>\n<br>\n')

            # write help
            msg = adata.uns['conga_results'].get(content_tag+HELP_SUFFIX, '')
            msg = _add_html_newlines(msg) # reformat for html
            if msg:
                out.write(f'{msg}<br>\n')

            # now the image
            out.write(f'Image source: {pngfile}<br>\n')
            out.write(f'<img src="{pngfile}" {extra_tag} />\n')

        else: # results is a pandas DataFrame
            out.write(f'<h1>{content_tag}</h1>\n<br>\n')
            msg = adata.uns['conga_results'].get(content_tag+HELP_SUFFIX, '')
            msg = _add_html_newlines(msg)
            if msg:
                out.write(f'{msg}<br>\n')
            out.write(results.head(max_dataframe_lines).to_html(
                index=show_table_indices))
            if results.shape[0] > max_dataframe_lines:
                out.write(f'Omitted {results.shape[0]-max_dataframe_lines} '
                          'lines\n')
            out.write('<br>\n')


    out.write("""
    </body>
    </html>
    """)

    out.close()

