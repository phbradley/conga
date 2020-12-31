########################################################################################################################
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import sys
from sys import exit
import pandas as pd
import math
from sklearn.decomposition import PCA
from . import svg_basic
from . import util
from . import preprocess
from . import pmhc_scoring
from . import correlations
from . import tcr_scoring
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

def get_integers_color_dict( num_categories, cmap=plt.get_cmap('tab20') ):
    C = {}
    for i in range(num_categories):
        C[i] = cmap( float(i)/(num_categories-1))
    return C

def add_categorical_legend( ax, categories, colors, legend_fontsize=None, ncol=None ):
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


def make_single_rank_genes_logo( ranks, svgfile,
                                 margin = 10,
                                 logo_width = 300,
                                 logo_max_height = 400,
                                 top_pval_for_max_height = 1e-50,
                                 min_pval_for_scaling = 1e-300,
                                 num_genes_to_show = 10,
                                 signcolors=False,
                                 create_png = True,
                                 create_html = False ):
    cmds = []

    upper_left=  [margin, margin]
    cmds.extend( make_rank_genes_logo_stack( ranks, upper_left, logo_width, logo_max_height,
                                             top_pval_for_max_height=top_pval_for_max_height,
                                             min_pval_for_scaling=min_pval_for_scaling,
                                             signcolors=signcolors,
                                             num_genes_to_show=num_genes_to_show ) )

    svg_basic.create_file( cmds, logo_width+2*margin, logo_max_height + 2*margin, svgfile,
                           create_png=create_png, create_html=create_html, background_color='white' )

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
            if util.is_vdj_gene(gene, organism, include_constant_regions=True): # NEWNEWNEW
                continue
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
        clusters_gex= None,
        clusters_tcr= None,
        clusters_gex_names=None,
        clusters_tcr_names=None,
        ignore_tcr_cluster_colors = False,
        good_bicluster_tcr_scores=None,
        rank_genes_uns_tag = 'rank_genes_good_biclusters',
        include_alphadist_in_tcr_feature_logos=False,
        max_expn_for_gene_logo = 2.5, # or max over the clps, whichever is larger
        show_pmhc_info_in_logos = False,
        nocleanup = False, # dont delete temporary image files (useful for debugging)
        conga_scores = None,
        good_score_mask = None,
        make_batch_bars = None,
        batch_keys = None,
        make_cluster_gex_logos = True,
        draw_edges_between_conga_hits = True,
        add_conga_scores_colorbar = False,
        add_gex_logos_colorbar = False,
        spacey = False, # add extra space to make everything harder to read

        ## controls for the gene expression thumbnails that come before the actual logos:
        gex_header_genes=None,
        make_gex_header=True,
        make_gex_header_raw=True,
        make_gex_header_nbrZ=True,
        gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'], # was alphadist
        include_full_tcr_cluster_names_in_logo_lines=False,

):
    ''' need:
    * gex/tcr clusters: (obsm)
    * clones files for cluster-pairs with enough "good" clones
    * clone sizes (obs: clone_sizes)
    * tcrs (obs)
    * 2d projections: gex and tcr (obsm: X_gex_2d, X_tcr_2d)
    * conga scores (obs: conga_scores) UNLESS passed in
    * X_igex (obsm: X_igex) and X_igex_genes (uns)
    * rank genes info for each cluster-pair

    This code is a bit of a mess. OK, actually a huge mess. It started out simple but then when through
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
            clusters_gex_names = adata.uns.get('clusters_gex_names', [str(x) for x in range(np.max(clusters_gex)+1)])
    elif clusters_gex_names is None:
        clusters_gex_names = [str(x) for x in range(np.max(clusters_gex)+1)]

    if clusters_tcr is None:
        clusters_tcr = np.array(adata.obs['clusters_tcr'])
        if clusters_tcr_names is None:
            clusters_tcr_names = adata.uns.get('clusters_tcr_names', [str(x) for x in range(np.max(clusters_tcr)+1)])
    elif clusters_tcr_names is None:
        clusters_tcr_names = [str(x) for x in range(np.max(clusters_tcr)+1)]

    if good_score_mask is None:
        good_score_mask = np.array(list(adata.obs['good_score_mask']))
    X_gex_2d = adata.obsm['X_gex_2d']
    X_tcr_2d = adata.obsm['X_tcr_2d']
    if conga_scores is None:
        conga_scores = np.array(adata.obs['conga_scores'])
    X_igex = adata.obsm['X_igex']
    X_igex_genes = list(adata.uns['X_igex_genes']) #for .index
    organism = adata.uns['organism']

    if show_pmhc_info_in_logos:
        if 'X_pmhc' not in adata.obsm_keys():
            print('ERROR: include_pmhc_info_in_dendrogram=True but no X_pmhc info in adata.obsm_keys()')
            sys.exit()
        pmhc_var_names = adata.uns['pmhc_var_names']
        X_pmhc = adata.obsm['X_pmhc']
        assert X_pmhc.shape == ( adata.shape[0], len(pmhc_var_names))

    tcrs = preprocess.retrieve_tcrs_from_adata(adata)
    ################################################################# no more unpacking below here...

    num_clones = adata.shape[0]
    num_good_clones = np.sum(good_score_mask)
    assert X_igex.shape == (num_clones, len(X_igex_genes))


    if logo_genes is None:
        logo_genes = default_logo_genes[organism]

    if gex_header_genes is None:
        header2_genes = default_gex_header_genes[organism]
    else:
        header2_genes = gex_header_genes[:]
    if gex_header_tcr_score_names:
        header2_tcr_scores = tcr_scoring.make_tcr_score_table(adata, gex_header_tcr_score_names)
        assert header2_tcr_scores.shape == ( adata.shape[0], len(gex_header_tcr_score_names))
    else:
        header2_tcr_scores = None

    if 'clone_sizes' in header2_genes:
        X_igex = np.hstack( [X_igex, np.log(np.array(clone_sizes)[:,np.newaxis]) ] )
        X_igex_genes.append('clone_sizes')
        assert X_igex.shape == (num_clones, len(X_igex_genes))
    if 'nndists_gex' in header2_genes:
        if 'nndists_gex' in adata.obs_keys():
            X_igex = np.hstack( [X_igex, np.array(adata.obs['nndists_gex'])[:,np.newaxis] ] )
            X_igex_genes.append('nndists_gex')
            assert X_igex.shape == (num_clones, len(X_igex_genes))
        else:
            print('WARNING nndists_gex not found in adata.obs')
            header2_genes.remove('nndists_gex')

    gene_width = 6


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
        gex_nbrhood_X_igex.append( np.sum( X_igex[nbrhood_mask,:], axis=0 )/num )
        if header2_tcr_scores is not None:
            gex_nbrhood_tcr_scores.append( np.sum( header2_tcr_scores[nbrhood_mask,:], axis=0 )/num )

    gex_nbrhood_X_igex = np.array(gex_nbrhood_X_igex)
    assert gex_nbrhood_X_igex.shape == X_igex.shape
    if header2_tcr_scores is not None:
        gex_nbrhood_tcr_scores = np.array(gex_nbrhood_tcr_scores)
        assert gex_nbrhood_tcr_scores.shape == header2_tcr_scores.shape

    # set max expression levels for gex logo


    ## read info on louvains and tcr
    ## node2cluspair only includes good nodes
    node2cluspair = { i:(x,y) for i,(x,y,m) in enumerate(zip(clusters_gex, clusters_tcr, good_score_mask)) if m}
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
        all_dbl_nbrs[ii] = [ x for x in nbrs_tcr[ii] if good_score_mask[x] and x in gex_nbrs ]


    ## parse the X_igex matrix
    all_scores = {} # map from clp to [ num_clusters_tcr, fracs, means]


    logo_gene_indices = [ X_igex_genes.index(x) if x in X_igex_genes else None for x in logo_genes ]

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


    clps = [ x[1] for x in reversed( sorted( [ (y[0],x) for x,y in all_scores.items() ]))]

    # aim for 10in width
    dendro_width = 1.0
    single_batch_bar_width = 0.25
    batch_bars_width = 0 if not make_batch_bars else single_batch_bar_width * len(batch_keys)
    title_logo_width = 0.75
    rg_logo_width = 1.5
    score_logo_width = 0.0 if good_bicluster_tcr_scores is None else 0.5
    tcr_logo_width = 8
    logo_height = 0.5
    frac_scale = 130
    margin = 0.2 # left and right margins actually
    top_margin = 0.2
    bottom_margin = 1.3*margin
    yspacer_above_logos = 0.1
    yspacer_below_header = 0.15 if spacey else 0.0
    xspacer_within_header = 0.15 if spacey else 0.0

    gex_logo_width = logo_height * ( gene_width / (math.sqrt(3)+1) ) if make_cluster_gex_logos else 0.0

    fig_width = ( 2*margin + dendro_width + title_logo_width + batch_bars_width + rg_logo_width + tcr_logo_width +
                  score_logo_width + gex_logo_width )
    conga_scores_colorbar_width = 0.035*fig_width if add_conga_scores_colorbar else 0.0
    header_height = (fig_width-conga_scores_colorbar_width-xspacer_within_header-2*margin)/6.

    num_header2_rows = int(make_gex_header_nbrZ)+int(make_gex_header_raw)
    gex_logo_key_width = 2 if make_cluster_gex_logos else 0 # in header2 cols

    if make_gex_header:
        num_header2_tcr_score_cols = max(gex_logo_key_width,
                                         (gex_logo_key_width+len(gex_header_tcr_score_names))//num_header2_rows)
        gap_between_gex_and_tcr_panels = logo_height/10. if num_header2_tcr_score_cols else 0
        header2_height = (fig_width-2*margin-gap_between_gex_and_tcr_panels)/\
                         (len(header2_genes) + num_header2_tcr_score_cols)

    else:
        header2_height = 0.
    fig_height = ( top_margin + bottom_margin + header_height + num_header2_rows * header2_height +
                   yspacer_below_header + (len(clps)>0)*yspacer_above_logos + len(clps)*logo_height )

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


    # num_clusters_gex = int( popen('grep -c ^gex_louvain_chisq '+logfile).readlines()[0])
    # num_clusters_tcr = int( popen('grep -c ^tcr_louvain_chisq '+logfile).readlines()[0])

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
        if spacey:
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
        colors = np.array( [ max(0, -1*np.log10(conga_scores[x])) for x in ks ] )
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
        if spacey:
            plt.title('CoNGA scores',fontsize=9,pad=1)
            plt.xlabel(f'{proj_tag} UMAP1', fontsize=7, labelpad=1)
        else:
            plt.title('CoNGA scores ({}2D)'.format(proj_tag),fontsize=9,pad=1)
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
                plt.plot( xvals, yvals, zorder=2, marker='o', linestyle='None',
                          color=gex_colors[clp[0]], markeredgecolor='none', markeredgewidth=0 )
            else:
                for score_threshold, markersize in [ [1., 6], [1e-1, 7], [1e-2, 8] ]:
                    mask = cscores <= score_threshold
                    if np.sum(mask)==0:
                        continue # nothing to plot
                    print('score_threshold:', score_threshold, 'num:', np.sum(mask),
                          'markersize:', markersize)
                    reorder = np.argsort(cscores[mask])[::-1] # most significant last
                    plt.plot( xvals[mask][reorder], yvals[mask][reorder], zorder=2,
                              marker='o', linestyle='None', markersize=markersize,
                              color=gex_colors[clp[0]], markerfacecoloralt=tcr_colors[clp[1]],
                              fillstyle='left', markeredgecolor='none', markeredgewidth=0 )
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.ylabel('{} UMAP2'.format(TAG), fontsize=12,labelpad=1)
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        if spacey:
            plt.title('CoNGA hits', fontsize=9, pad=1) #pad=8)
            plt.xlabel(f'{proj_tag} UMAP1', fontsize=7, labelpad=1)
        else:
            plt.title('CoNGA hits ({}2D)'.format(proj_tag),fontsize=9,pad=1)

        if add_conga_scores_colorbar and icol==1:
            ax_width = conga_scores_colorbar_width/8.
            ax_height = 0.72*header_height
            bottom = (fig_height-top_margin-header_height+0.08*header_height)/fig_height
            height = ax_height/fig_height
            left = (margin+6*header_height+xspacer_within_header+0.5*ax_width)/fig_width
            width = ax_width/fig_width
            plt.axes( [left,bottom,width,height] )
            sm = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap='viridis')
            sm.set_array(np.sqrt( np.array( [ max(0, -1*np.log10(conga_scores[x])) for x in ks ] )))

            #plt.scatter([],[],c=[],cmap='viridis', vmin=vmin, vmax=vmax)
            #plt.gca().set_visible(False)
            cbar = plt.colorbar(mappable=sm, cax=plt.gca(),
                                ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
            # cbar = plt.colorbar(
            #     cax=plt.gca(), ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
            cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'], fontsize=7)
            #cbar.ax.set_ylabel('CoNGA score')
            cbar.ax.set_title('CoNGA\nscore', loc='left', ha='left', fontsize=9)



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
            #cbar.ax.set_ylabel('CoNGA score')
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
        plt.text(0.0,0.0, 'Clusters (size>{:d})'.format(int(min_cluster_size)-1),
                 ha='left', va='top', transform=plt.gca().transAxes)
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
        print('clp:',irow,len(leaves),logo_pngfile)
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
            plt.plot( [0.0], [0.0], marker='o', linestyle='None',
                      color=gex_colors[clp[0]], markeredgecolor='none', markeredgewidth=0, markersize=30)
        else:
            plt.plot( [0.0], [0.0], marker='o', linestyle='None',
                      color=gex_colors[clp[0]], markerfacecoloralt=tcr_colors[clp[1]],
                      fillstyle='left', markeredgecolor='none', markeredgewidth=0, markersize=30)
        # if gex_cluster_names is None:
        #     plt.text(-0.1,0,str(clp[0]),va='center',ha='right')
        # else:
        name = clusters_gex_names[clp[0]]
        short_name_fontsize, long_name_fontsize = 11, 6
        if len(name) <= 2:
            plt.text(-0.1,0,name,va='center',ha='right',fontsize=short_name_fontsize)
        else:
            plt.text(-0.1,0,name,va='center',ha='right',fontsize=long_name_fontsize,
                     bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=1))

        if not ignore_tcr_cluster_colors:
            if include_full_tcr_cluster_names_in_logo_lines:
                name = clusters_tcr_names[clp[1]]
                if len(name)<=2:
                    plt.text( 0.1, 0, name, va='center', ha='left', fontsize=short_name_fontsize)
                else:
                    plt.text( 0.1, 0, name, va='center', ha='left', fontsize=long_name_fontsize,
                              bbox=dict(facecolor='white', alpha=0.5, edgecolor='white', pad=1))
            else:
                plt.text( 0.1,0,str(clp[1]),va='center',ha='left', fontsize=short_name_fontsize)
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
                if num_batch_key_choices <= 10:
                    colors = plt.get_cmap('tab10').colors[:num_batch_key_choices]
                elif num_batch_key_choices <= 20:
                    colors = plt.get_cmap('tab20').colors[:num_batch_key_choices]
                else:
                    cmap = plt.get_cmap('tab20')
                    colors = [cmap(0.99*x/(num_batch_key_choices-1)) for x in range(num_batch_key_choices)]
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


        if make_cluster_gex_logos: # make the gex logo
            left = (margin+dendro_width+title_logo_width+batch_bars_width+rg_logo_width+tcr_logo_width+score_logo_width)/fig_width
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




    plt.savefig(logo_pngfile, dpi=300)

    if not nocleanup:
        for tmpfile in tmpfiles:
            if exists(tmpfile):
                os.remove(tmpfile)
            svgfile = tmpfile[:-3]+'svg'
            if exists(svgfile):
                os.remove(svgfile)



def make_n_pseudopoints( n, xy, radius_in_sdevs=0.25 ):
    ''' Find the centroid of the 2D xy array, pick n points around there
    '''
    assert xy.shape[1] == 2
    center = np.mean(xy, axis=0)
    sdev = np.sqrt( np.sum( np.square(xy-center[np.newaxis,:]))/xy.shape[0])
    radius = radius_in_sdevs * sdev
    print('make_n_pseudopoints:', n, sdev, radius, center)
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
):
    assert ax or pngfile

    if results_df.shape[0] == 0: # no results to plot
        print('plot_ranked_strings_on_cells: empty results')
        return


    # new wrinkle: clone_index column may be -1, if we included the cluster-level results
    # figure out if these are gex or tcr clusters
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
            if DONT_SHOW_LOW_MAIT_SCORES and gene=='mait' and 'gex' in xy_tag and row[4]<0:
                mask.append(False)
            elif pval > ranking_threshold or gene in exclude_strings or (ii!=-1 and ii in seen) or \
                 (gene,clp) in seen or clp_count[clp] >= strings_per_cluster:
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
    max_color_value = np.sqrt( -1*np.log10( max_color_pval/ranking_threshold) )

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

    if pngfile:
        plt.savefig(pngfile)



def make_summary_figure(
        adata,
        tcr_genes_results,
        gex_scores_results,
        pngfile,
        markersize=20,
        pval_threshold_for_tcr_genes_results=0.05, # but bonferroni is harder on tcrs since there are so many lame genes
        pval_threshold_for_gex_scores_results=0.05,
):
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
        add_categorical_legend( plt.gca(), cluster_names, [C[x] for x in range(num_clusters)] )
        if irow==0:
            plt.text(0.01,0.01,'{} clonotypes'.format(adata.shape[0]), ha='left', va='bottom', fontsize=8,
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
        cbar = plt.colorbar(ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
        cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'])
        cbar_axs.append(cbar.ax)
        #xmin,xmax = plt.xlim()
        #ymin,ymax = plt.ylim()
        plt.title('CoNGA scores')



        ## now a plot of the nbrhood gene/score enrichments
        if xy_tag == 'tcr':
            plot_ranked_strings_on_cells(adata, tcr_genes_results, 'X_tcr_2d', 'clone_index',
                                         'mwu_pvalue_adj', pval_threshold_for_tcr_genes_results, 'feature',
                                         ax=axs[irow,2] )
            plt.sca(axs[irow,2])
            plt.title("Differential gene expression in TCR nbrhoods")
            plt.xlabel('{} UMAP1'.format(XY_TAG))
            plt.ylabel('')
        else:
            plot_ranked_strings_on_cells(adata, gex_scores_results, 'X_gex_2d', 'clone_index',
                                         'mwu_pvalue_adj', pval_threshold_for_gex_scores_results, 'feature',
                                         direction_column='ttest_stat',
                                         ax=axs[irow,2])
            plt.title('TCR sequence feature bias in GEX nbrhoods')
            plt.xlabel('{} UMAP1'.format(XY_TAG))
            plt.ylabel('')

    #plt.tight_layout()
    plt.subplots_adjust(left=0.05, hspace=0.16, wspace=0.19, bottom=0.05, right = 0.97, top=0.97)
    # fiddle with sizes
    for irow in range(2):
        ax = axs[irow,1]
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 1.17, box.height])
        ax = cbar_axs[irow]
        box = ax.get_position()
        ax.set_position([box.x0+0.04, box.y0, box.width, box.height])
    plt.savefig(pngfile)





def make_clone_gex_umap_plots(
        adata,
        outfile_prefix,
        max_clones = 16,
        dpi=200,
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
        plt.scatter(xy[mask,0], xy[mask,1], s=16, c=clusters_gex[mask], #c='blue',
                    alpha=0.5)
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

    if 'batch_keys' not in adata.uns_keys():
        print('make_clone_batch_clustermaps: no batch_keys in adata.uns')
        return

    cmap_for_row_scores = plt.get_cmap(cmap_for_row_scores)

    num_clones = adata.shape[0]
    batch_keys = adata.uns['batch_keys']

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
        print('saving', pngfile)
        cm.savefig(pngfile, dpi=200)
        #assert False




def make_feature_panel_plots(
        adata,
        xy_tag,
        all_nbrs,
        results_df,
        pngfile,
        max_panels_per_bicluster=3,
        max_pvalue=0.05,
        nrows=6,
        ncols=4,
        use_nbr_frac=None
):
    '''
    '''
    xy = adata.obsm['X_{}_2d'.format(xy_tag)]

    # first let's figure out how many panels we could have
    # sort the results by pvalue
    df = results_df.sort_values('mwu_pvalue_adj')# makes a copy

    nbr_frac_for_cluster_results = max(all_nbrs.keys()) if use_nbr_frac is None else use_nbr_frac

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
        print(f)
        if f in var_names:
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
            feature_score_table = tcr_scoring.make_tcr_score_table(adata, [f])
            feature_to_raw_values[f] = feature_score_table[:,0]


    plt.figure(figsize=(ncols*3,nrows*3))
    plotno=0
    for row in df.itertuples():
        plotno+=1
        plt.subplot(nrows, ncols, plotno)

        feature = row.feature

        scores = feature_to_raw_values[feature]

        row_nbr_frac = use_nbr_frac if use_nbr_frac is not None else nbr_frac_for_cluster_results if row.nbr_frac==0.0 \
                       else row.nbr_frac
        nbrs = all_nbrs[row_nbr_frac][0] if xy_tag=='gex' else all_nbrs[row_nbr_frac][1]
        assert nbrs.shape[0] == adata.shape[0]
        num_neighbors = nbrs.shape[1] # this will not work for ragged nbr arrays (but we could change it to work)

        nbr_averaged_scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

        plt.scatter(xy[:,0], xy[:,1], c=nbr_averaged_scores, cmap='coolwarm', s=15)
        plt.title('{} ({:d},{:d}) {:.1e}'.format(feature, int(row.gex_cluster+.1), int(row.tcr_cluster+.1),
                                                 row.mwu_pvalue_adj))
        plt.xticks([],[])
        plt.yticks([],[])
        if (plotno-1)//ncols == nrows-1:
            plt.xlabel('{} UMAP1'.format(xy_tag.upper()))
        if plotno%ncols == 1:
            plt.ylabel('{} UMAP2'.format(xy_tag.upper()))

    plt.tight_layout()
    plt.savefig(pngfile)


def get_raw_feature_scores( feature, adata, feature_type=None):
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
    elif feature_type == 'gex' or (feature_type==None and feature in adata.raw.var_names):
        # will this be slow? creating list every time... probably OK for plotting routines
        return adata.raw.X[:, list(adata.raw.var_names).index(feature)].toarray()[:,0]
    else:
        assert feature_type in [ None, 'tcr']
        return tcr_scoring.make_tcr_score_table(adata,[feature])[:,0].astype(float)



def plot_hotspot_umap(
        adata,
        xy_tag,
        results_df,
        pngfile,
        nbrs = None,
        compute_nbr_averages=False,
        nrows=6,
        ncols=4,
):
    """
    xy_tag: use 'gex' or 'tcr' to set umap space used for plotting the hotspot features
    results_df : pandas df of hotspot features to plot. Expected columns are pvalue_adj feature feature_type
    where feature_type is either 'gex' or 'tcr'. We need that since some feature strings can be both. Output
    from correlations.find_hotspots can be fed in directly
    """
    if results_df.shape[0]==0:
        print('no results to plot:', pngfile)
        return

    df = results_df.sort_values('pvalue_adj') # ensure sorted
    df = df.iloc[:nrows*ncols,:]

    xy = adata.obsm['X_{}_2d'.format(xy_tag)]
    var_names = list(adata.raw.var_names)
    if df.shape[0] < nrows*ncols:
        # make fewer panels
        nrows = max(1, int(np.sqrt(df.shape[0])))
        ncols = (df.shape[0]-1)//nrows + 1

    plt.figure(figsize=(ncols*3,nrows*3))
    plotno=0
    for row in df.itertuples():
        plotno+=1
        plt.subplot(nrows, ncols, plotno)

        scores = get_raw_feature_scores( row.feature, adata , row.feature_type )
        # if row.feature_type == 'gex':
        #     if row.startswith
        #     assert row.feature in var_names
        #     scores = adata.raw.X[:, var_names.index(row.feature)].toarray()[:,0]
        # else:
        #     scores = tcr_scoring.make_tcr_score_table(adata,[row.feature])[:,0]
        if compute_nbr_averages:
            assert nbrs.shape[0] == adata.shape[0]
            num_neighbors = nbrs.shape[1] # this will not work for ragged nbr arrays (but we could change it to work)
            scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

        reorder = np.argsort(scores)
        plt.scatter(xy[reorder,0], xy[reorder,1], c=scores[reorder], cmap='viridis', s=15)
        plt.title('{} {:.1e}'.format(row.feature, row.pvalue_adj))
        plt.xticks([],[])
        plt.yticks([],[])
        if (plotno-1)//ncols == nrows-1:
            plt.xlabel('{} UMAP1'.format(xy_tag.upper()))
        if plotno%ncols == 1:
            plt.ylabel('{} UMAP2'.format(xy_tag.upper()))

    plt.tight_layout()
    plt.savefig(pngfile)



def plot_interesting_features_vs_clustermap(
        adata,
        features,
        pngfile,
        dist_tag, # gex or tcr
        feature_types = None, # list of either 'gex','tcr'
        nbrs = None,
        compute_nbr_averages=False,
        feature_labels = None,
        feature_scores = None,
        feature_scores_cmap = 'viridis',
        show_gex_cluster_colorbar = None,
        show_tcr_cluster_colorbar = None,
        show_VJ_gene_segment_colorbars = None,
        min_vmax = 0.45,
        max_vmax = 0.45,
        max_redundant_features = None,
        redundancy_threshold = 0.9, # correlation
        extra_feature_colors = None,
        rescale_factor_for_self_features = 0.33,
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


    assert dist_tag in ['gex','tcr']

    if show_gex_cluster_colorbar is None:
        show_gex_cluster_colorbar = dist_tag=='gex'

    if show_tcr_cluster_colorbar is None:
        show_tcr_cluster_colorbar = dist_tag=='tcr'

    if show_VJ_gene_segment_colorbars is None:
        show_VJ_gene_segment_colorbars = dist_tag=='tcr'

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

    # compute linkage of cells based on tcr
    X = adata.obsm['X_pca_'+dist_tag] # the kernal pca components
    print(f'computing pairwise X_pca_{dist_tag} distances, {X.shape}')
    Y = distance.pdist(X, metric='euclidean')

    print(f'computing linkage matrix from pairwise X_pca_{dist_tag} distances')
    cells_linkage = hierarchy.linkage(Y, method='ward')

    A = np.zeros((nrows, ncols))

    print('filling the score array for', len(features), 'features')
    for ii,feature in enumerate(features):
        scores = get_raw_feature_scores( feature, adata, feature_types if feature_types is None else feature_types[ii])
        # is_gex_feature = feature in var_names and ( feature_types is None or feature_types[ii]=='gex')
        # if is_gex_feature:
        #     scores = adata.raw.X[:, var_names.index(feature)].toarray()[:,0]
        # else:
        #     # could end up here if for some reason it's a gex feature but not present in var_names (shouldnt happen!)
        #     scores = tcr_scoring.make_tcr_score_table(adata, [feature])[:,0].astype(float)

        mn, std = np.mean(scores), np.std(scores)
        scores = (scores-mn)
        if std!=0:
            scores /= std

        if compute_nbr_averages:
            num_neighbors = nbrs.shape[1] # this will not work for ragged nbr arrays (but we could change it to work)
            scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

        if feature_types is not None:
            if feature_types[ii] == dist_tag:
                scores *= rescale_factor_for_self_features # lighten these up a bit since they will be nbr correlated

        A[ii,:] = scores

    tiny_lines = None

    if max_redundant_features is not None:
        feature_nbrs = {}

        # filter features by correlation
        # will subset: features, feature_types, feature_labels, A, nrows
        C = 1-distance.squareform(distance.pdist(A, metric='correlation'), force='tomatrix')
        feature_nbr_counts = [0]*len(features)
        feature_mask = np.full(len(features), True)
        for ii,f1 in enumerate(features):
            # am I too close to a previous feature?
            for jj in range(ii-1):
                if feature_types is None or feature_types[ii] == feature_types[jj]:
                    if C[ii,jj] > redundancy_threshold and feature_mask[jj]:
                        print('close:', C[ii,jj], f1, features[jj])
                        feature_nbr_counts[jj] += 1
                        if feature_nbr_counts[jj] > max_redundant_features:
                            feature_mask[ii] = False
                            feature_nbrs.setdefault(features[jj],[]).append(f1)
                            print('too many close:', feature_nbr_counts[jj], features[jj])
                            break
        if np.sum(feature_mask)<len(features): # have to exclude some
            # write out the feature neighbors for inspection
            dfl = []
            for f, nbrs in feature_nbrs.items():
                dfl.append(OrderedDict(feature=f, redundant_neighbors=','.join(nbrs), num_redundant_neighbors=len(nbrs)))
            df = pd.DataFrame(dfl)
            df.sort_values('num_redundant_neighbors', inplace=True, ascending=False)
            df.to_csv(pngfile+'_filtered_feature_info.tsv', sep='\t', index=False)
            tiny_lines = []
            for l in df.itertuples():
                tiny_lines.extend( [x+'\n' for x in [l.feature+':']+ l.redundant_neighbors.split(',')+[''] ] )

            #
            A = A[feature_mask,:]
            features = np.array(features)[feature_mask]
            if feature_types is not None:
                feature_types = np.array(feature_types)[feature_mask]
            new_feature_labels = []
            for old_label, nbr_count, m in zip(feature_labels, feature_nbr_counts, feature_mask):
                if m:
                    if nbr_count>max_redundant_features:
                        new_feature_labels.append('{} [+{}]'.format(old_label, nbr_count-max_redundant_features))
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

    if show_gex_cluster_colorbar:
        clusters_gex = np.array(adata.obs['clusters_gex'])
        num_clusters_gex = np.max(clusters_gex)+1
        cluster_color_dict = get_integers_color_dict(num_clusters_gex)
        colorbar_names.append( 'gex_cluster')
        colorbar_colors.append( [cluster_color_dict[x] for x in clusters_gex] )
        colorbar_sorted_tuples.append( [ ('gexclus'+str(x), cluster_color_dict[x]) for x in range(num_clusters_gex)])

    if show_tcr_cluster_colorbar:
        clusters_tcr = np.array(adata.obs['clusters_tcr'])
        num_clusters_tcr = np.max(clusters_tcr)+1
        cluster_color_dict = get_integers_color_dict(num_clusters_tcr)
        colorbar_names.append( 'tcr_cluster')
        colorbar_colors.append( [cluster_color_dict[x] for x in clusters_tcr] )
        colorbar_sorted_tuples.append( [ ('tcrclus'+str(x), cluster_color_dict[x]) for x in range(num_clusters_tcr)])

    if show_VJ_gene_segment_colorbars:
        tcrs = preprocess.retrieve_tcrs_from_adata(adata)
        gene_colors, sorted_gene_colors = assign_colors_to_conga_tcrs(tcrs, organism, return_sorted_color_tuples=True)
        colorbar_names.extend('VA JA VB JB'.split())
        colorbar_colors.extend(gene_colors)
        num_to_show = 10 # in the text written at the top
        colorbar_sorted_tuples.extend( [ x[:num_to_show] for x in sorted_gene_colors ] )

    # add some row colors
    colorbar_row_colors = []
    if feature_scores is not None:
        cm = plt.get_cmap(feature_scores_cmap)
        vmin, vmax = min(feature_scores), max(feature_scores)
        if vmax==vmin: vmax +=.001
        normed_scores = [ (x-vmin)/(vmax-vmin) for x in feature_scores]
        colorbar_row_colors.append( [ cm(x) for x in normed_scores])

    if feature_types is not None:
        type2color = {'tcr':'C0', 'gex':'C1'}
        colorbar_row_colors.append( [ type2color.get(x,'black') for x in feature_types] )

    if extra_feature_colors is not None:
        colorbar_row_colors.append(extra_feature_colors)

    if not colorbar_row_colors:
        colorbar_row_colors = None

    if feature_types is None:
        mx = np.max(np.abs(A))
    else:
        # try to only use other-type features for setting max
        mask = [ x=='gex' for x in feature_types]
        if np.sum(mask):
            mx = np.max(np.abs(A[mask,:]))
        else:
            mx = np.max(np.abs(A))

    mx = min(max_vmax, max(min_vmax, mx))

    print('making clustermap; num_features=', len(features))
    xmargin = 0.1
    ymargin = 0.1
    col_dendro_height = 2.5 #inches
    col_colors_height = len(colorbar_colors)*0.2
    heatmap_height = len(features)*0.15
    row_dendro_width = 2.5 #inches
    row_colors_width = 0.1 if colorbar_row_colors is None else len(colorbar_row_colors)*0.3
    heatmap_width = 7.5
    fig_width = 2*xmargin + row_dendro_width + row_colors_width + heatmap_width
    fig_height = 2*ymargin + col_dendro_height + col_colors_height + heatmap_height

    x0 = xmargin/fig_width
    x1 = (xmargin+row_dendro_width)/fig_width
    x2 = (xmargin+row_dendro_width+row_colors_width)/fig_width
    x3 = (fig_width-xmargin)/fig_width

    y0 = ymargin/fig_height
    y1 = (ymargin+heatmap_height)/fig_height
    y2 = (ymargin+heatmap_height+col_colors_height)/fig_height
    y3 = (fig_height-ymargin)/fig_height

    # not present in older version of seaborn:
    #dendrogram_inches = 2.
    #dendrogram_ratio = (dendrogram_inches/fig_height, dendrogram_inches/fig_width)
    cm = sns.clustermap(A, col_linkage=cells_linkage, metric='correlation', cmap='coolwarm', vmin=-mx, vmax=mx,
                        col_colors=colorbar_colors, row_colors=colorbar_row_colors,
                        figsize=(fig_width,fig_height))#, dendrogram_ratio=dendrogram_ratio)

    # print(dir(cm))
    # def show_pos(ax, tag):
    #     if ax is None:
    #         return
    #     pos = ax.get_position()
    #     print(f'axis lower_left: {pos.x0:7.3f} {pos.y0:7.3f} upper right: {pos.x1:7.3f} {pos.y1:7.3f}  {tag}')

    # show_pos(cm.ax_heatmap, 'heatmap')
    # show_pos(cm.ax_col_dendrogram, 'col_dendro')
    # show_pos(cm.ax_col_colors, 'col_colors')
    # show_pos(cm.ax_row_dendrogram, 'row_dendro')
    # show_pos(cm.ax_row_colors, 'row_colors')
    # show_pos(cm.cax, 'cax')

    cm.ax_row_dendrogram.set_position([x0,y0,x1-x0,y1-y0])
    if cm.ax_row_colors is not None:
        cm.ax_row_colors.set_position([x1,y0,x2-x1,y1-y0])
    cm.ax_heatmap.set_position([x2,y0,x3-x2,y1-y0])
    if cm.ax_col_colors is not None:
        cm.ax_col_colors.set_position([x2,y1,x3-x2,y2-y1])
    cm.ax_col_dendrogram.set_position([x2,y2,x3-x2,y3-y2])
    cm.cax.set_position([x0,(fig_height-ymargin-1.5)/fig_height,0.3/fig_width,1.5/fig_height])

    row_indices = cm.dendrogram_row.reordered_ind
    row_labels = [feature_labels[x] for x in row_indices]
    cm.ax_heatmap.set_yticks( 0.5+np.arange(nrows))
    cm.ax_heatmap.set_yticklabels( row_labels, rotation='horizontal', fontdict={'fontsize':8} )
    cm.ax_heatmap.set_xticks([])
    cm.ax_heatmap.set_xticklabels([])

    # show the top couple genes and their colors
    ax = cm.ax_col_dendrogram
    col_dendro_width = fig_width * (ax.get_position().x1 - ax.get_position().x0 )
    col_dendro_height = fig_height * (ax.get_position().y1 - ax.get_position().y0 )

    line_height_inches = 0.09
    fontsize=7
    for ii, sorted_tuples in enumerate(colorbar_sorted_tuples):
        ax.text(-0.01, 1.0 - (ii*line_height_inches+0.05)/col_dendro_height,
                colorbar_names[ii], va='top', ha='right', fontsize=fontsize, color='k', transform=ax.transAxes)
        for jj in range(len(sorted_tuples)):
            gene, color = sorted_tuples[jj]
            ax.text(float(jj)/len(sorted_tuples), 1.0 - (ii*line_height_inches+0.05)/col_dendro_height,
                    gene, va='top', ha='left', fontsize=fontsize, color=color, transform=ax.transAxes)

    # show the tiny_lines
    if tiny_lines is not None:
        ax = cm.ax_row_dendrogram
        row_dendro_width = fig_width * (ax.get_position().x1 - ax.get_position().x0 )
        row_dendro_height = fig_height * (ax.get_position().y1 - ax.get_position().y0 )

        max_lines = int(row_dendro_height / 0.08)
        col_width_frac = 0.33 / row_dendro_width
        x_offset = 0
        tiny_fontsize=5
        ax.text(0.0, 1.0, 'Correlated\nfeatures\nomitted:', va='bottom', ha='left', fontsize=tiny_fontsize,
                transform=ax.transAxes, zorder=-1, color='black')
        while tiny_lines:
            ax.text(x_offset, 1.0, ''.join(tiny_lines[:max_lines]), va='top', ha='left', fontsize=tiny_fontsize,
                    transform=ax.transAxes, zorder=-1, color='blue')
            x_offset += col_width_frac
            tiny_lines = tiny_lines[max_lines:]
            if tiny_lines and x_offset+col_width_frac/2>1.0:
                ax.text(x_offset, 1.0, '...', va='top', ha='left', fontsize=tiny_fontsize,transform=ax.transAxes)
                break

    print('saving', pngfile)
    cm.savefig(pngfile, dpi=200)


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

    plt.figure(figsize=(10,4))

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
        **kwargs # passed to make_logo_plots
):
    cluster_tcr_score_names = tcr_scoring.all_tcr_scorenames

    good_mask = (scores <= max_good_score)

    bic_counts = Counter( (x,y) for x,y,m in zip(clusters_gex, clusters_tcr, good_mask) if m )

    num_good_biclusters = sum( 1 for x,y in bic_counts.items() if y>=min_cluster_size )

    if num_good_biclusters:
        # calc tcr sequence features of good cluster pairs
        good_bicluster_tcr_scores = correlations.calc_good_cluster_tcr_features(
            adata, good_mask, clusters_gex, clusters_tcr, cluster_tcr_score_names, min_count=min_cluster_size)

        # run rank_genes on most common bics
        if 'rank_genes_uns_tag' in kwargs:
            rank_genes_uns_tag = kwargs['rank_genes_uns_tag']
            del kwargs['rank_genes_uns_tag']# since we pass it explicitly below
        else:
            rank_genes_uns_tag = f'rg_{num_good_biclusters}_{np.sum(good_mask)}_biclusters'

        correlations.run_rank_genes_on_good_biclusters(
            adata, good_mask, clusters_gex, clusters_tcr, min_count=min_cluster_size, key_added= rank_genes_uns_tag)

        make_logo_plots(
            adata, nbrs_gex, nbrs_tcr, min_cluster_size, pngfile,
            conga_scores = scores,
            clusters_gex = clusters_gex,
            clusters_tcr = clusters_tcr,
            good_score_mask = good_mask,
            good_bicluster_tcr_scores = good_bicluster_tcr_scores,
            rank_genes_uns_tag = rank_genes_uns_tag,
            **kwargs)


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
    #unpacking
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    num_clones = adata.shape[0]

    #
    nbrhood_pvals = { 'gex':np.full((num_clones,), num_clones).astype(float),
                      'tcr':np.full((num_clones,), num_clones).astype(float) }
    for l in hotspot_nbrhood_results_df.itertuples():
        assert l.feature_type[3:] == '_nbrs_vs_graph'
        tag = l.feature_type[:3]
        nbrhood_pvals[tag][l.clone_index] = min(l.pvalue_adj, nbrhood_pvals[tag][l.clone_index])

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
            kwargs['rank_genes_uns_tag'] = f'rg_hotspot_{tag}_nbrhood_biclusters'

        make_cluster_logo_plots_figure(adata, pvals, pvalue_threshold, fake_clusters_gex, fake_clusters_tcr,
                                       nbrs_gex, nbrs_tcr, min_cluster_size, pngfile, **kwargs)



def make_tcr_clumping_plots(
        adata,
        results_df, # generated by tcr_clumping.py: assess_tcr_clumping
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size_for_logos,
        pvalue_threshold_for_logos,
        outfile_prefix,
        max_color_pvalue = 1e-16, # in the UMAP figure; no brighter beyond there
        **logo_plot_args,
):
    num_clones = adata.shape[0]
    fake_clusters_gex = np.zeros((num_clones,)).astype(int)
    fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
    clumping_pvals = np.full( (num_clones,), num_clones).astype(float)
    for l in results_df.itertuples():
        clumping_pvals[ l.clone_index] = min( l.pvalue_adj, clumping_pvals[l.clone_index])
        fake_clusters_tcr[l.clone_index] = l.clumping_group

    pngfile = f'{outfile_prefix}_tcr_clumping_logos.png'

    if 'rank_genes_uns_tag' not in logo_plot_args:
        logo_plot_args['rank_genes_uns_tag'] = 'rg_tcr_clumping_biclusters'

    make_cluster_logo_plots_figure(
        adata, clumping_pvals, pvalue_threshold_for_logos, fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile, **logo_plot_args)

    # make umaps colored by clumping pvals
    plt.figure(figsize=(12,6))
    for icol, tag in enumerate(['gex','tcr']):
        plt.subplot(1,2,icol+1)
        colors = np.sqrt( np.maximum(0.0, -1*np.log10(np.maximum(1e-100, clumping_pvals)))) # no log10 of 0.0
        #print('colors:', tag, np.max(colors), list(colors[:100]))
        reorder = np.argsort(colors)
        xy = adata.obsm['X_{}_2d'.format(tag)] # same umap as feature nbr-type
        vmax = np.sqrt(-1*np.log10(max_color_pvalue))
        #vmax = np.sqrt(-1*np.log10(1e-5))
        plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=vmax)
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
        batch_bias_pvals[ l.clone_index] = min( l.pvalue_adj, batch_bias_pvals[l.clone_index])
        fake_clusters_tcr[l.clone_index] = l.cluster_group

    pngfile = f'{outfile_prefix}_batch_bias_logos.png'

    make_cluster_logo_plots_figure(
        adata, batch_bias_pvals, pvalue_threshold_for_logos, fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
        min_cluster_size_for_logos, pngfile, rank_genes_uns_tag='rg_batch_bias_biclusters')

    # make umaps colored by batch_bias pvals
    plt.figure(figsize=(12,6))
    for icol, tag in enumerate(['gex','tcr']):
        plt.subplot(1,2,icol+1)
        colors = np.sqrt( np.maximum(0.0, -1*np.log10(np.maximum(1e-100, batch_bias_pvals)))) # no log10 of 0.0
        #print('colors:', tag, np.max(colors), list(colors[:100]))
        reorder = np.argsort(colors)
        xy = adata.obsm['X_{}_2d'.format(tag)] # same umap as feature nbr-type
        vmax = np.sqrt(-1*np.log10(max_color_pvalue))
        #vmax = np.sqrt(-1*np.log10(1e-5))
        plt.scatter( xy[reorder,0], xy[reorder,1], c=colors[reorder], vmin=0, vmax=vmax)
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



def make_tcrdist_trees( adata, output_prefix, group_by = None):
    """
    generate TCRdist trees by gene expression cluster
    group_by: use 'clusters_gex' or 'clusters_tcr' to generate by trees by respective cluster assignments.
    """

    if group_by is None or group_by == 'clusters_gex':
        group_by = 'clusters_gex'
        tag = 'GEX'
    elif group_by == 'clusters_tcr':
        group_by = 'clusters_tcr'
        tag = 'TCR'

    width = 800
    height = 1000
    xpad = 25
    organism = adata.uns['organism']
    precomputed = False
    clusters = np.array(adata.obs[group_by])

    num_clusters = np.max(clusters)+1
    tcrs = preprocess.retrieve_tcrs_from_adata(adata)

    num_clones = adata.shape[0]
    if 'conga_scores' in adata.obs_keys():
        conga_scores = np.maximum( 1e-100, np.array(adata.obs['conga_scores']) ) # no zeros!
        scores = np.sqrt( np.maximum( 0.0, -1*np.log10( 100*conga_scores/num_clones)))
    else:
        scores = np.zeros((adata.shape[0],))

    tcrdist = TcrDistCalculator(organism)

    x_offset = 0
    all_cmds = []

    #color_score_range = [-1*np.log(10), -1*np.log(1e-5)]
    color_score_range = [0, 3.0]
    print('color_score_range:', color_score_range)

    for clust in range(num_clusters):
        cmask = (clusters==clust)
        csize = np.sum(cmask)
        #cinds = np.nonzero(cmask)[0]

        ctcrs   = [x for x,y in zip(  tcrs, cmask) if y]
        cscores = [x for x,y in zip(scores, cmask) if y]

        if not precomputed:
            print('computing tcrdist distances:', clust, csize)
            cdists = np.array([ tcrdist(x,y) for x in ctcrs for y in ctcrs]).reshape(csize,csize)
        else:
            assert False # tmp hack

        cmds = make_tcr_tree_svg_commands(
            ctcrs, organism, [x_offset,0], [width,height], cdists, max_tcrs_for_trees=400, tcrdist_calculator=tcrdist,
            color_scores=cscores, color_score_range = color_score_range, title=f'{tag} cluster {clust}')

        x_offset += width + xpad

        all_cmds.extend(cmds)

    svgfile = f'{output_prefix}_{tag}_cluster_tcrdist_trees.svg'
    print('making:', svgfile[:-3]+'png')
    svg_basic.create_file(all_cmds, x_offset-xpad, height, svgfile, create_png= True )
