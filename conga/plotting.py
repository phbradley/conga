import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
from sys import exit
import pandas as pd
import math
from . import preprocess as pp
from . import svg_basic
from . import util
from . import tcr_scoring
from .tcrdist.make_tcr_logo import make_tcr_logo_for_tcrs
from .tcrdist.tcr_distances import TcrDistCalculator
import os
from os.path import exists
from collections import Counter
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
              'Ccl5', 'Cxcr3', 'Zbtb16', 'Nkg7', 'Klrd1']
}


default_gex_header_genes = {
    'human': ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
    'mouse': ['clone_sizes','Cd4', 'Cd8a', 'Cd8b1', 'Sell', 'Itgal', 'Gzma', 'Ccl5', 'Il2rb', 'Ikzf2', 'Pdcd1', 'Zbtb16']
}


# this is a bit of a hack
DONT_SHOW_LOW_MAIT_SCORES = True

def get_integers_color_dict( num_categories, cmap=plt.get_cmap('tab20') ):
    C = {}
    for i in range(num_categories):
        C[i] = cmap( float(i)/(num_categories-1))
    return C

def add_categorical_legend( ax, categories, colors, legend_fontsize=None ):
    for idx, label in enumerate(categories):
        color = colors[idx]
        # use empty scatter to set labels
        ax.scatter([], [], c=[color], label=label)
    ax.legend(
        frameon=False, loc='center left',
        bbox_to_anchor=(0.98, 0.5),
        ncol=(1 if len(categories) <= 20 # was 14
              else 2 if len(categories) <= 30 else 3),
        fontsize=legend_fontsize,
        borderpad=0.0,
        handletextpad=0.0)
    # Shrink current axis by 10% to fit legend and match
    # size of plots that are not categorical
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.91, box.height])

def add_integers_legend( ax, color_dict ):
    categories = [ str(x) for x in sorted( color_dict.keys() ) ]
    colors = [ color_dict[int(x)] for x in categories ]
    add_categorical_legend( ax, categories, colors )

def make_rank_genes_logo_stack( ranks, upper_left, logo_width, max_logo_height,
                                top_pval_for_max_height = 1e-30,
                                min_pval_for_scaling = 1e-300,
                                signcolors=False,
                                num_genes_to_show = 10 ):
    ''' ranks is a list of (gene,l2r,pval)
    '''

    def pval_factor( pval, min_pval=top_pval_for_max_height ):
        return math.sqrt( max(1e-3, -1 * math.log10( max(min_pval,pval) ) ))
        #return -1 * math.log10( max(min_pval,pval) )

    top_pval = ranks[0][2]
    logo_height = max_logo_height * pval_factor(top_pval) / pval_factor(top_pval_for_max_height)

    if logo_height<1e-3:
        return []

    height_factors = [ pval_factor( x[2], min_pval_for_scaling ) for x in ranks[:num_genes_to_show] ]

    height_scale = logo_height / sum(height_factors)

    x0,y_offset = upper_left
    y_offset += ( max_logo_height - logo_height )

    cmds = []
    for gene,l2r,pval in ranks[:num_genes_to_show]:
        height = height_scale * pval_factor(pval,min_pval_for_scaling)
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

def _parse_clusterpair_rank_genes( adata, uns_tag = 'rank_genes_good_cluster_pairs' ):
    ''' returns all_ranks dict mapping from cluspair to list [ (gene, l2r, pval), ... ]
    '''

    names = pd.DataFrame(adata.uns[uns_tag]['names'])

    all_ranks = {}
    for clptag in names.columns:
        assert clptag.startswith('clp_') and clptag.count('_') == 2
        clp = ( int(clptag.split('_')[1]), int(clptag.split('_')[2]) ) # gex_cluster, tcr_cluster
        ranks = []
        for igene, gene in enumerate( adata.uns[uns_tag]['names'][clptag] ):
            if pd.isnull(gene):
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
        good_cluster_pair_tcr_scores=None,
        rank_genes_uns_tag = 'rank_genes_good_cluster_pairs',
        include_alphadist_in_tcr_feature_logos=False,
        max_expn_for_gene_logo = 2.5, # or max over the clps, whichever is larger
        show_pmhc_info_in_logos = False,

        ## controls for the gene expression thumbnails that come before the actual logos:
        gex_header_genes=None,
        make_gex_header=True,
        make_gex_header_raw=True,
        make_gex_header_nbrZ=True,
        gex_header_tcr_score_names = ['mhci2', 'cdr3len', 'cd8', 'nndists_tcr'], # was alphadist
        include_full_tcr_cluster_names_in_logo_lines=False,

):
    ''' need:
    * gex/tcr clusters: (obsm)
    * clones files for cluster-pairs with enough "good" clones
    * clone sizes (obs: clone_sizes)
    * tcrs (obs)
    * 2d projections: gex and tcr (obsm: X_gex_2d, X_tcr_2d)
    * conga scores (aka all_max_nbr_chisq): (obsm: conga_scores)
    * X_igex (obsm: X_igex) and X_igex_genes (uns)
    * rank genes info for each cluster-pair
    '''

    if not make_gex_header:
        make_gex_header_raw, make_gex_header_nbrZ = False, False

    ## unpack data from adata arrays ##################################
    clone_sizes = adata.obs['clone_sizes']
    if clusters_gex is None:
        clusters_gex = adata.obs['clusters_gex']
        if clusters_gex_names is None:
            clusters_gex_names = adata.uns.get('clusters_gex_names', [str(x) for x in range(np.max(clusters_gex)+1)])
    elif clusters_gex_names is None:
        clusters_gex_names = [str(x) for x in range(np.max(clusters_gex)+1)]

    if clusters_tcr is None:
        clusters_tcr = adata.obs['clusters_tcr']
        if clusters_tcr_names is None:
            clusters_tcr_names = adata.uns.get('clusters_tcr_names', [str(x) for x in range(np.max(clusters_tcr)+1)])
    elif clusters_tcr_names is None:
        clusters_tcr_names = [str(x) for x in range(np.max(clusters_tcr)+1)]

    good_score_mask = np.array(list(adata.obs['good_score_mask']))
    X_gex_2d = adata.obsm['X_gex_2d']
    X_tcr_2d = adata.obsm['X_tcr_2d']
    conga_scores = adata.obsm['conga_scores'] # (num_clones,3)
    X_igex = adata.obsm['X_igex']
    X_igex_genes = list(adata.uns['X_igex_genes']) #for .index
    organism = adata.uns['organism']
    # nndists_gex = np.array(adata.obs['nndists_gex']) # for choosing a representative clonotype
    # nndists_tcr = np.array(adata.obs['nndists_tcr']) # to label

    if show_pmhc_info_in_logos:
        if 'X_pmhc' not in adata.obsm_keys():
            print('ERROR: include_pmhc_info_in_dendrogram=True but no X_pmhc info in adata.obsm_keys()')
            sys.exit()
        pmhc_var_names = adata.uns['pmhc_var_names']
        X_pmhc = adata.obsm['X_pmhc']
        assert X_pmhc.shape == ( adata.shape[0], len(pmhc_var_names))

    tcrs = pp.retrieve_tcrs_from_adata(adata)
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
        X_igex = np.hstack( [X_igex, np.array(adata.obs['nndists_gex'])[:,np.newaxis] ] )
        X_igex_genes.append('nndists_gex')
        assert X_igex.shape == (num_clones, len(X_igex_genes))

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
    all_max_nbr_chisq = { i:np.max(x) for i,x in enumerate(conga_scores) }


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
        means = np.array([ 0 if x is None else means[x] for x in logo_gene_indices ])
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
        all_ranks = _parse_clusterpair_rank_genes( adata, uns_tag=rank_genes_uns_tag)


    clps = [ x[1] for x in reversed( sorted( [ (y[0],x) for x,y in all_scores.items() ]))]

    # aim for 10in width
    dendro_width = 1.0
    title_logo_width = 0.75
    rg_logo_width = 1.5
    score_logo_width = 0.0 if good_cluster_pair_tcr_scores is None else 0.5
    tcr_logo_width = 8
    gex_logo_width = 2
    logo_height = 0.5
    frac_scale = 130
    margin = 0.2 # left and right margins actually
    top_margin = 0.2
    bottom_margin = 1.3*margin

    gex_logo_width = logo_height * ( gene_width / (math.sqrt(3)+1) )

    fig_width = 2*margin + dendro_width + title_logo_width + rg_logo_width + tcr_logo_width + score_logo_width + \
                gex_logo_width
    header_height = (fig_width-2*margin)/6.

    num_header2_rows = int(make_gex_header_nbrZ)+int(make_gex_header_raw)
    gex_logo_key_width = 2 # in header2 cols

    if make_gex_header:
        num_header2_tcr_score_cols = max(gex_logo_key_width,
                                         (gex_logo_key_width+len(gex_header_tcr_score_names))//num_header2_rows)
        gap_between_gex_and_tcr_panels = logo_height/10. if num_header2_tcr_score_cols else 0
        header2_height = (fig_width-2*margin-gap_between_gex_and_tcr_panels)/\
                         (len(header2_genes) + num_header2_tcr_score_cols)

    else:
        header2_height = 0.
    fig_height = top_margin + bottom_margin + header_height + num_header2_rows * header2_height + len(clps)*logo_height

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

    # store distances between clusterpairs
    D = np.full( (len(clps),len(clps)), 1.0 ) # dist is 1.0 if no edges between clusterpair
    for i in range(len(clps)):
        D[i,i] = 0.0

    ##############################################
    # make a header: 6 square scatter plots across the top
    # left to right:
    # 1 GEX UMAP, colored by GEX clusters
    # 2 TCR UMAP, colored by TCR clusters
    # 3 GEX UMAP, colored by conga score
    # 4 GEX UMAP, only good cells, colored by cluster-pairs
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
        left = (margin+3*icol*header_height)/fig_width # make these plots square, so height=width
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
        plt.title('{} clusters ({}2D)'.format(clusters_tag, proj_tag), fontsize=9, pad=1) #pad=8)
        if icol==0:
            plt.text(0.01,0.01,'{} clones'.format(len(all_xy)), ha='left', va='bottom', fontsize=8,
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
        left = (margin+(3*icol+1)*header_height)/fig_width # make these plots square, so height=width
        width = header_height/fig_width
        plt.axes( [left,bottom,width,height] )

        ## make a plot colored by pval
        colors = np.array( [ max(0,all_max_nbr_chisq[x]) for x in ks ] )
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
        #cbar = plt.colorbar(ticks=[ math.sqrt(-1*math.log10(x)) for x in [1.0,1e-1,1e-2,1e-3,1e-4,1e-5] ])
        #cbar.ax.set_yticklabels(['1','0.1','0.01','0.001','1e-4','1e-5'])
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        plt.title('CoNGA scores ({}2D)'.format(proj_tag),fontsize=9,pad=1)
        #plt.title('kNN overlap E-values {}-UMAP'.format('GEX' if icol==0 else 'TCR'),fontsize=12,pad=1)


        ##############################################
        ## now just showing the top scoring clones
        bottom = (fig_height-top_margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+(3*icol+2)*header_height)/fig_width # make these plots square, so height=width
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

        # the inter-cluspair lines
        light_gray = '#D3D3D3'
        line_segments = LineCollection( inter_lines, linewidths=0.25, zorder=1, colors=light_gray, alpha=0.25)
        plt.gca().add_collection(line_segments)

        # the intra-cluspair lines
        dark_gray = '#888888'
        line_segments = LineCollection( intra_lines, linewidths=0.25, zorder=1, colors='k')#dark_gray, alpha=0.5)
        plt.gca().add_collection(line_segments)


        for clp,nodes in cluspair2nodes.items():
            xvals = [ all_xy[x][0] for x in nodes ]
            yvals = [ all_xy[x][1] for x in nodes ]
            if ignore_tcr_cluster_colors:
                plt.plot( xvals, yvals, zorder=2, marker='o', linestyle='None',
                          color=gex_colors[clp[0]], markeredgecolor='none', markeredgewidth=0 )
            else:
                plt.plot( xvals, yvals, zorder=2, marker='o', linestyle='None',
                          color=gex_colors[clp[0]], markerfacecoloralt=tcr_colors[clp[1]],
                          fillstyle='left', markeredgecolor='none', markeredgewidth=0 )
        plt.xticks([],[])
        plt.yticks([],[])
        #plt.ylabel('{} UMAP2'.format(TAG), fontsize=12,labelpad=1)
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        plt.title('Top-scoring clones ({}2D)'.format(proj_tag),fontsize=9,pad=1)


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
                bottom = (fig_height-top_margin-header_height-header2_height)/fig_height
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
                plt.text( xmx, ymx, gene, va='top', ha='right', fontsize=8)

            if make_gex_header_nbrZ: ## make another plot with the gex-nbrhood averaged scores
                bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height)/fig_height
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
                plt.text( xmn+pad, ymx, 'nbrZ', va='top', ha='left', fontsize=8)
                plt.text( xmx, ymx, gene, va='top', ha='right', fontsize=8)

        if gex_header_tcr_score_names:
            for ii, scoretag in enumerate(gex_header_tcr_score_names):
                irow=ii//num_header2_tcr_score_cols
                icol=ii%num_header2_tcr_score_cols
                ## make another plot with the gex-nbrhood averaged TCR scores
                bottom = (fig_height-top_margin-header_height-(irow+1)*header2_height)/fig_height
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
                          scoretag.replace('alphadist','ad').replace('nndists_tcr','NNd').replace('cdr3len','len'),
                          va='top', ha='right', fontsize=8)
                plt.text( xmn+pad, ymx, 'tcr', va='top', ha='left', fontsize=8)


    if len(clps)>1:

        if True: ## make a key for the gex logo
            bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height)/fig_height
            height = header2_height/fig_height
            left = (fig_width-margin-2*header2_height)/fig_width # make these plots square, so height=width
            width = 2*header2_height/fig_width
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
            plt.scatter(xy[:,0], xy[:,1], c=color, s=350, zorder=1)
            for gene,(x,y) in zip( logo_genes, xy ):
                plt.text(x,y,gene,fontsize=8,ha='center',va='center',rotation=30, zorder=2)
            plt.xlim( [-0.5, gene_width-0.5] )
            plt.ylim( [-yshift-0.5, yshift+0.5 ] )
            plt.xticks([],[])
            plt.yticks([],[])
            #plt.title('key for GEX logo',fontsize=8,va='top',pad=1)

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
        plt.text(0.0,0.0, 'Biclusters (size>{:d})'.format(min_cluster_size-1),
                 ha='left', va='top', transform=plt.gca().transAxes)
        leaves = R['leaves'][:] #list( hierarchy.leaves_list( Z ) )
        leaves.reverse() # since we are drawing them downward, but the leaf-order increases upward
    elif clps:
        leaves = [0]
    else:
        leaves = []

    tmpfiles = []
    for irow, iclp in enumerate(leaves):
        last_row = ( irow == len(leaves)-1 )
        print('clp:',irow,len(leaves),logo_pngfile)
        clp = clps[iclp]
        bottom = (fig_height-top_margin-header_height-num_header2_rows*header2_height-(irow+1)*logo_height)/fig_height
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

        # make the rank genes logo
        if all_ranks is not None:
            clp_rank_genes = all_ranks[clp]
            for r in range(3):
                pngfile = '{}.tmp_{}_{}_{}.png'.format(logo_pngfile, clp[0], clp[1], r)
                start_rank, stop_rank = 3*r, 3*r+3
                left = (margin+dendro_width+title_logo_width+r*0.333*rg_logo_width)/fig_width
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
            left = (margin+dendro_width+title_logo_width+rg_logo_width+iab*0.5*tcr_logo_width)/fig_width
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
        if good_cluster_pair_tcr_scores is not None:
            clp_rank_scores = good_cluster_pair_tcr_scores[clp]
            if DONT_SHOW_LOW_MAIT_SCORES:
                clp_rank_scores = [ x for x in clp_rank_scores if x[0] != 'mait' or x[1]>0 ]
            if not include_alphadist_in_tcr_feature_logos:
                clp_rank_scores = [ x for x in clp_rank_scores if x[0] != 'alphadist' ]
            pngfile = '{}.tmpsc_{}_{}.png'.format(logo_pngfile, clp[0], clp[1])
            left = (margin+dendro_width+title_logo_width+rg_logo_width+tcr_logo_width)/fig_width
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
                plt.text(1.0, -0.05, 'TCRseq\nfeatures', ha='right', va='top', fontsize=8,
                         transform=plt.gca().transAxes)


        # make the gex logo
        left = (margin+dendro_width+title_logo_width+rg_logo_width+tcr_logo_width+score_logo_width)/fig_width
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
            plt.text(0.5,-0.05, 'GEX logo', ha='center', va='top', transform=plt.gca().transAxes)
            # msg = 'GEX logo genes: row1= {} row2= {} row3= {}'\
            #     .format(','.join(logo_genes[:gene_width-1]),
            #             ','.join(logo_genes[gene_width-1:2*gene_width-1]),
            #             ','.join(logo_genes[2*gene_width-1:]) )
            # plt.text(1.0,-0.5, msg, ha='right', va='top', transform=plt.gca().transAxes, fontsize=7)




    plt.savefig(logo_pngfile, dpi=300)


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
        clusters = adata.obs['clusters_gex']
        clp_index_for_minus1 = 0
    else:
        assert 'tcr' in xy_tag
        clusters = adata.obs['clusters_tcr']
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
    conga_scores = adata.obsm['conga_scores'] # (num_clones,3)

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
            plt.text(0.01,0.01,'{} clones'.format(adata.shape[0]), ha='left', va='bottom', fontsize=8,
                     transform=plt.gca().transAxes)


        ############################
        ## now a plot colored by pval
        plt.sca(axs[irow,1])
        colors = np.maximum( conga_scores.max(axis=1), 0.0 )
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
            exclude_strings = ['5830405F06Rik'] # bad mouse gene, actually a tcr v gene
            plot_ranked_strings_on_cells(adata, tcr_genes_results, 'X_tcr_2d', 'clone_index',
                                         'mwu_pvalue_adj', pval_threshold_for_tcr_genes_results, 'feature',
                                         exclude_strings=exclude_strings, ax=axs[irow,2] )
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





def make_clone_plots(adata, num_clones_to_plot, pngfile):
    ''' This is called before we've condensed to a single cell per clone
    So we don't have PCA or UMAP yet
    '''


    print('make_clone_plots: cluster_and_tsne_and_umap')
    adata = pp.cluster_and_tsne_and_umap( adata, skip_tcr=True )

    tcrs = pp.retrieve_tcrs_from_adata(adata)

    unique_tcrs = sorted(set(tcrs))

    tcr2clone_id = { x:i for i,x in enumerate(unique_tcrs)}
    clone_ids = np.array( [ tcr2clone_id[x] for x in tcrs ])

    counts = Counter(clone_ids)

    nrows = int(np.sqrt(num_clones_to_plot))
    ncols = (num_clones_to_plot-1)//nrows + 1

    plt.figure(figsize=(3*ncols,3*nrows))
    plotno=0

    xy = adata.obsm['X_gex_2d']
    for (clone_id, clone_size) in counts.most_common(num_clones_to_plot):
        plotno += 1
        plt.subplot(nrows, ncols, plotno)
        plt.scatter(xy[:,0], xy[:,1], s=16, c='gray',alpha=0.2)
        mask = clone_ids==clone_id
        plt.scatter(xy[mask,0], xy[mask,1], s=16, c='blue', alpha=0.5)
        plt.xlabel('GEX UMAP1')
        plt.ylabel('GEX UMAP2')
        plt.xticks([], [])
        plt.yticks([], [])
        plt.text(0, 0, '{} cells'.format(clone_size), ha='left', va='bottom', transform=plt.gca().transAxes)
    plt.tight_layout()
    print('making:', pngfile)
    plt.savefig(pngfile)





def make_feature_panel_plots(
        adata,
        xy_tag,
        all_nbrs,
        results_df,
        pngfile,
        max_panels_per_cluster_pair=3,
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
        if clp_counts[clp] >= max_panels_per_cluster_pair:
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

        if use_nbr_frac is None:
            nbrs = all_nbrs[row.nbr_frac][0] if xy_tag=='gex' else all_nbrs[row.nbr_frac][1]
        else:
            nbrs = all_nbrs[use_nbr_frac][0] if xy_tag=='gex' else all_nbrs[use_nbr_frac][1]
        assert nbrs.shape[0] == adata.shape[0]
        num_neighbors = nbrs.shape[1] # this will not work for ragged nbr arrays (but we could change it to work)

        nbr_averaged_scores = ( scores + scores[ nbrs ].sum(axis=1) )/(num_neighbors+1)

        plt.scatter(xy[:,0], xy[:,1], c=nbr_averaged_scores, cmap='coolwarm', s=15)
        plt.title('{} ({},{}) {:.1e}'.format(feature, row.gex_cluster, row.tcr_cluster, row.mwu_pvalue_adj))
        plt.xticks([],[])
        plt.yticks([],[])
        if (plotno-1)//ncols == nrows-1:
            plt.xlabel('{} UMAP1'.format(xy_tag.upper()))
        if plotno%ncols == 1:
            plt.ylabel('{} UMAP2'.format(xy_tag.upper()))

    plt.tight_layout()
    plt.savefig(pngfile)
