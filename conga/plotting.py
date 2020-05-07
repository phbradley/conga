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
import os
from os.path import exists
from collections import Counter
import matplotlib.image as mpimg
from matplotlib.collections import LineCollection
from scipy.cluster import hierarchy
from scipy.spatial import distance

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
        bbox_to_anchor=(1, 0.5),
        ncol=(1 if len(categories) <= 20 # was 14
              else 2 if len(categories) <= 30 else 3),
        fontsize=legend_fontsize)
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
                                num_genes_to_show = 10 ):
    ''' ranks is a list of (gene,l2r,pval)
    '''

    def pval_factor( pval, min_pval=top_pval_for_max_height ):
        return math.sqrt( -1 * math.log10( max(min_pval,pval) ) )
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
        cmds.append( svg_basic.text_in_box( [x0,y0], [x1,y1], gene, color='black' ) )

        y_offset+= height
    return cmds


def make_single_rank_genes_logo( ranks, svgfile,
                                 margin = 10,
                                 logo_width = 300,
                                 logo_max_height = 400,
                                 top_pval_for_max_height = 1e-50,
                                 min_pval_for_scaling = 1e-300,
                                 num_genes_to_show = 10,
                                 create_png = True,
                                 create_html = False ):
    cmds = []

    upper_left=  [margin, margin]
    cmds.extend( make_rank_genes_logo_stack( ranks, upper_left, logo_width, logo_max_height,
                                             top_pval_for_max_height=top_pval_for_max_height,
                                             min_pval_for_scaling=min_pval_for_scaling,
                                             num_genes_to_show=num_genes_to_show ) )

    svg_basic.create_file( cmds, logo_width+2*margin, logo_max_height + 2*margin, svgfile,
                           create_png=create_png, create_html=create_html, background_color='white' )

def _parse_clusterpair_rank_genes( adata, uns_tag = 'rank_genes_good_cluster_pairs' ):
    ''' returns all_ranks dict mapping from cluspair to list [ (gene, l2r, pval), ... ]
    '''

    rgtag = 'rank_genes_good_cluster_pairs'

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



def make_tcr_logo( tcrs, organism, ab, pngfile):
    clonesfile = pngfile+'_clones.tsv'
    util.make_clones_file( tcrs, clonesfile )

    exe = '{}make_tcr_logo.py'.format( util.TCRDIST_REPO )
    if not exists(exe):
        print( 'ERROR: unable to locate python script in tcr-dist repo:', exe)
        exit()

    #> /dev/null 2> /dev/null'\
    cmd = '{} {} --organism {} --ABs {} --outfile_prefix {} --clones_file {}'\
        .format( util.PYTHON2_EXE, exe, organism, ab, pngfile[:-4], clonesfile )
    util.run_command(cmd)
    if not exists(pngfile):
        print('make_tcr_logo:: cmd failed:', cmd)
        exit()
    os.remove(clonesfile)



def make_logo_plots(
        adata,
        nbrs_gex,
        nbrs_tcr,
        min_cluster_size,
        logo_pngfile,
        clusters_gex= None,
        clusters_tcr= None,
        ignore_tcr_cluster_colors = False,
        gex_cluster_names = None,
        rank_genes_uns_tag = 'rank_genes_good_cluster_pairs'
):
    ''' need:
    * gex/tcr clusters: (obsm)
    * clones files for cluster-pairs with enough "good" clones
    * clone sizes (obs: clone_size)
    * tcrs (passed in)
    * 2d projections: gex and tcr (obsm: X_gex_2d, X_tcr_2d)
    * conga scores (aka all_max_nbr_chisq): (obsm: conga_scores)
    * X_igex (obsm: X_igex) and X_igex_genes (uns)
    * rank genes info for each cluster-pair
    '''

    ## unpack data from adata arrays ##################################
    clone_sizes = adata.obs['clone_sizes']
    if clusters_gex is None:
        clusters_gex = adata.obs['clusters_gex']
    if clusters_tcr is None:
        clusters_tcr = adata.obs['clusters_tcr']
    good_score_mask = np.array(list(adata.obs['good_score_mask']))
    X_gex_2d = adata.obsm['X_gex_2d']
    X_tcr_2d = adata.obsm['X_tcr_2d']
    conga_scores = adata.obsm['conga_scores'] # (num_clones,3)
    X_igex = adata.obsm['X_igex']
    X_igex_genes = list(adata.uns['X_igex_genes']) #for .index
    organism = adata.uns['organism']

    tcrs = pp.retrieve_tcrs_from_adata(adata)
    tcr_scoretags = ['mhci_score', 'cdr3len_score', 'cd8_score', 'alphadist_score']
    all_tcr_scores = np.array( [np.array(adata.obs[x]) for x in tcr_scoretags] ).transpose()
    assert all_tcr_scores.shape == ( adata.shape[0], len(tcr_scoretags))
    ################################################################# no more unpacking below here...

    num_clones = adata.shape[0]
    num_good_clones = np.sum(good_score_mask)
    assert X_igex.shape == (num_clones, len(X_igex_genes))


    logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL',
                  'GNLY','PRF1','GZMA','GZMB','GZMK','GZMH',
                  'CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ]

    header2_genes = ['CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1']
    #header2_genes = ['CD4','CD8A','CD8B','SELL','GNLY','GZMA','GZMB','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'] #MME
    num_header2_tcr_score_cols = 2

    gene_width = 6


    assert len(logo_genes) == 3*gene_width - 2

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
        gex_nbrhood_tcr_scores.append( np.sum( all_tcr_scores[nbrhood_mask,:], axis=0 )/num )

    gex_nbrhood_X_igex = np.array(gex_nbrhood_X_igex)
    assert gex_nbrhood_X_igex.shape == X_igex.shape
    gex_nbrhood_tcr_scores = np.array(gex_nbrhood_tcr_scores)
    assert gex_nbrhood_tcr_scores.shape == all_tcr_scores.shape


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


    if 0:
        all_nbrs = {} # sets for faster membership checking
        for ii in node2cluspair:
            all_nbrs[ii] = [ frozenset(x for x in nbrs_gex[ii] if good_score_mask[x]),
                             frozenset(x for x in nbrs_tcr[ii] if good_score_mask[x]) ]

        all_sym_nbrs = {} # symmetrized versions of nbrs
        for ii, ii_nbrs in all_nbrs.items():
            new_gex_nbrs = frozenset( x for x in ii_nbrs[0] if ii in all_nbrs[x][0] )
            new_tcr_nbrs = frozenset( x for x in ii_nbrs[1] if ii in all_nbrs[x][1] )
            all_sym_nbrs[ii] = [ new_gex_nbrs, new_tcr_nbrs ]

        all_nbrs = all_sym_nbrs # use the symmetrized versions


    # setup dbl nbrs
    all_dbl_nbrs = {}
    for ii in node2cluspair:
        gex_nbrs = frozenset( nbrs_gex[ii])
        all_dbl_nbrs[ii] = [ x for x in nbrs_tcr[ii] if good_score_mask[x] and x in gex_nbrs ]


    ## parse the X_igex matrix
    all_scores = {} # map from clp to [ num_clusters_tcr, fracs, means]

    for g in logo_genes:
        if g not in X_igex_genes:
            print('X_igex is missing logo_gene:', g, X_igex_genes)
    logo_gene_indices = [ X_igex_genes.index(x) if x in X_igex_genes else None for x in logo_genes ]

    all_scores = {}
    for clp, nodes in cluspair2nodes.items():
        if len(nodes) < min_cluster_size:
            continue
        X = X_igex[nodes,:]
        means = np.mean(X,axis=0)
        fracs = np.sum(X>1e-3, axis=0) / float(len(nodes))

        fracs = [ 0 if x is None else fracs[x] for x in logo_gene_indices ]
        means = [ 0 if x is None else means[x] for x in logo_gene_indices ]
        all_scores[clp] = [ len(nodes), fracs, means ]



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
    tcr_logo_width = 8
    gex_logo_width = 2
    logo_height = 0.5
    frac_scale = 130
    margin = 0.2

    gex_logo_width = logo_height * ( gene_width / (math.sqrt(3)+1) )

    fig_width = 2*margin + dendro_width + title_logo_width + rg_logo_width + tcr_logo_width + gex_logo_width
    header_height = (fig_width-2*margin)/5.
    header2_height = (fig_width-2*margin)/(len(header2_genes) + num_header2_tcr_score_cols)
    fig_height = 2*margin + header_height + 2*header2_height + len(clps)*logo_height

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

    # make a header
    for icol,all_xy in enumerate( [all_xy_gex, all_xy_tcr] ):
        proj_tag = 'GEX TCR'.split()[icol]

        bottom = (fig_height-margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+2*icol*header_height)/fig_width # make these plots square, so height=width
        width = header_height/fig_width
        plt.axes( [left,bottom,width,height] )

        ## make a plot colored by pval
        ks = sorted( all_xy.keys() )
        xy = np.array( [ all_xy[x] for x in ks ] )
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
        plt.title('kNN overlap E-values ({}2D)'.format(proj_tag),fontsize=11,pad=1)
        #plt.title('kNN overlap E-values {}-UMAP'.format('GEX' if icol==0 else 'TCR'),fontsize=12,pad=1)

        plt.text(0.01,0.01,'{} clones'.format(len(all_xy)), ha='left', va='bottom', fontsize=8,
                 transform=plt.gca().transAxes)

        bottom = (fig_height-margin-header_height)/fig_height
        height = header_height/fig_height
        left = (margin+(2*icol+1)*header_height)/fig_width # make these plots square, so height=width
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
        plt.title('interesting clones ({}2D)'.format(proj_tag),fontsize=11,pad=1)


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


    # show the clusters on the two different UMAP projections
    for icol,(xy, clusters) in enumerate( zip( [X_gex_2d, X_tcr_2d], [clusters_gex, clusters_tcr] ) ):
        if icol==1 and ignore_tcr_cluster_colors:
            clusters = clusters_gex
        bottom = (fig_height-margin-0.5*header_height)/fig_height
        height = 0.5*header_height/fig_height
        left = (margin+4*header_height+0.5*icol*header_height)/fig_width # make these plots square, so height=width
        width = 0.5*header_height/fig_width
        plt.axes( [left,bottom,width,height] )
        proj_tag = 'GEX TCR'.split()[icol]
        num_clusters = np.max(clusters)+1
        C = get_integers_color_dict( num_clusters)
        colors = [ C[x] for x in clusters]
        plt.scatter( xy[:,0], xy[:,1], s=small_markersize, c=colors )
        plt.xticks([],[])
        plt.yticks([],[])
        clusters_tag = proj_tag if not ignore_tcr_cluster_colors else 'combo'
        plt.title('{} clusters ({}2D)'.format(clusters_tag, proj_tag),fontsize=7,pad=8)
        xmn,xmx = plt.xlim()
        ymn,ymx = plt.ylim()
        # show the colors
        for i in range(num_clusters):
            plt.text((i+1)/(num_clusters+1), 1.005, str(i), color=C[i], ha='center', va='bottom', fontsize=6,
                     transform=plt.gca().transAxes)
            # plt.scatter( [(i+0.5)/(num_clusters+1)], [1.0], s=small_markersize, color=C[i],
            #              transform=plt.gca().transAxes)
        plt.xlim((xmn,xmx))
        plt.ylim((ymn,ymx))



    ## make a key for the gex logo
    bottom = (fig_height-margin-header_height)/fig_height
    height = 0.5*header_height/fig_height
    left = (margin+4*header_height)/fig_width
    width = header_height/fig_width
    plt.axes( [left,bottom,width,height] )
    yshift = math.sqrt(3)/2 # 30 60 90 triangle
    xy = []
    for r in range(3):
        xoff = 0 if r==1 else 0.5
        yoff = yshift * ( 1-r )
        ng = gene_width if r==1 else gene_width-1
        for ii in range(ng):
            xy.append( [xoff+ii, yoff] )

    for gene,(x,y) in zip( logo_genes, xy ):
        plt.text(x,y,gene,fontsize=11,ha='center',va='center',rotation=30)
    plt.xlim( [-0.5, gene_width-0.5] )
    plt.ylim( [-yshift-0.5, yshift+0.5 ] )
    plt.xticks([],[])
    plt.yticks([],[])
    plt.title('key for GEX logo',fontsize=8,va='top',pad=1)

    # make a line of GEX plots, right now for the logo_genes
    for ig, gene in enumerate(header2_genes):
        if gene not in X_igex_genes:
            print('missing header2_gene:', gene,'from X_igex_genes')
            continue
        index = X_igex_genes.index(gene)
        bottom = (fig_height-margin-header_height-header2_height)/fig_height
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
            plt.ylabel('GEX UMAP2')
        plt.text( xmn+pad, ymx, 'raw', va='top', ha='left', fontsize=8)
        plt.text( xmx, ymx, gene, va='top', ha='right', fontsize=8)

        ## make another plot with the gex-nbrhood averaged scores
        bottom = (fig_height-margin-header_height-2*header2_height)/fig_height
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
        plt.scatter( X_gex_2d[:,0], X_gex_2d[:,1], s=small_markersize, c=colors,cmap='coolwarm', vmin=-mx, vmax=mx )
        plt.xticks([],[])
        plt.yticks([],[])
        xmn,xmx = plt.xlim()
        _,ymx = plt.ylim()
        pad = 0.02*(xmx-xmn)
        if ig==0:
            plt.ylabel('GEX UMAP2')
        plt.text( xmn+pad, ymx, 'nbrZ', va='top', ha='left', fontsize=8)
        plt.text( xmx, ymx, gene, va='top', ha='right', fontsize=8)

    for ii, scoretag in enumerate(tcr_scoretags):
        irow=ii//2
        icol=ii%2
        ## make another plot with the gex-nbrhood averaged scores
        bottom = (fig_height-margin-header_height-(irow+1)*header2_height)/fig_height
        height = header2_height/fig_height
        left = (fig_width-margin-(icol+1)*header2_height)/fig_width # make these plots square, so height=width
        width = header2_height/fig_width
        plt.axes( [left,bottom,width,height] )
        vals = all_tcr_scores[:,ii]
        mean, std = np.mean(vals), np.std(vals)
        if std<0.01:
            print('low_std:', scoretag, std)
            std = 0.01
        colors = (gex_nbrhood_tcr_scores[:,ii]-mean)/std

        mx = max( 0.2, np.max(np.abs(colors)))
        plt.scatter( X_gex_2d[:,0], X_gex_2d[:,1], s=small_markersize, c=colors,cmap='coolwarm', vmin=-mx, vmax=mx )
        plt.xticks([],[])
        plt.yticks([],[])
        xmn,xmx = plt.xlim()
        _,ymx = plt.ylim()
        pad = 0.02*(xmx-xmn)
        plt.text( xmx-pad, ymx, scoretag[:-6].replace('alphadist','ad').replace('cdr3len','len'),
                  va='top', ha='right', fontsize=8)
        plt.text( xmn+pad, ymx, 'tcr', va='top', ha='left', fontsize=8)

        # make nbr-averaged plot

    if len(clps)>1:
        # make a dendrogram along the LHS
        bottom = margin/fig_height
        height = (fig_height-2*margin-header_height-2*header2_height)/fig_height
        left = (margin)/fig_width
        width = dendro_width/fig_width
        ax = plt.axes( [left,bottom,width,height] )
        tree_method = 'average'
        tree_optimal_ordering = False

        Z = hierarchy.linkage( distance.squareform(D,force='tovector'),
                               method=tree_method, optimal_ordering=tree_optimal_ordering )

        hierarchy.dendrogram( Z, orientation='left', ax=ax, link_color_func= lambda x:'k' )
        plt.xlim((1.03,0.0))
        plt.axis('off')
        plt.text(0.0,0.0, 'Clusters (size>{:d})'.format(min_cluster_size),
                 ha='left', va='top', transform=plt.gca().transAxes)
        leaves = list( hierarchy.leaves_list( Z ) )
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
        bottom = (fig_height-margin-header_height-2*header2_height-(irow+1)*logo_height)/fig_height
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
        if gex_cluster_names is None:
            plt.text(-0.1,0,str(clp[0]),va='center',ha='right')
        else:
            plt.text(-0.1,0,gex_cluster_names[clp[0]],va='center',ha='right',fontsize=6)
        if not ignore_tcr_cluster_colors:
            plt.text( 0.1,0,str(clp[1]),va='center',ha='left')
        plt.text(-1,1,'{}'.format(num_nodes),va='top',ha='left',fontsize=8)
        plt.text( 1,1,'{}'.format(num_cells),va='top',ha='right',fontsize=8)
        # if agbt:
        #     # figure out the top two pmhcs and their umis
        #     pmhc_umis = Counter()
        #     for node in nodes:
        #         for pmhc,umi in all_pmhc_umis[node].items():
        #             pmhc_umis[pmhc] += umi
        #     topl = pmhc_umis.most_common(2)
        #     plt.text(-1,0,'{}\n{}\n{:.2f}'.format(topl[0][0][:3], topl[0][0][4:], topl[0][1]/len(nodes)),
        #              va='center',ha='left',fontsize=6)
        #     plt.text( 1,0,'{}\n{}\n{:.2f}'.format(topl[1][0][:3], topl[1][0][4:], topl[1][1]/len(nodes)),
        #               va='center',ha='right',fontsize=6)

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
            make_tcr_logo( [ tcrs[x] for x in cluspair2nodes[clp] ], organism, ab, pngfile )
            image = mpimg.imread(pngfile)
            plt.imshow(image)
            plt.axis('off')
            if last_row:
                plt.text(0.5,-0.05, 'TCR{} logo'.format(ab), ha='center', va='top', transform=plt.gca().transAxes)

        # make the gex logo
        left = (margin+dendro_width+title_logo_width+rg_logo_width+tcr_logo_width)/fig_width
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

        plt.scatter( xy[:,0], xy[:,1], s=sizes, c=colors, cmap = plt.get_cmap('Reds'), vmin=0, vmax=3.5 )
        plt.xlim( [-0.5, gene_width-0.5] )
        plt.ylim( [-yshift-0.5, yshift+0.5 ] )
        plt.xticks([],[])
        plt.yticks([],[])
        if last_row:
            plt.text(0.5,-0.05, 'GEX logo', ha='center', va='top', transform=plt.gca().transAxes)
    plt.savefig(logo_pngfile, dpi=300)


    for tmpfile in tmpfiles:
        if exists(tmpfile):
            os.remove(tmpfile)
        svgfile = tmpfile[:-3]+'svg'
        if exists(svgfile):
            os.remove(svgfile)

