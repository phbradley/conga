import numpy as np
#import pandas as pd

#from . import basic
from .basic import *
from .amino_acids import amino_acids
#from . import all_genes
from . import tcrdist_svg_basic # hacky code duplication
from . import score_trees_devel
from . import tcr_distances
from . import make_tcr_logo
#from . import util
#from . import tcr_sampler ## for analyze_junction


def make_tcr_tree_svg_commands(
        tcrs,
        organism,
        upper_left,
        width_height,
        tcrdist_matrix, # symmetric numpy matrix
        tcrdist_clustering_threshold=77.1, # try to emulate classic tcr-dist paired threshold (?)
        tcrdist_calculator=None,
        max_tcrs_for_trees=300, # otherwise the actual svg tree gets too tall/dense
        color_scores=None,
        color_score_range=None,
        title=''
):
    ''' tcrs is a list of tuples: tcrs = [ (atcr1,btcr1), (atcr2,btcr2), ....
    atcr1 = (va1, ja1, cdr3a1, cdr3a_nucseq1, *)
    btcr1 = (vb1, jb1, cdr3b1, cdr3b_nucseq1, *)

    upper_left and width_height are tuples (x,y)

    returns the svg cmds

    untested, currently under development
    '''

    assert tcrdist_matrix.shape == (len(tcrs), len(tcrs))
    if color_scores is None:
        color_scores = [1.0]*len(tcrs)
    if color_score_range is None:
        color_score_range = [np.min(color_scores), np.max(color_scores)+1e-6] # add epsilon

    if tcrdist_calculator is None:
        tcrdist_calculator = tcr_distances.TcrDistCalculator(organism)

    font_family = MONOSPACE_FONT_FAMILY # in basic, from ../convert_svg_to_png.py
    gap_character = '-' ## different from some other places
    min_cluster_size = 1
    min_cluster_size_for_glyphs = 5
    min_cluster_fraction_for_glyphs = 0.03

    max_covered_fraction_for_glyphs = 0.75

    max_glyphs = 10

    xmargin = 0.01 * width_height[0]
    ymargin = 0.01 * width_height[1]
    glyphs_to_tree_spacer = 0.01 * width_height[0]
    glyph_size_text_width = 0.01 * width_height[0]

    glyph_region_width = width_height[0] * 1200.0/2200
    ab_glyphs_spacer = 0.03 * glyph_region_width
    tree_width = width_height[0] - glyph_region_width-glyph_size_text_width-glyphs_to_tree_spacer-2*xmargin
    glyph_height = width_height[1] * 135.0/2000
    tree_height = width_height[1] - 2*ymargin

    tree_x0 = upper_left[0] + xmargin + glyph_size_text_width + glyph_region_width + glyphs_to_tree_spacer
    tree_y0 = upper_left[1] + ymargin

    branch_width_fraction = 0.3

    log10_of_zero = -100

    N = len(tcrs)

    all_nbrs = []
    all_dists = tcrdist_matrix #tmphack
    radius = tcrdist_clustering_threshold

    # note that all_nbrs[i] includes i
    for i in range(N):
        all_nbrs.append(np.nonzero(all_dists[i,:] <= radius)[0] )

    deleted = np.full((N,), False)

    centers = []
    all_members = []

    while True:
        clusterno = len(centers)

        nbr_counts = np.array([ np.sum((~deleted)[x]) for x in all_nbrs ]) # all_nbrs is ragged: list of arrays
        nbr_counts[deleted] = 0
        best_nbr_count = np.max(nbr_counts)
        center = np.argmax(nbr_counts)

        if best_nbr_count < min_cluster_size:
            break

        centers.append( center )

        members = [center]
        deleted[center] = True
        for nbr in all_nbrs[center]:
            if not deleted[nbr]:
                deleted[nbr] = True
                members.append( nbr )

        assert len(members) == best_nbr_count
        all_members.append( frozenset(members) )

    ## possibly subsample
    tree_indices = list(range(len(tcrs)))
    if len(tree_indices) > max_tcrs_for_trees:
        tree_indices = random.sample( tree_indices, max_tcrs_for_trees )

    ## now get some info together for plotting
    all_center_dists = {}
    all_scores = []
    sizes = []
    names = []

    real_cluster_number2fake_cluster_number = {}
    fake_cluster_number2real_cluster_number = {}
    fake_ic=-1
    real_sizes = []
    for ic,center in enumerate(centers):
        ok_members = []
        for m in all_members[ic]:
            if m in tree_indices:
                ok_members.append(m)
        if not ok_members:
            real_cluster_number2fake_cluster_number[ic] = -1
            continue
        fake_ic += 1
        real_cluster_number2fake_cluster_number[ic] = fake_ic
        fake_cluster_number2real_cluster_number[fake_ic] = ic
        assert fake_ic == len(sizes)
        size = len(ok_members)
        real_size = len(all_members[ic])
        all_scores.append( [ color_scores[x] for x in ok_members ] )
        sizes.append( size )
        real_sizes.append( real_size )
        names.append('')

    ## now do cluster center distances now that we know which ones are in the plot
    for ic,center in enumerate(centers):
        fake_ic = real_cluster_number2fake_cluster_number[ic]
        if fake_ic<0: continue
        for jc,other_center in enumerate(centers):
            fake_jc = real_cluster_number2fake_cluster_number[jc]
            if fake_jc<0: continue
            all_center_dists[ (fake_ic,fake_jc) ] = all_dists[center][other_center]


    print('num_tcrs:',len(tcrs),'num_clusters:',len(centers),'fake_num_tcrs',sum(sizes),\
        'fake_num_clusters:',len(sizes))

    cmds = []
    toplabel_fontsize = 15
    toplabel_xpos=upper_left[0]+xmargin
    toplabel_ypos=upper_left[1]+ymargin+0.9*toplabel_fontsize # x,y passed to make_text is the lower_left corner
    if title:
        cmds.append( tcrdist_svg_basic.make_text(title, (toplabel_xpos, toplabel_ypos),
                                                 toplabel_fontsize, font_family=font_family))
        toplabel_ypos += toplabel_fontsize
    cmds.append( tcrdist_svg_basic.make_text('{} xcrs {} clusters'.format(len(tcrs), len(centers)),
                                             (toplabel_xpos, toplabel_ypos),
                                             toplabel_fontsize, font_family=font_family))
    toplabel_ypos += toplabel_fontsize
    cmds.append( tcrdist_svg_basic.make_text('color_range: {:.2f}-{:.2f}'\
                                             .format(color_score_range[0], color_score_range[1]),
                                             (toplabel_xpos, toplabel_ypos),
                                             toplabel_fontsize, font_family=font_family))


    percentile = -1 # I believe this means average the scores when coloring a branch of the tree

    tree = score_trees_devel.Make_tree( all_center_dists, len(names),
                                        score_trees_devel.Update_distance_matrix_AL,
                                        all_scores, percentile )


    plotter = tcrdist_svg_basic.SVG_tree_plotter()

    tree_p0 = [tree_x0, tree_y0]
    tree_p1 = [tree_x0 + tree_width, tree_y0 + tree_height ]


    ## node_position tells us where the different clusters are located, vertically
    ##
    force_min_rmsd = min(74.0, max(0.0, tcrdist_clustering_threshold-3.0))
    node_position, Transform, canvas_tree_min_rmsd, canvas_tree_w_factor = score_trees_devel.Canvas_tree(
        tree, names, sizes, tree_p0, tree_p1, branch_width_fraction,
        plotter, label_internal_nodes = False,
        score_range_for_coloring = color_score_range,
        rmsd_bar_label_stepsize=75, force_min_rmsd=force_min_rmsd )

    max_rmsd_for_glyphs = 3.0*radius

    my_min_cluster_size_for_glyphs = max( min_cluster_size_for_glyphs,
                                          int(floor(0.5 + min_cluster_fraction_for_glyphs*len(tcrs))))

    while True:


        ## get all internal horizontal edges (ie subtrees) that are merged at an rmsd below a threshold
        ## and have at least my_min_cluster_size_for_glyphs
        def get_good_edges( subtree, node_position, sizes, real_sizes ):
            if score_trees_devel.IsALeaf( subtree ):
                return []
            else:
                big_rmsd = subtree[2]
                edges = []
                for ii in range(2):
                    iitree = subtree[ii]
                    little_rmsd = iitree[2]
                    center= score_trees_devel.Center( iitree, node_position, sizes, use_sizes_as_weights=True )
                    real_size = score_trees_devel.Size( iitree, real_sizes )
                    if little_rmsd <= max_rmsd_for_glyphs and real_size >= my_min_cluster_size_for_glyphs:
                        ## this tree is OK
                        edges.append( tuple( ( little_rmsd, big_rmsd, center, real_size,
                                               tuple( sorted( score_trees_devel.Node_members(iitree) ) ) ) ) )
                    edges.extend( get_good_edges( iitree, node_position, sizes, real_sizes ) )
            return edges


        good_edges = get_good_edges( tree, node_position, sizes, real_sizes )

        # if True:
        #     for edge in good_edges:
        #         print 'good:',epitope, ab, edge
            #exit()

        ## figure out which edges we should draw glyphs for
        ## take all the true clusters
        ## try to maximize coverage?
        ##
        glyph_cmds = []
        glyph_edges = []
        while len(glyph_edges) < max_glyphs and len(glyph_edges)<len(good_edges):

            ## figure out who is covered
            covered = set()
            for e in glyph_edges:
                for fake_ic in e[4]:
                    covered.add( fake_ic )

            new_edge = ()
            best_score = -1e6
            for e in good_edges:
                if e in glyph_edges: continue
                (little_rmsd,big_rmsd,center,real_size,clusters) = e
                new_size = real_size
                for fake_ic in clusters:
                    if fake_ic in covered:
                        new_size -= real_sizes[fake_ic]
                if new_size<my_min_cluster_size_for_glyphs: continue
                if (real_size-new_size) > max_covered_fraction_for_glyphs * real_size: continue

                if little_rmsd<1e-3: ## true cluster
                    sortscore = 1000 + new_size
                else:
                    sortscore = new_size - ( little_rmsd / max_rmsd_for_glyphs ) * len(tcrs)
                if sortscore > best_score:
                    new_edge = e[:]
                    best_score = sortscore
            if not new_edge: break
            #print 'new_edge:',new_edge
            glyph_edges.append( new_edge )

            ## let's make a box around this edge in the tree
            (little_rmsd,big_rmsd,center,real_size,clusters) = new_edge
            fake_size = sum( [sizes[x] for x in clusters] )
            line_width = max(1,int(floor(0.5+ fake_size*canvas_tree_w_factor )))
            box_x0 = Transform(max(canvas_tree_min_rmsd,little_rmsd)) ; box_x1 = Transform(big_rmsd)
            assert little_rmsd<big_rmsd
            glyph_cmds.append( tcrdist_svg_basic.rectangle( ( box_x0, center - line_width/2.0 ),
                                                            ( box_x1, center + line_width/2.0 ),
                                                            'none', 'black', stroke_width=3, dashed=True ) )
            glyph_cmds.append( tcrdist_svg_basic.make_text( '%d'%real_size, ( box_x1+4, center + 0.75*15.0/2) ,
                                                            15, font_family=font_family))

        if len( glyph_edges ) < max_glyphs and my_min_cluster_size_for_glyphs>1:
            my_min_cluster_size_for_glyphs -= 1
        else:
            break

    cmds.extend(glyph_cmds)

    min_glyph_loc = toplabel_ypos + glyph_height/2
    max_glyph_loc = tree_p1[1] - glyph_height/2

    #glyph_location = dict( [ (e,max(min_glyph_loc,min(max_glyph_loc,e[2]))) for e in glyph_edges ] )
    glyph_location = dict( [ (e,e[2]) for e in glyph_edges ] )



    while True: ## keep looping until we are bump-free
        l = [(y,x) for x,y in glyph_location.items() ]

        l.sort()
        bump = False
        stepsize = 5 ## pixels
        for (loc1,edge1),(loc2,edge2) in zip( l[:-1], l[1:] ):
            sep = loc2-loc1
            if sep < 1.25*glyph_height:
                bump = True
                if glyph_location[edge1]-stepsize > min_glyph_loc:
                    glyph_location[edge1] -= stepsize
                if glyph_location[edge2]+stepsize < max_glyph_loc:
                    glyph_location[edge2] += stepsize
                #print('bump1', loc1, loc2)
        if bump:
            continue

        # first try to move UP if necessary (ie lower loc)
        # move the highest one first
        # these are sorted from smallest to largest location, ie from highest to lowest since y-location increase down
        for (loc,edge) in l:
            if glyph_location[edge] > max_glyph_loc:
                glyph_location[edge] -= stepsize
                #print('bump2', loc)
                bump = True
        if bump:
            continue

        for (loc,edge) in reversed(l): # from bottom to top
            if glyph_location[edge] < min_glyph_loc:
                glyph_location[edge] += stepsize
                print('bump3', loc)
                bump = True

        if not bump:
            break

    #exit()

    ## now let's draw some cluster summary glyphs next to the cluster tree


    ## write out the glyph-cluster sizes
    for glyph_edge in glyph_edges:
        real_size = glyph_edge[3]
        loc = glyph_location[glyph_edge]
        ## silly xloc, was hard-coded to 5
        cmds.append( tcrdist_svg_basic.make_text('%3d'%real_size, ( upper_left[0]+xmargin-5, loc+glyph_height/2.) , 15,
                                                 font_family=font_family))



    for ii_ab2, ab2 in enumerate('AB'):
        single_glyph_width = float(glyph_region_width-ab_glyphs_spacer) / 2 # A or B

        for glyph_edge in glyph_edges:
            (little_rmsd, big_rmsd, y_center, real_size, clusters ) = glyph_edge

            members = []
            for fake_ic in clusters:
                members.extend( all_members[ fake_cluster_number2real_cluster_number[ fake_ic ] ] )

            assert len(members) == real_size

            y0 = glyph_location[glyph_edge] - glyph_height/2.
            x0 = upper_left[0] + xmargin + glyph_size_text_width + ii_ab2 * ( single_glyph_width + ab_glyphs_spacer )

            glyph_tcrs = [tcrs[x] for x in members]
            cmds.extend(make_tcr_logo.make_tcr_logo_svg_commands_for_tcrs(
                glyph_tcrs, ab2, organism, [x0,y0], [single_glyph_width, glyph_height], tcrdist_calculator))


    return plotter.cmds + cmds




