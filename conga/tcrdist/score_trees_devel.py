## trees from distances -- this was all written ages ago before scipy/sklearn

import string
from os import popen,system,getcwd
import sys
from math import floor,log,exp
#from popen2 import popen2
#from os.path import exists
#from glob import glob

## for making a small glyph that shows the protein-dna interactions
## in a cluster


def IsALeaf(node): return (node[0] == node[1])


def Average_score (leaf_list, leaf_scores, percentile):
    ## negative percentile is signal to use straight average...
    ls = []
    for leaf in leaf_list: ## some of leaf_scores[leaf] could be empty lists...
        ls = ls + leaf_scores[leaf]
    if not ls: ## no score information for this set of leaves
        return None
    elif percentile >= 0:
        assert type(percentile) is int
        ls.sort()
        pos = (percentile * len(ls) ) / 100
        if pos==len(ls): pos = len(ls)-1

        return ls[pos]
    else:
        return sum( ls ) / float( len( ls ) )

class CallAverageScore:
    def __init__( self, percentile ):
        self.percentile = percentile
    def __call__( self, leaf_list, leaf_scores ):
        return Average_score( leaf_list, leaf_scores, self.percentile )



def Make_tree_new(
        distance,
        num_leaves,
        Update_distance_matrix,
        leaf_scores,
        Compute_average_score
):
    # Compute_average_score has to be callable with two args: leaf_list and leaf_scores
    N = num_leaves

    nodes = []
    for i in range(N):
        nodes.append( (i,i,0.0,Compute_average_score([i],leaf_scores) ) )

    for i in range(N): ## initialize distance matrix
        for j in range(N):
            distance[(nodes[i],nodes[j])] = distance[(i,j)]


    while N>1:

        ## find two closest nodes and join them
        min_d = 100000

        for ii, i in enumerate(nodes):
            for jj, j in enumerate(nodes):
                #if i<=j:continue # used to work
                if ii<=jj:continue
                if distance[(i,j)] < min_d:
                    min_d = distance[(i,j)]
                    n1 = i
                    n2 = j

        ##         print "num_nodes: %d   Joining: %s and %s   Distance: %7.3f\n"\
        ##               %(N,Show_small(n1),Show_small(n2),min_d)

        new_node_score = Compute_average_score( Node_members(n1)+Node_members(n2),leaf_scores)

        new_node = (n1,n2,min_d,new_node_score)


        ## update the distances
        Update_distance_matrix (new_node,nodes,distance)

        ## update the node_list
        nodes.append(new_node)
        del nodes[ nodes.index(n1)]
        del nodes[ nodes.index(n2)]

        N = N-1

    return nodes[0]

# def Make_tree_zero_scores( distance, num_leaves, Update_distance_matrix ):
#     N = num_leaves

#     nodes = []
#     for i in range(N):
#         nodes.append( (i,i,0.0,0.0) )

#     for i in range(N): ## initialize distance matrix
#         for j in range(N):
#             distance[(nodes[i],nodes[j])] = distance[(i,j)]


#     while N>1:

#         ## find two closest nodes and join them
#         min_d = 100000

#         for i in nodes:
#             for j in nodes:
#                 if i<=j:continue
#                 if distance[(i,j)] < min_d:
#                     min_d = distance[(i,j)]
#                     n1 = i
#                     n2 = j

#         new_node = (n1,n2,min_d,0.0)

#         ## update the distances
#         Update_distance_matrix (new_node,nodes,distance)

#         ## update the node_list
#         nodes.append(new_node)
#         del nodes[ nodes.index(n1)]
#         del nodes[ nodes.index(n2)]

#         N = N-1
#     return nodes[0]

## the old way
def Make_tree(distance,num_leaves,Update_distance_matrix,leaf_scores,percentile):
    func = CallAverageScore( percentile )
    return Make_tree_new( distance, num_leaves, Update_distance_matrix, leaf_scores, func )

def Copy_tree_update_scores( old_tree, leaf_scores, Compute_average_score ):
    members = Node_members( old_tree )
    score = Compute_average_score( members, leaf_scores )
    if IsALeaf( old_tree ):
        return ( old_tree[0], old_tree[1], old_tree[2], score )
    else:
        return ( Copy_tree_update_scores( old_tree[0], leaf_scores, Compute_average_score ),
                 Copy_tree_update_scores( old_tree[1], leaf_scores, Compute_average_score ),
                 old_tree[2], score )


def Show_tree(tree,names):
    if IsALeaf(tree):
        return names [tree[0]]
    else:
        return '('+Show_tree(tree[0],names)+':'+str(float(tree[2])/2)+','+\
               Show_tree(tree[1],names)+':'+str(float(tree[2])/2)+')'

def Show_small(tree):
    if IsALeaf(tree):
        return repr(tree[0])
    else:
        return '('+Show_small(tree[0])+','+Show_small(tree[1])+')'

def Node_members(node):
    if IsALeaf(node):
        return [node[0]]
    else:
        l1 = Node_members( node[0] )
        l2 = Node_members( node[1] )
        if min(l1)<min(l2):
            return l1+l2
        else:
            return l2+l1


def Update_distance_matrix_AL(new_node,old_nodes,distances): ## average linkage
    n1 = new_node[0]
    n2 = new_node[1]

    l1 = Node_members(new_node)
    distances [ (new_node,new_node)] = 0.0

    for n in old_nodes:
        if n==n1 or n==n2:continue
        l2 = Node_members(n)

        avg = 0.0
        count = 0
        for i in l1:
            for j in l2:
                assert i!=j
                avg = avg+ distances[(i,j)]
                count = count + 1

        distances[(n,new_node)] = avg/count
        distances[(new_node,n)] = avg/count

    return

def Update_distance_matrix_AL_GEOM(new_node,old_nodes,distances):
    dl = list(distances.values())
    dl.sort()
    for i in dl:
        if i!=0:
            min_log = log(i) - 3 ## closer than the closest non-id pair
            break

    n1 = new_node[0]
    n2 = new_node[1]

    l1 = Node_members(new_node)
    distances [ (new_node,new_node)] = 0.0

    for n in old_nodes:
        if n==n1 or n==n2:continue
        l2 = Node_members(n)

        avg = 0.0
        count = 0
        for i in l1:
            for j in l2:
                d = distances[(i,j)]
                assert i!=j
                count = count + 1
                if d == 0.0:
                    avg = avg + min_log
                else:
                    avg = avg + log( d )

        distances[(n,new_node)] = exp( avg / count )
        distances[(new_node,n)] = exp( avg / count )

    return

def Update_distance_matrix_SL(new_node,old_nodes,distances): ## single linkage
    n1 = new_node[0]
    n2 = new_node[1]

    l1 = Node_members(new_node)
    distances [ (new_node,new_node)] = 0.0

    for n in old_nodes:
        if n==n1 or n==n2:continue
        l2 = Node_members(n)

        min_d = 1000
        count = 0
        for i in l1:
            for j in l2:
                assert i!=j
                min_d = min(min_d, distances[(i,j)])

        distances[(n,new_node)] = min_d
        distances[(new_node,n)] = min_d

    return

def Center(tree,node_position,sizes,use_sizes_as_weights=False):
    l = Node_members(tree)
    pos = 0.0
    total_weight=0.0
    for i in l:
        if use_sizes_as_weights:
            pos = pos + sizes[i]*node_position[i]
            total_weight +=sizes[i]
        else:
            pos = pos + node_position[i]
            total_weight+= 1
    pos = pos / total_weight
    return pos


def Size(tree,sizes):
    if IsALeaf(tree):
        return sizes[tree[0]]
    else:
        return Size(tree[0],sizes)+Size(tree[1],sizes)


def Fig_tree(tree,node_position,sizes,use_sizes_as_weights=False): ## edge = [ [x0,y0], [x1,y1], score, size]
    if IsALeaf(tree):
        return []
    else:

        rmsd = tree[2]
        center = Center(tree,node_position,sizes,use_sizes_as_weights)

        c0 = Center(tree[0],node_position,sizes,use_sizes_as_weights)
        r0 = tree[0][2]
        score0 = tree[0][3]
        size0 = Size(tree[0],sizes)
        if IsALeaf(tree[0]):
            cluster0 = tree[0][0]
        else:
            cluster0 = -1
        e0_horizontal = [ [rmsd, c0], [r0,c0], score0, size0, cluster0]
        e0_vertical   = [ [rmsd, c0], [rmsd,center], score0, 1, cluster0]

        c1 = Center(tree[1],node_position,sizes,use_sizes_as_weights)
        r1 = tree[1][2]
        score1 = tree[1][3]
        size1 = Size(tree[1],sizes)
        if IsALeaf(tree[1]):
            cluster1 = tree[1][0]
        else:
            cluster1 = -1
        e1_horizontal = [ [rmsd, c1], [r1, c1], score1, size1 , cluster1]
        e1_vertical   = [ [rmsd, c1], [rmsd,center], score1, 1, cluster1]

        return [ e0_vertical,e0_horizontal,e1_vertical,e1_horizontal] + \
               Fig_tree( tree[0],node_position,sizes,use_sizes_as_weights ) + \
               Fig_tree( tree[1],node_position,sizes,use_sizes_as_weights )

def Node_labels(tree,sizes,node_position,use_sizes_as_weights=False):
    if IsALeaf(tree):return []
    else:
        pos = [tree[2],Center(tree,node_position,sizes,use_sizes_as_weights)]
        size = 0
        for leaf in Node_members(tree):
            size = size+sizes[leaf]
        return [ [ repr(size), pos] ] + \
               Node_labels(tree[0],sizes,node_position,use_sizes_as_weights) + \
               Node_labels(tree[1],sizes,node_position,use_sizes_as_weights)



## return the y-coordinates of the different clusters

def Canvas_tree(
        tree,
        names,
        sizes,
        upper_left,
        lower_right,
        branch_width_fraction,
        plotter,
        label_singletons = False,
        label_internal_nodes = True,
        font=None,
        score_range_for_coloring=None,
        vertical_line_width = 1,
        show_colorful_rmsd_bar = False,
        rmsd_bar_label_stepsize=1,
        rmsd_bar_label_fontsize=10,
        force_min_rmsd=None,
        verbose=False,
):
    ## score_range_for_coloring is a tuple:(mn,mx)
    if score_range_for_coloring: assert len(score_range_for_coloring) == 2

    assert upper_left[0] < lower_right[0]
    assert upper_left[1] < lower_right[1]

    x0,y0 = upper_left
    x1,y1 = lower_right

    #plot_width  = x1-x0
    plot_height = y1-y0

    ## now we are assuming the origin is the top-left corner, like in svg

    ## plot_width and plot_height in pixels

    ## plotter has methods:
    ## .make_line ( [x0,y0], [x1,y1], line_width, normalized_score, extra_tag)
    ## .make_text (text,  [x,y], font)

    branch_width_pixels = plot_height * branch_width_fraction
    #branch_width_pixels = min(100,plot_height/5)

    ## allocate widths for branches; widths measure cluster sizes
    total = sum(sizes)
    w_factor = float( branch_width_pixels) / total
    total = 0
    for s in sizes:
        width = max(1,int(floor(0.5+ s*w_factor))) ## in pixels
        total = total+width
    remainder = plot_height - total
    cluster_width = float(remainder)/len(names) ## padding alotted to each cluster

    if verbose:
        print('branch_width_pixels: {:.2f} plot_height: {:.2f} cluster_padding: {:.3f} w_factor: {:.3f} num_clusters: {} total_members: {}'\
              .format( branch_width_pixels,plot_height,cluster_width,w_factor,len(sizes),sum(sizes)))

    ## position nodes vertically on tree
    nodes = Node_members(tree)
    node_position = {}
    mark = y1
    for i,node in enumerate(nodes):
        width = max(1,int(floor(0.5+ sizes[node] * w_factor)))
        node_position[ node ] = mark-width/2
        mark = mark - cluster_width - width

    edges = Fig_tree(tree,node_position,sizes,use_sizes_as_weights=True) ## each edge = [[x0,y0],[x1,y1],score,size,cluster]

    ## set fontsize: is this still right??

    if not font:
        font = 2*min(25, max (15, int(floor( 0.5 + (cluster_width+7.5)/10))))

    ## rescale the x-positions
    max_rmsd = tree[2]
    min_rmsd = tree[2]
    for e in edges:
        if e[0][0]>0: min_rmsd = min(min_rmsd,e[0][0])
        if e[1][0]>0: min_rmsd = min(min_rmsd,e[1][0])
    min_rmsd = max(0,min_rmsd-0.5)

    if force_min_rmsd is not None:
        min_rmsd = force_min_rmsd

    def Transform(rmsd,min_rmsd = min_rmsd, max_rmsd = max_rmsd):
        return int (floor ( 0.5 + x0 + (x1-x0) * (rmsd - min_rmsd) / (max_rmsd - min_rmsd)))


    ## rescale colors
    scores = []
    for e in edges:
        edge_score = e[2]
        if edge_score != None:
            scores.append(edge_score)

    min_score = min(scores)
    max_score = max(scores)
    if verbose:
        print('min_score:',min_score,'max_score:',max_score,score_range_for_coloring)
    if max_score == min_score:
        max_score = max_score + 1
    if score_range_for_coloring:
        min_score,max_score = score_range_for_coloring

    ## write the edges
    for e in edges:
        start = [ Transform (max(e[0][0],min_rmsd)), e[0][1]] ## rescale x-position
        stop = [ Transform (max(e[1][0],min_rmsd)), e[1][1]]
        edge_score = e[2]
        if edge_score == None:
            normalized_score = None
        else:
            normalized_score = max(0.0, min(1.0, float( e[2] - min_score)/(max_score-min_score) ) )
        line_width = max(1,int(floor(0.5+ e[3]*w_factor)))
        if e[4]>=0: ## it's a real cluster edge
            cluster = e[4]
            extra_tag = 'cluster%02d.%03d'%(cluster,sizes[cluster])
        else:
            extra_tag = 'dummy'

        if start[0] == stop[0]:
            line_width = vertical_line_width

        plotter.make_line(start,stop,line_width,normalized_score,extra_tag)


    ## show scale
    line_x0 = Transform(min_rmsd)
    line_x1 = Transform(max_rmsd)
    if show_colorful_rmsd_bar:
        line_steps = 512
        line_step_size = (line_x1-line_x0)/float(line_steps)
        for i in range( line_steps ):
            norm_score = (float(i)+0.5)/(line_steps)
            #print 'norm_score:',norm_score
            plotter.make_line([ line_x0 + i*line_step_size,y0+5],
                              [ line_x0 +(i+1)*line_step_size,y0+5],3,norm_score )
    else:
        line_ypad = 5
        line_ypad = 2
        plotter.make_line([ line_x0,y0+line_ypad],[ line_x1,y0+line_ypad],2,0.0,color='black' )


    for i in range(int(floor(min_rmsd+1)), 1+int(floor(tree[2])), rmsd_bar_label_stepsize):
        plotter.make_text( str(i), [Transform(i),y0], rmsd_bar_label_fontsize)


    #plotter.make_text( 'Colors: from blue (%7.2f) to red (%7.2f)'%(min_score,max_score),
    #                   [x0,y0+25],25)


    ## label leaves
    for i in range(len(names)):
        if sizes[i] == 1 and not label_singletons: continue
        extra_tag = 'cluster%02d.%03d'%(i,sizes[i])
        plotter.make_text(names[i],
                          [Transform(min_rmsd),node_position[i]],
                          font,extra_tag)

    ## label internal vertices with sizes
    if label_internal_nodes:
        for l in Node_labels (tree,sizes,node_position,use_sizes_as_weights=True):
            plotter.make_text(l[0], [Transform(l[1][0]),l[1][1]], font)

    return ( node_position, Transform, min_rmsd, w_factor ) ## nuisance



