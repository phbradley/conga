import numpy as np
import pandas as pd

from . import basic
from .basic import *
from .amino_acids import amino_acids
from . import all_genes
from . import tcrdist_svg_basic # hacky code duplication
from . import tcr_distances
from . import util
from . import tcr_sampler ## for analyze_junction

junction_bars = True

font_family = MONOSPACE_FONT_FAMILY # in basic, from ../convert_svg_to_png.py

greek_alpha = 'A'
greek_beta  = 'B'

# these don't always seem to work
#greek_alpha = '&#x3b1;'
#greek_beta  = '&#x3b2;'

junction_bars_color = { 'V':  'silver',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'black',
                        'J':  'dimgray' }

gap_character = '-' ## different from some other places

xmargin = 10
ymargin = 10

default_vj_logo_width = 200
default_pwmplusgaps_width = 770
default_xpad = 30

default_pwm_height = 100
default_junction_bars_height = 35.0 if junction_bars else 0.0
default_ypad = 3

default_ab_glyphs_spacer = 15

#used for rescaling below
default_width = 2*default_vj_logo_width + 2*default_xpad + default_pwmplusgaps_width
default_height = default_pwm_height + default_ypad + default_junction_bars_height


boxpad = 2 #OK to keep this fixed independent of scaling??

## returns
def make_tcr_logo(
        upper_left, # [x0,y0] location of the upper left corner
        tcrs, # list of tuples with tcr info in specific order, see docstring below
        members, # list of integers, indices into tcrs
        all_dists, # matrix of distances between tcrs
        ab,# A or B
        organism,
        rep_colors, # dictionary of colors indexed by the V and J gene names in the tcrs list
        vj_logo_width,
        pwmplusgaps_width,
        xpad,
        pwm_height,
        junction_bars_height,
        ypad,
        show_full_cdr3
):
    """ Returns a list of commands for svg output via tcrdist_svg_basic.create_file

    width = vj_logo_width + xpad + pwmplusgaps_width + xpad + vj_logo_width
    height = pwm_height + ypad + junction_bars_height

    tcrs = [( va_rep, ja_rep, vb_rep, jb_rep, cdr3a, cdr3b, cdr3a_nucseq_src, cdr3b_nucseq_src ), ...]


    """
    cmds = []

    def trim_gene_name_for_logo(gene, vj, ab=ab, organism=organism):
        if vj=='J':
            return gene[4:]
        if ( '_gd' in organism and ab=='B' or # could be TRAV or TRDV
             '_ig' in organism and ab=='A' ): # could be IGKV or IGLV
            return gene[2]+gene[4:]
        return gene[4:]


    assert ab in ['A','B']
    assert len(all_dists) == len(tcrs) and len(all_dists[0]) == len(tcrs)

    if ab== 'A':
        v_index, j_index, cdr3_index, nucseq_src_index = 0, 1, 4, 6
        junction_bars_order = ['V','N','J']
        greek_letter = greek_alpha
    else:
        assert ab=='B'
        v_index, j_index, cdr3_index, nucseq_src_index = 2, 3, 5, 7
        junction_bars_order = ['V','N1','D','N2','J']
        greek_letter = greek_beta

    real_size = len(members)

    distl = []
    for m1 in members:
        avgdis = 0.
        for m2 in members:
            avgdis += all_dists[m1][m2]
        avgdis/=len(members)
        distl.append( ( avgdis, m1 ) )
    distl.sort()

    center = distl[0][1]
    #print('center_tcr:', len(tcrs), distl[0], tcrs[center][v_index], tcrs[center][cdr3_index] )

    ## count v,j gene reps
    v_count = {}
    j_count = {}
    for rep in [tcrs[x][v_index] for x in members ]: v_count[rep] = v_count.get(rep,0)+1
    for rep in [tcrs[x][j_index] for x in members ]: j_count[rep] = j_count.get(rep,0)+1

    center_cdr3 = tcrs[center][ cdr3_index ]
    if not show_full_cdr3: center_cdr3 = center_cdr3[3:-2]
    L = len(center_cdr3)

    pwm = {}
    junction_pwm = {}
    gap_count = {}
    for i in range(L):
        pwm[i] = dict(list(zip(amino_acids+[gap_character],[0]*21)))
        gap_count[i]=0
    for i in range(3*L):
        junction_pwm[i] = dict( list(zip( junction_bars_order+[gap_character],
                                     [0.]*(1+len(junction_bars_order)))))

    for member in members:
        member_cdr3 = tcrs[member][cdr3_index]
        member_junction = tcrs[member][nucseq_src_index] ## a list
        if not show_full_cdr3:
            member_cdr3 = member_cdr3[3:-2]
            member_junction = member_junction[9:-6]
        assert len(member_junction) == 3*len(member_cdr3)
        a,b = tcr_distances.align_cdr3_regions( center_cdr3, member_cdr3, gap_character )

        for i in range(len(a)):
            if a[i] == gap_character:
                if i and a[i-1]==gap_character:continue
                gap_count[i-1] += 1

            else: ## b[i] could be gap_character or an amino_acid
                pwmpos = i - a[:i].count(gap_character)
                pwm[pwmpos][ b[i] ] += 1
                for j in range(3*pwmpos,3*pwmpos+3):
                    if b[i] == gap_character:
                        junction_pwm[j]['-'] += 1
                    else:
                        bpos = i-b[:i].count(gap_character)
                        assert b[i] == member_cdr3[bpos] ## sanity
                        bsrc = member_junction[3*bpos + j-3*pwmpos ]
                        junction_pwm[j][bsrc] += 1

    ## normalize the pwms
    for i in range(L):
        tot = float(sum(pwm[i].values()))
        for aa in pwm[i]:
            pwm[i][aa] /= tot
    for i in range(3*L):
        tot = float(sum(junction_pwm[i].values()))
        for aa in junction_pwm[i]:
            junction_pwm[i][aa] /= tot

    ## now we want to make a pwm
    total_gaps = sum(gap_count.values())


    column_width = pwmplusgaps_width / ( L + float(total_gaps)/len(members) )

    y0 = upper_left[1]
    y1 = y0 + pwm_height

    ## make a v-gene logo
    vl = [(y, trim_gene_name_for_logo(x, 'V'), rep_colors[x]) for x,y in v_count.items()]
    jl = [(y, trim_gene_name_for_logo(x, 'J'), rep_colors[x]) for x,y in j_count.items()]

    #single_glyph_width = 2* vj_logo_width + pwmplusgaps_width + 2*xpad
    #x0 = xmargin + glyph_size_text_width + ii_ab2 * ( single_glyph_width + ab_glyphs_spacer )
    x0 = upper_left[0]

    cmds.append( tcrdist_svg_basic.make_stack( (x0,y0), (x0 + vj_logo_width,y1), vl ) )

    ## a box around the V-logo
    cmds.append( tcrdist_svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + vj_logo_width+boxpad,y1+boxpad ),
                                      'none', 'black', stroke_width=1 ) )

    if junction_bars: ## label V-logo down below
        text = 'V'+greek_letter
        fontsize = junction_bars_height*0.9
        p0 = [ x0 + 0.5* vj_logo_width - 0.6*fontsize, y1+junction_bars_height ]
        cmds.append( tcrdist_svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

    x0 += vj_logo_width + xpad

    cmds.append( tcrdist_svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + pwmplusgaps_width+boxpad,y1+boxpad ),
                                      'none', 'black', stroke_width=1 ) )



    #print ic,L,size,y0,y1

    ## now show each column
    prev_gap_column_width = 0.0
    jb_rights = [] #debugging
    for pos in range(L):
        ## first the column of aas
        colpwm={}
        colpwm[0] = pwm[pos]
        #if verbose:
        #    print 'colpwm:',pos,pwm[pos]
        cmds.append( tcrdist_svg_basic.protein_logo( (x0,y0), (x0+column_width,y1), colpwm ) )

        save_x0 = x0 ## for junction_bars

        x0 += column_width

        ## any gaps?
        if gap_count[pos]:
            gap_column_width = float( column_width * gap_count[pos] ) / len(members)
            cmds.append( tcrdist_svg_basic.text_in_box( (x0,y0), (x0+gap_column_width,y1), gap_character, 'black' ) )

            x0+= gap_column_width
        else:
            gap_column_width = 0.0


        if junction_bars:
            junction_bar_width = ( column_width + gap_column_width/2. + prev_gap_column_width/2. )/3.
            junction_bar_x0 = save_x0 - prev_gap_column_width/2.
            # print 'left:',junction_bar_x0,'right:',junction_bar_x0 + 3.*junction_bar_width,\
            #     'prev_gap_column_width:',prev_gap_column_width,'gap_column_width:',gap_column_width,\
            #     'save_x0:',save_x0,'column_width:',column_width,'junction_bar_width:',junction_bar_width
            if jb_rights:
                assert abs( junction_bar_x0 - jb_rights[-1] )<1e-3
            jb_rights.append( junction_bar_x0 + 3.*junction_bar_width )

            y2 = y1+junction_bars_height
            for j in range(3):
                col = junction_pwm[3*pos+j]
                lcol = [ ( col[x],x) for x in junction_bars_order ]
                # lcol = [ (y,x) for x,y in col.iteritems()]
                # lcol.sort()
                # lcol.reverse()
                y1shift = y1+ ypad
                ## largest at the top
                for frac,a in lcol:
                    if a==gap_character: continue
                    y1shift_next = y1shift + frac * junction_bars_height
                    color = junction_bars_color[ a ]
                    p0 = [ junction_bar_x0+ j   *junction_bar_width, y1shift]
                    p1 = [ junction_bar_x0+(j+1)*junction_bar_width, y1shift_next ]
                    cmds.append( tcrdist_svg_basic.rectangle( p0, p1, fill=color, stroke=color ) )
                    y1shift = y1shift_next


        prev_gap_column_width = gap_column_width


    ## now the J-logo
    cmds.append( tcrdist_svg_basic.make_stack( (x0,y0), (x0+vj_logo_width,y1), jl ) )

    cmds.append( tcrdist_svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + vj_logo_width+boxpad,y1+boxpad ),
                                      'none', 'black', stroke_width=1 ) )

    if junction_bars: ## label V-logo down below
        text = 'J'+greek_letter
        fontsize = junction_bars_height * 0.9
        p0 = [ x0 + 0.5* vj_logo_width - 0.6*fontsize, y1+junction_bars_height ]
        cmds.append( tcrdist_svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )


    return cmds #####################################





def make_default_logo_svg_cmds(
        upper_left, width, height, organism, tcr_infos, chain,
        tcrdist_calculator = None,
        add_fake_alleles = False,
        show_full_cdr3 = False
):
    # right now single-chain only
    # returns cmds

    assert chain in 'AB'

    if tcrdist_calculator is None:
        tcrdist_calculator = tcr_distances.TcrDistCalculator(organism)

    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( tcr_infos, organism )

    rep_colors = {}
    for info in tcr_infos:
        for vj in 'vj':
            for abl in 'ab':
                rep   = info[vj+abl+'_label_rep']
                color =  info[vj+abl+'_label_rep_color']
                rep_colors[rep] = color

    tcrs = []

    dist_tcrs = []

    def add_fake_allele_info(x):
        if '*' not in x:
            return x+'*01'
        else:
            return x

    for l in tcr_infos:

        mouse = None #l['subject']
        epitope = None #l['epitope']
        cdr3a = l['cdr3a']
        cdr3b = l['cdr3b']

        ## for computing distances
        va_gene = l['va_gene']
        ja_gene = l['ja_gene']
        va_genes = l['va_genes'].split(';')

        vb_gene = l['vb_gene']
        jb_gene = l['jb_gene']
        vb_genes = l['vb_genes'].split(';')


        # add '*01' -- hacky!
        if add_fake_alleles:
            va_genes = list(map( add_fake_allele_info, va_genes ))
            vb_genes = list(map( add_fake_allele_info, vb_genes ))
            va_gene = add_fake_allele_info( va_gene )
            ja_gene = add_fake_allele_info( ja_gene )
            vb_gene = add_fake_allele_info( vb_gene )
            jb_gene = add_fake_allele_info( jb_gene )

        va_reps = frozenset( ( all_genes.all_genes[organism][x].rep for x in va_genes ))
        vb_reps = frozenset( ( all_genes.all_genes[organism][x].rep for x in vb_genes ))

        if chain == 'A':
            dist_tcrs.append( (va_gene, ja_gene, cdr3a) )
        else:
            dist_tcrs.append( (vb_gene, jb_gene, cdr3b) )

        #dist_tcrs.append( [ va_reps, vb_reps, cdr3a, cdr3b ] )
        #all_info.append( l )


        ## note that we are using mm1 reps here that also dont have allele info
        va_rep = l['va_label_rep']
        ja_rep = l['ja_label_rep']
        vb_rep = l['vb_label_rep']
        jb_rep = l['jb_label_rep']

        cdr3a_nucseq_src = ['V']*(3*len(cdr3a)) ## hack, unused
        cdr3b_nucseq_src = ['V']*(3*len(cdr3b))
        if junction_bars:

            if chain == 'A':
                a_junction_results = tcr_sampler.analyze_junction( organism, va_gene, ja_gene,
                                                                   cdr3a, l['cdr3a_nucseq'].lower(),
                                                                   return_cdr3_nucseq_src=True )
                cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,\
                    a_trims, a_inserts, cdr3a_nucseq_src = a_junction_results
            elif chain == 'B':
                b_junction_results = tcr_sampler.analyze_junction( organism, vb_gene, jb_gene,
                                                                   cdr3b, l['cdr3b_nucseq'].lower(),
                                                                   return_cdr3_nucseq_src=True )

                cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,\
                    b_trims, b_inserts, cdr3b_nucseq_src = b_junction_results
                ## try to distinguish between N before D and N after D
                for i in range(len(cdr3b_nucseq_src)):
                    if cdr3b_nucseq_src[i] == 'N':
                        if cdr3b_nucseq_src[:i].count('D')==0:
                            cdr3b_nucseq_src[i] = 'N1'
                        else:
                            cdr3b_nucseq_src[i] = 'N2'


        assert len(cdr3a_nucseq_src) == 3*len(cdr3a)
        assert len(cdr3b_nucseq_src) == 3*len(cdr3b)
        #print cdr3b, cdr3b_nucseq_src
        tcrs.append( ( va_rep, ja_rep, vb_rep, jb_rep, cdr3a, cdr3b,
                       cdr3a_nucseq_src, cdr3b_nucseq_src ) )


    ## compute distances, used in logo construction for picking the center tcr for aligning against
    #print 'computing distances:',len(dist_tcrs)
    all_dists = np.zeros( ( len(dist_tcrs), len(dist_tcrs)) )
    for i,t1 in enumerate( dist_tcrs ):
        for j in range(i+1,len(dist_tcrs)):
            dist = tcrdist_calculator.single_chain_distance(t1, dist_tcrs[j])
            all_dists[i][j] = dist
            all_dists[j][i] = dist

    # now make the logo
    members = list(range(len(tcrs)))


    scale_w = float( width ) / default_width
    scale_h = float( height ) / default_height

    # scale everything by our desired height, width

    return make_tcr_logo( upper_left, tcrs, members, all_dists, chain, organism, rep_colors,
                          scale_w * default_vj_logo_width,
                          scale_w * default_pwmplusgaps_width,
                          scale_w * default_xpad,
                          scale_h * default_pwm_height,
                          scale_h * default_junction_bars_height,
                          scale_h * default_ypad,
                          show_full_cdr3 )




def make_tcr_logo_svg_commands_for_tcrs(
        tcrs,
        chain,
        organism,
        upper_left,
        width_height,
        tcrdist_calculator=None # so we can avoid recalculating V gene distances over and over again
):
    ''' tcrs is a list of tuples: tcrs = [ (atcr1,btcr1), (atcr2,btcr2), ....
    atcr1 = (va1, ja1, cdr3a1, cdr3a_nucseq1, *)
    btcr1 = (vb1, jb1, cdr3b1, cdr3b_nucseq1, *)

    upper_left and width_height are tuples (x,y)

    returns the svg cmds

    untested, currently under development
    '''
    assert chain in 'AB'

    add_fake_alleles = False
    show_full_cdr3 = False

    # SILLY backwards compatibility, fill out info that would have been in clones file
    infos = []
    for tcr in tcrs:
        atcr, btcr = tcr
        info = {}
        for iab, ab in enumerate('ab'):
            info[f'cdr3{ab}'] = tcr[iab][2]
            info[f'cdr3{ab}_nucseq'] = tcr[iab][3]
            for ivj, vj in enumerate('vj'):
                gene = tcr[iab][ivj]
                assert gene in all_genes.all_genes[organism]
                for tag in 'gene genes rep reps'.split():
                    info[f'{vj}{ab}_{tag}'] = gene
                info[f'{vj}{ab}_countreps'] = all_genes.all_genes[organism][gene].count_rep
        infos.append(info)

    if tcrdist_calculator is None:
        tcrdist_calculator = tcr_distances.TcrDistCalculator(organism)

    cmds = make_default_logo_svg_cmds(
        upper_left, width_height[0], width_height[1], organism, infos, chain,
        tcrdist_calculator = tcrdist_calculator,
        add_fake_alleles = add_fake_alleles, show_full_cdr3 = show_full_cdr3 )


    return cmds


def make_tcr_logo_for_tcrs(
        tcrs,
        chain,
        organism,
        pngfile,
        tcrdist_calculator=None # so we can avoid recalculating V gene distances over and over again
):
    ''' tcrs is a list of tuples: tcrs = [ (atcr1,btcr1), ....
    atcr1 = (va1, ja1, cdr3a1,*)
    btcr1 = (vb1, jb1, cdr3b1,*)
    '''
    assert chain in 'AB'
    assert pngfile.endswith('.png')

    add_fake_alleles = False
    show_full_cdr3 = False

    # backwards compatibility, fill out info that would have been in clones file
    infos = []
    for tcr in tcrs:
        atcr, btcr = tcr
        info = {}
        for iab, ab in enumerate('ab'):
            info[f'cdr3{ab}'] = tcr[iab][2]
            info[f'cdr3{ab}_nucseq'] = tcr[iab][3]
            for ivj, vj in enumerate('vj'):
                gene = tcr[iab][ivj]
                assert gene in all_genes.all_genes[organism]
                for tag in 'gene genes rep reps'.split():
                    info[f'{vj}{ab}_{tag}'] = gene
                info[f'{vj}{ab}_countreps'] = all_genes.all_genes[organism][gene].count_rep
        infos.append(info)

    if tcrdist_calculator is None:
        tcrdist_calculator = tcr_distances.TcrDistCalculator(organism)

    cmds = make_default_logo_svg_cmds(
        [xmargin,ymargin], default_width, default_height, organism, infos,
        chain,tcrdist_calculator = tcrdist_calculator,
        add_fake_alleles = add_fake_alleles, show_full_cdr3 = show_full_cdr3 )


    svg_width = 2*xmargin + default_width
    svg_height = 2*ymargin + default_height

    tcrdist_svg_basic.create_file( cmds, svg_width, svg_height, pngfile[:-4]+'.svg', create_png=True )

    #print default_width,default_height


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if __name__ == '__main__':

    # with Parser(locals()) as p:
    #     #p.str('args').unspecified_default().multiple().required()
    #     p.str('clones_file').required()
    #     p.str('organism').required()
    #     p.str('outfile_prefix')
    #     p.flag('add_fake_alleles')       # --flag_arg  (no argument passed)
    #     p.flag('show_full_cdr3')       # --flag_arg  (no argument passed)
    #     #p.flag('verbose')       # --flag_arg  (no argument passed)
    #     p.multiword('ABs').cast(lambda x: x.split())
    #     p.multiword('epitopes').cast(lambda x: x.split())
    clones_file = 'tmpclones.tsv'
    organism = 'human'
    add_fake_alleles = False
    show_full_cdr3 = False
    AB = 'A'
    outfile_prefix = None

    junction_bars = True

    if outfile_prefix is None:
        outfile_prefix = clones_file[:-4]

    df = pd.read_csv( clones_file, sep='\t' )
    infos = [ x._asdict() for x in df.itertuples() ]

    for l in infos:
        for vj in 'vj':
            for ab in 'ab':
                countreps_tag = '{}{}_countreps'.format(vj, ab)
                if countreps_tag not in l:
                    genes_tag = '{}{}_genes'.format(vj, ab)
                    genes = l[genes_tag].split(';')
                    l[countreps_tag] = ';'.join(util.countreps_from_genes(genes, organism))

    tcrdist_calculator = tcr_distances.TcrDistCalculator(organism)

    cmds = make_default_logo_svg_cmds( [xmargin,ymargin], default_width, default_height, organism, infos, AB,
                                       tcrdist_calculator = tcrdist_calculator,
                                       add_fake_alleles = add_fake_alleles, show_full_cdr3 = show_full_cdr3 )


    svg_width = 2*xmargin + default_width
    svg_height = 2*ymargin + default_height

    tcrdist_svg_basic.create_file( cmds, svg_width, svg_height, outfile_prefix+'.svg', create_png=True )

    #print default_width,default_height
