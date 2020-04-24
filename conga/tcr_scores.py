#from phil import *
import math
from os.path import exists

def read_cd8_score_params():
    # setup the scoring params
    # made by read_flunica_gene_usage_clustermaps*py
    infofile = '/home/pbradley/gitrepos/conga/conga/data/cd48_score_params_nomait.txt'

    scoretags = 'cdr3_len cdr3_aa gene'.split()

    all_scorevals = {'A':{}, 'B':{} }
    for tag in scoretags:
        all_scorevals['A'][tag] = {}
        all_scorevals['B'][tag] = {}


    min_l2e = -0.6
    max_l2e = 0.6

    # scores reflect cd8 versus cd4 preference

    for line in open(infofile,'r'):
        l = line.split()
        tag = l[0]
        assert tag.endswith('_cdx_pval:')
        tag = tag[:-10]
        key = l[10]
        if tag[:4] == 'cdr3':
            ab = tag[4]
            tag = tag[:4]+tag[5:]
        else:
            assert key[:2] == 'TR'
            ab = key[2]
        assert tag in scoretags
        assert l[6] == 'cd4:'
        cd4_mean = float(l[7])
        cd8_mean = float(l[9])
        tot = 0.5*(cd4_mean + cd8_mean)
        if not tot:
            l2e = 0.0
        elif cd8_mean==0:
            l2e = min_l2e
        else:
            l2e = max( min_l2e, min( max_l2e, math.log( cd8_mean / tot, 2.0 ) ) )

        all_scorevals[ab][tag][key] = l2e
        #print '{:7.3f} {} {} {}'.format(l2e,ab,tag,key)
    return all_scorevals

all_scorevals = read_cd8_score_params()

def cd8_score_tcr_chain( ab, v, j, cdr3 ):
    global all_scorevals

    cdr3_len_ranges = { 'A': ( min( int(x) for x in all_scorevals['A']['cdr3_len'] ),
                               max( int(x) for x in all_scorevals['A']['cdr3_len'] ) ),
                        'B': ( min( int(x) for x in all_scorevals['B']['cdr3_len'] ),
                               max( int(x) for x in all_scorevals['B']['cdr3_len'] ) ) }

    assert ab in 'AB'
    if '*' in v:
        v = v[:v.index('*')]
    if '*' in j:
        j = j[:j.index('*')]
    v_score = all_scorevals[ab]['gene'].get(v,0.0)
    j_score = all_scorevals[ab]['gene'].get(j,0.0)
    mn,mx = cdr3_len_ranges[ab]
    L = max( mn, min( mx, len(cdr3) ) )
    cdr3_len_score = all_scorevals[ab]['cdr3_len'][str(L)]
    cdr3_aa_score = 0.0
    if L>8:
        for aa in cdr3[4:-4]:
            cdr3_aa_score += all_scorevals[ab]['cdr3_aa'][aa]
    score = v_score + j_score + cdr3_len_score + cdr3_aa_score

    return ( score, v_score, j_score, cdr3_len_score, cdr3_aa_score )

def mhci_score_cdr3( cdr3 ): # cdr3 is untrimmed!
    if len(cdr3) <= 8:
        return 0
    fgloop = cdr3[4:-4]
    return ( len(fgloop) + 3.0 * fgloop.count('C') + 2.0 * fgloop.count('W') + fgloop.count('R') + fgloop.count('K')
             + 0.5*fgloop.count('H') - fgloop.count('D') - fgloop.count('E') )

def mhci_score_tcr( tcr ):
    return mhci_score_cdr3( tcr[0][2] ) + 2*mhci_score_cdr3( tcr[1][2] ) # double-weight the beta CDR3 score


def cd8_score_tcr( tcr ):
    atcr, btcr = tcr
    return ( cd8_score_tcr_chain( 'A', atcr[0], atcr[1], atcr[2] )[0] +
             cd8_score_tcr_chain( 'B', btcr[0], btcr[1], btcr[2] )[0] )


def read_locus_order( remove_slashes_from_gene_names= False ):
    ''' returns  all_locus_order:  all_locus_order[ab][gene] = int(l[1])
    '''
    # read the gene order from imgt
    all_locus_order = {'A':{}, 'B':{}}
    for ab in 'AB':
        filename = '/home/pbradley/gitrepos/conga/conga/data/imgt_tr{}_locus_order.txt'.format(ab.lower())
        assert exists(filename)
        for line in open(filename,'r'):
            l = line.split()
            if l[0] == 'IMGT':continue
            assert len(l) in [2,3]
            gene = l[0]
            if '/' in gene and remove_slashes_from_gene_names:
                print( 'remove /',gene)
                gene = gene.replace('/','')
            all_locus_order[ab][gene] = int(l[1])
    return all_locus_order

all_locus_order = read_locus_order()

trav_list = [ x[1] for x in sorted( (y,x) for x,y in all_locus_order['A'].items() if x[:4] == 'TRAV' ) ]
traj_list = [ x[1] for x in sorted( (y,x) for x,y in all_locus_order['A'].items() if x[:4] == 'TRAJ' ) ]


def alphadist_score_tcr( tcr ):
    global trav_list
    global traj_list
    va, ja = tcr[0][:2]
    va = va[:va.index('*')]
    ja = ja[:ja.index('*')]
    if va in trav_list:
        va_dist = len(trav_list)-1 -trav_list.index(va)
    else:
        print('alphadist_score_tcr: unrecognized va:', va)
        va_dist = 0.5*(len(trav_list)-1)

    if ja in traj_list:
        ja_dist = traj_list.index(ja)
    else:
        print('alphadist_score_tcr: unrecognized ja:', ja)
        ja_dist = 0.5*(len(traj_list)-1)
    return va_dist + ja_dist

def cdr3len_score_tcr(tcr):
    ''' double-weight the beta chain
    '''
    return len(tcr[0][2]) + 2*len(tcr[1][2])
