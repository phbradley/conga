#from phil import *
from sys import exit
import math
from os.path import exists
from pathlib import Path
import pandas as pd
import numpy as np
from . import util
from . import preprocess as pp
from . import imhc_scoring
from . import cd8_scoring
from .tcrdist.all_genes import all_genes


cdr3_score_FG = 'fg'
cdr3_score_CENTER = 'cen'

cdr3_score_modes = [ cdr3_score_FG, cdr3_score_CENTER ]
default_cdr3_score_mode = cdr3_score_FG

fg_trim = 4
center_len = 5

aa_props_file = Path.joinpath( Path(util.path_to_data), 'aa_props.tsv')
assert exists(aa_props_file)
aa_props_df = pd.read_csv(aa_props_file, sep='\t')
aa_props_df.set_index('aa', inplace=True)

# all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait', 'inkt'] +\
#                      [ '{}_{}'.format(x,y) for x in aa_props_df.columns for y in cdr3_score_modes ]

#tmp hacking SIMPLIFY -- dont include info on which version of the loop is used for scoring
all_tcr_scorenames = ['alphadist', 'cd8', 'cdr3len', 'imhc', 'mait', 'inkt', 'nndists_tcr'] + list(aa_props_df.columns)

amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def read_cd8_score_params():
    # setup the scoring params
    # made by read_flunica_gene_usage_clustermaps*py
    infofile = Path.joinpath( Path(util.path_to_data), 'cd48_score_params_nomait.txt')

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

def old_imhc_score_cdr3( cdr3 ): # cdr3 is untrimmed!
    '''This is the old version of the score, "fit" by staring at the sequences a little bit
    '''
    if len(cdr3) <= 8:
        return 0
    fgloop = cdr3[4:-4]
    return ( len(fgloop) + 3.0 * fgloop.count('C') + 2.0 * fgloop.count('W') + fgloop.count('R') + fgloop.count('K')
             + 0.5*fgloop.count('H') - fgloop.count('D') - fgloop.count('E') )

def old_imhc_score_tcr( tcr ):
    '''This is the old version of the score, "fit" by staring at the sequences a little bit
    '''
    return old_imhc_score_cdr3( tcr[0][2] ) + 2*old_imhc_score_cdr3( tcr[1][2] ) # double-weight the beta CDR3


def cd8_score_tcr( tcr ):
    atcr, btcr = tcr
    return ( cd8_score_tcr_chain( 'A', atcr[0], atcr[1], atcr[2] )[0] +
             cd8_score_tcr_chain( 'B', btcr[0], btcr[1], btcr[2] )[0] )

def is_human_mait_alpha_chain(atcr):
    return ( atcr[0].startswith('TRAV1-2') and
             ( atcr[1].startswith('TRAJ33') or
               atcr[1].startswith('TRAJ20') or
               atcr[1].startswith('TRAJ12') ) and
             len(atcr[2]) == 12 )

def is_mouse_mait_alpha_chain(atcr):
    return ( atcr[0].startswith('TRAV1*') and atcr[1].startswith('TRAJ33') and len(atcr[2]) == 12 )

def is_mouse_inkt_alpha_chain(atcr):
    return ( atcr[0].startswith('TRAV11') and atcr[1].startswith('TRAJ18') and len(atcr[2]) == 15 )

def is_human_inkt_tcr(tcr):
    # could also put some limits on cdr3s?
    return ( tcr[0][0].startswith('TRAV10') and
             tcr[0][1].startswith('TRAJ18') and
             len(tcr[0][2]) in [14,15,16] and  # 15 seems to be the consensus
             tcr[1][0].startswith('TRBV25') )


def mait_score_tcr(tcr, organism):
    if 'human' in organism:
        return float(is_human_mait_alpha_chain(tcr[0]))

    elif 'mouse' in organism:
        return float(is_mouse_mait_alpha_chain(tcr[0]))

    else:
        print('unrecognized organism:', organism)
        exit()
        return 0.


def inkt_score_tcr(tcr, organism):
    if 'human' in organism:
        return float(is_human_inkt_tcr(tcr))

    elif 'mouse' in organism:
        return float(is_mouse_inkt_alpha_chain(tcr[0]))

    else:
        print('unrecognized organism:', organism)
        exit()
        return 0.


def read_locus_order( remove_slashes_from_gene_names= False ):
    ''' returns  all_locus_order:  all_locus_order[ab][gene] = int(l[1])
    '''
    # read the gene order from imgt
    all_locus_order = {'A':{}, 'B':{}}
    for ab in 'AB':
        fn = f'imgt_tr{ab.lower()}_locus_order.txt'
        filename = Path.joinpath( Path( util.path_to_data ), fn)
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
        #print('alphadist_score_tcr: unrecognized va:', va)
        va_dist = 0.5*(len(trav_list)-1)

    if ja in traj_list:
        ja_dist = traj_list.index(ja)
    else:
        #print('alphadist_score_tcr: unrecognized ja:', ja)
        ja_dist = 0.5*(len(traj_list)-1)
    return va_dist + ja_dist

def cdr3len_score_tcr(tcr):
    ''' double-weight the beta chain
    '''
    return len(tcr[0][2]) + 2*len(tcr[1][2])



def property_score_cdr3(cdr3, score_name, score_mode):
    ''' Currently averages the score over the number of aas scores
    '''
    global aa_props_df
    global cdr3_score_CENTER
    global cdr3_score_FG
    global fg_trim
    global center_len

    col = aa_props_df[score_name]

    if score_mode == cdr3_score_CENTER:
        cdr3 = cdr3[1:] # trim off the first 'C', makes structurally symmetric
        if len(cdr3)<center_len:
            return np.mean(col)
        else:
            ntrim = (len(cdr3)-center_len)//2
            return sum(col[aa] for aa in cdr3[ntrim:ntrim+center_len])/center_len

    elif score_mode == cdr3_score_FG:
        fgloop = cdr3[fg_trim:-fg_trim]
        if fgloop:
            return sum( col[aa] for aa in fgloop )/len(fgloop)
        else:
            return np.mean(col) # was returning 0 here; prob should use aa-frequency-weighted average...
    else:
        print( 'property_score_cdr3:: unrecognized score_mode:', score_mode)
        exit()
        return 0.0

def property_score_tcr(tcr, score_name, score_mode, alpha_weight=1.0, beta_weight=1.0):
    return ( alpha_weight * property_score_cdr3(tcr[0][2], score_name, score_mode ) +
             beta_weight  * property_score_cdr3(tcr[1][2], score_name, score_mode ) )


def make_tcr_score_table(adata, scorenames):
    ''' Returns an array of the tcr scores of shape: (adata.shape[0], len(scorenames))
    '''
    global aa_props_df
    organism = adata.uns['organism']

    tcrs = pp.retrieve_tcrs_from_adata(adata)
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    organism_genes = all_genes[organism]
    genes = frozenset( organism_genes.keys())
    count_reps = frozenset( [ x.count_rep for x in organism_genes.values() ] )

    cols = []
    for name in scorenames:
        if name == 'cdr3len':
            cols.append( [ cdr3len_score_tcr(x) for x in tcrs ])
        elif name.startswith('tcr_cluster'):
            num = int(name[11:])
            cols.append( [ float(x==num) for x in clusters_tcr])
        elif name == 'alphadist':
            cols.append( [ alphadist_score_tcr(x) for x in tcrs ])
        elif name == 'oldcd8':# the 'old' cd8 score
            cols.append( [ cd8_score_tcr(x) for x in tcrs ])
        elif name == 'cd8': # see comparison between old/new in cd8_scoring.py
            cols.append(cd8_scoring.make_cd8_score_table_column(tcrs))
        elif name == 'old_imhc':
            cols.append( [ old_imhc_score_tcr(x) for x in tcrs ])
        elif name == 'imhc':
            cols.append( imhc_scoring.make_imhc_score_table_column(tcrs, aa_props_df))
        elif name == 'mait':
            organism = adata.uns['organism']
            cols.append( [ mait_score_tcr(x, organism) for x in tcrs ])
        elif name == 'inkt':
            organism = adata.uns['organism']
            cols.append( [ inkt_score_tcr(x, organism) for x in tcrs ])
        elif name == 'nndists_tcr':
            if 'nndists_tcr' not in adata.obs_keys():
                print('WARNING nndists_tcr score requested but not present in adata.obs!!!!')
                cols.append( np.zeros( adata.shape[0] ) )
            else:
                cols.append( np.array(adata.obs['nndists_tcr']) )
        elif name == 'N_ins': # number of N insertions
            if 'N_ins' not in adata.obs_keys():
                print('WARNING N_ins requested but not present in adata.obs!!!!')
                cols.append( np.zeros( adata.shape[0] ) )
            else:
                cols.append( np.array(adata.obs['N_ins']).astype(float) )
        elif name in genes:
            matched = False
            for i_ab,ab in enumerate('AB'):
                for i_vj,vj in enumerate('VJ'):
                    ii_genes = set([x for x,y in organism_genes.items() if y.chain==ab and y.region==vj ])
                    if name in ii_genes:
                        assert not matched
                        matched = True
                        cols.append( [float(x[i_ab][i_vj]==name) for x in tcrs])
            assert matched

        elif name in count_reps:
            matched = False
            for i_ab,ab in enumerate('AB'):
                for i_vj,vj in enumerate('VJ'):
                    ii_count_reps = set([x.count_rep for x in organism_genes.values() if x.chain==ab and x.region==vj ])
                    if name in ii_count_reps:
                        assert not matched
                        matched = True
                        cols.append( [float(organism_genes[x[i_ab][i_vj]].count_rep==name) for x in tcrs])
            assert matched

        else:
            score_mode = name.split('_')[-1]
            score_name = '_'.join( name.split('_')[:-1])
            if score_mode not in cdr3_score_modes:
                score_mode = default_cdr3_score_mode
                score_name = name
            cols.append( [ property_score_tcr(x, score_name, score_mode) for x in tcrs ] )
    table = np.array(cols).transpose()#[:,np.newaxis]
    #print( table.shape, (adata.shape[0], len(scorenames)) )
    assert table.shape == (adata.shape[0], len(scorenames))

    return table

