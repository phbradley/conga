from .basic import *
from .all_genes import all_genes, gap_character
from .amino_acids import amino_acids
from .tcr_distances_blosum import blosum, bsd4

## see the TcrDistCalculator at the end for simple tcrdist calculations


##############################################33

GAP_PENALTY_V_REGION = 4
GAP_PENALTY_CDR3_REGION = 12 # same as GAP_PENALTY_V_REGION=4 since WEIGHT_CDR3_REGION=3 is not applied
WEIGHT_V_REGION = 1
WEIGHT_CDR3_REGION = 3
DISTANCE_MATRIX = bsd4
ALIGN_CDR3S = False
TRIM_CDR3S = True


def blosum_character_distance( a, b, gap_penalty ):
    if a== gap_character and b == gap_character:
        return 0
    elif a == '*' and b == '*':
        return 0
    elif a == gap_character or b == gap_character or a=='*' or b=='*':
        return gap_penalty
    else:
        # assert a in amino_acids
        # assert b in amino_acids
        # maxval = min( blosum[(a,a)], blosum[(b,b)] )
        # return maxval - blosum[(a,b)]
        return DISTANCE_MATRIX[ (a,b) ]

def blosum_sequence_distance( aseq, bseq, gap_penalty ):
    assert len(aseq) == len(bseq)
    dist = 0.0
    for a,b in zip(aseq,bseq):
        if a == ' ':
            assert b== ' '
        else:
            dist += blosum_character_distance( a, b, gap_penalty )
    return dist

def align_cdr3_regions( a, b, gap_character ):
    if len(a) == len(b):
        return (a[:],b[:])

    if len(a)<len(b): ## s0 is the shorter sequence
        s0,s1 = a,b
    else:
        s0,s1 = b,a

    lendiff = len(s1)-len(s0)

    best_score=-1000
    best_gappos=0 # in case len(s0) == 1

    # the gap comes after s0[gappos]

    for gappos in range(len(s0)-1):
        score=0
        for i in range(gappos+1):
            score += blosum[ (s0[i],s1[i]) ]
        for i in range(gappos+1,len(s0)):
            score += blosum[ (s0[i],s1[i+lendiff]) ]
        if score>best_score:
            best_score = score
            best_gappos = gappos
    ## insert the gap
    s0 = s0[:best_gappos+1] + gap_character*lendiff + s0[best_gappos+1:]

    assert len(s0) == len(s1)

    if len(a)<len(b): ## s0 is the shorter sequence
        return ( s0, s1 )
    else:
        return ( s1, s0 )

## align
##
##   shortseq[        ntrim: gappos   ] with longseq[ ntrim: gappos ] and
##   shortseq[ -1*remainder: -1*ctrim ] with longseq[ -1*remainder: -1*ctrim ]
##
## but be careful about negative indexing if ctrim is 0
##
## the gap comes after position (gappos-1) ie there are gappos amino acids before the gap
##
##
## DOES NOT INCLUDE THE GAP PENALTY
##
def sequence_distance_with_gappos( shortseq, longseq, gappos ):
    ntrim = 3 if TRIM_CDR3S else 0
    ctrim = 2 if TRIM_CDR3S else 0
    remainder = len(shortseq)-gappos
    dist = 0.0
    count =0
    if ntrim < gappos:
        for i in range(ntrim,gappos):
            dist += DISTANCE_MATRIX[ (shortseq[i], longseq[i] ) ]
            count += 1
    if ctrim < remainder:
        for i in range(ctrim, remainder):
            dist += DISTANCE_MATRIX[ (shortseq[-1-i], longseq[-1-i] ) ]
            count += 1
    return dist, count


def weighted_cdr3_distance( seq1, seq2 ):
    shortseq,longseq = (seq1,seq2) if len(seq1)<=len(seq2) else (seq2,seq1)

    ## try different positions of the gap
    lenshort = len(shortseq)
    lenlong = len(longseq)
    lendiff = lenlong - lenshort

#    assert lenshort>3 ##JCC testing
    assert lenshort > 1##JCC testing
    assert lendiff>=0
    if TRIM_CDR3S:
        # used to say >3+2 but I don't think that's necessary
        assert lenshort >= 3+2 ## something to align... NOTE: Minimum length of cdr3 protein carried into clones file is currently set in the read_sanger_data.py script!

    if not ALIGN_CDR3S:
        ## if we are not aligning, use a fixed gap position relative to the start of the CDR3
        ## that reflects the typically longer and more variable-length contributions to
        ## the CDR3 from the J than from the V. For a normal-length
        ## CDR3 this would be after the Cys+5 position (ie, gappos = 6; align 6 rsds on N-terminal side of CDR3).
        ## Use an earlier gappos if lenshort is less than 11.
        ##
        gappos = min( 6, 3 + (lenshort-5)//2 )
        best_dist, count = sequence_distance_with_gappos( shortseq, longseq, gappos )

    else:
        ## the CYS and the first G of the GXG are 'aligned' in the beta sheet
        ## the alignment seems to continue through roughly CYS+4
        ## ie it's hard to see how we could have an 'insertion' within that region
        ## gappos=1 would be a insertion after CYS
        ## gappos=5 would be a insertion after CYS+4 (5 rsds before the gap)
        ## the full cdr3 ends at the position before the first G
        ## so gappos of len(shortseq)-1 would be gap right before the 'G'
        ## shifting this back by 4 would be analogous to what we do on the other strand, ie len(shortseq)-1-4
        min_gappos = 5
        max_gappos = len(shortseq)-1-4
        while min_gappos>max_gappos:
            min_gappos -= 1
            max_gappos += 1
        for gappos in range( min_gappos, max_gappos+1 ):
            dist, count = sequence_distance_with_gappos( shortseq, longseq, gappos )
            if gappos>min_gappos:
                assert count==best_count
            if gappos == min_gappos or dist < best_dist:
                best_dist = dist
                best_gappos = gappos
                best_count = count
        #print 'align1:',shortseq[:best_gappos] + '-'*lendiff + shortseq[best_gappos:], best_gappos, best_dist
        #print 'align2:',longseq, best_gappos, best_dist


    ## Note that WEIGHT_CDR3_REGION is not applied to the gap penalty
    ##
    return  WEIGHT_CDR3_REGION * best_dist + lendiff * GAP_PENALTY_CDR3_REGION

def compute_all_v_region_distances( organism ):
    rep_dists = {}
    for chain in 'AB': # don't compute inter-chain distances
        repseqs = []
        for id,g in all_genes[organism].items():
            if g.chain == chain and g.region == 'V':
                merged_loopseq = ' '.join( g.cdrs[:-1])
                repseqs.append( ( id, merged_loopseq  ) )
                rep_dists[ id ] = {}

        for r1,s1 in repseqs:
            for r2,s2 in repseqs:
                #if r1[2] != r2[2]: continue
                rep_dists[r1][r2] = WEIGHT_V_REGION * \
                                    blosum_sequence_distance( s1, s2, GAP_PENALTY_V_REGION )

    return rep_dists

def compute_distance(t1, t2, chains, rep_dists):
    ''' t1/2 = [ va_reps, vb_reps, l['cdr3a'], l['cdr3b'] ]
    '''
    dist=0.0
    if 'A' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[0] for y in t2[0] ) ) +\
                weighted_cdr3_distance( t1[2], t2[2] )
    if 'B' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[1] for y in t2[1] ) ) +\
                weighted_cdr3_distance( t1[3], t2[3] )
    return dist

##################################################################################################################
##################################################################################################################

class TcrDistCalculator:
    ''' Simple tcr dist calculator
    precomputes the V-region distances at construction time

    use like this:

    tdist = TcrDistCalculator('human')

    d = tdist(tcr1, tcr2)

    tcr1 and tcr2 are both paired tcrs, represented as tuples of tuples:

    tcr1/2 = ( (va, ja, cdr3a,...), (vb, jb, cdr3b,...) )

    what we really need is:
    * tcr1[0][0] = V-alpha gene
    * tcr1[0][2] = CDR3-alpha
    * tcr1[1][0] = V-beta gene
    * tcr1[1][2] = CDR3-beta

    and same for tcr2
    '''
    def __init__(self, organism):
        self.rep_dists = compute_all_v_region_distances(organism)

    def __call__(self, tcr1, tcr2):

        ''' tcr1 and tcr2 are both paired tcrs, represented as tuples of tuples:

        tcr1/2 = ( (va, ja, cdr3a,...), (vb, jb, cdr3b,...) )

        what we really need is:
        * tcr1[0][0] = V-alpha gene
        * tcr1[0][2] = CDR3-alpha
        * tcr1[1][0] = V-beta gene
        * tcr1[1][2] = CDR3-beta

        and same for tcr2
        '''

        return ( self.rep_dists[tcr1[0][0]][tcr2[0][0]] + weighted_cdr3_distance(tcr1[0][2], tcr2[0][2]) +
                 self.rep_dists[tcr1[1][0]][tcr2[1][0]] + weighted_cdr3_distance(tcr1[1][2], tcr2[1][2]) )


    def single_chain_distance(self, chain1, chain2 ):
        ''' chain1 and chain2 are single-chain tcrs, represented as tuples:

        chain1 = (va, ja, cdr3a,...)  OR  chain1 = (vb, jb, cdr3b,...)

        what we really need is:
        * chain1[0] = V-alpha gene   (or V-beta)
        * chain1[2] = CDR3-alpha     (or CDR3-beta)

        and same for chain2
        '''
        return self.rep_dists[chain1[0]][chain2[0]] + weighted_cdr3_distance(chain1[2], chain2[2])




