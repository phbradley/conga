#from . import amino_acids
from .amino_acids import amino_acids
import sys

def get_alphabet( pwm ):
    alphabet = list(pwm[0].keys())[:]
    alphabet.sort()
    return alphabet

def check_pwm( pwm, tol= 0.001 ):
    L = len(pwm)
    alphabet = get_alphabet( pwm )
    for pos in range(L):
        for aa,val in pwm[pos].items(): assert val > -1e-6
        total = sum( ( pwm[pos][aa] for aa in alphabet ) )
        assert abs( total - 1.0 ) - tol

def create_protein_pwm_from_sequences( seqs, pseudocounts = 0.0 ):
    return create_pwm_from_sequences( seqs, amino_acids, pseudocounts=pseudocounts )

def create_dna_pwm_from_sequences( seqs, pseudocounts=0.0 ):
    return create_pwm_from_sequences( seqs, list('acgt'), pseudocounts=pseudocounts )

def create_pwm_from_sequences( seqs, alphabet, pseudocounts=0.0 ):
    pwm = {}
    if len(seqs) == 0: return pwm
    L = len( seqs[0] )

    for pos in range(L):
        pwm[ pos ] = dict( list(zip( alphabet, [pseudocounts]*len(alphabet) )) )

    for s in seqs:
        assert len(s) == L
        for pos in range(L):
            if s[pos] not in alphabet:
                sys.stderr.write('logo_tools.create_pwm_from_sequences: skipping bad character %s\n'%(s[pos]))
                continue
            pwm[ pos ][ s[pos] ] += 1

    for pos in range(L):
        norm = 1.0 / sum( pwm[pos].values() )
        for a in alphabet: pwm[ pos ][ a ] *= norm

    check_pwm( pwm )

    return pwm



base_partner = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',
                'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                'R': 'Y', 'Y': 'R',
                'S': 'S', 'W': 'W',
                'K': 'M', 'M': 'K',
                '.': '.' }

nucleotide_classes_lower_case = { 'a':'a',
                                  'c':'c',
                                  'g':'g',
                                  't':'t',
                                  'w':'at',
                                  's':'cg',
                                  'k':'gt',
                                  'm':'ac',
                                  'y':'ct',
                                  'r':'ag',
                                  'b':'cgt',
                                  'd':'agt',
                                  'h':'act',
                                  'v':'acg',
                                  'n':'acgt' }

nuc_match_lower_case = {}

for x in nucleotide_classes_lower_case:
    for y in nucleotide_classes_lower_case:
        nuc_match_lower_case[(x,y)] = False
        for x1 in nucleotide_classes_lower_case[x]:
            for y1 in nucleotide_classes_lower_case[y]:
                if x1==y1:
                    nuc_match_lower_case[(x,y)] = True
                    break

def nucleotide_symbols_match( a_in, b_in ):
    a = a_in.lower()
    b = b_in.lower()
    if a==b: return True
    return nuc_match_lower_case.get( (a,b), False )

def degnucs_contain(a,b):
    ''' Check if a contains b, allowing for degenerate nucleotide symbols
    '''
    #
    if len(b)>len(a):
        return False
    for offset in range(len(a)-len(b)):
        if all(nucleotide_symbols_match(x,y)
               for x,y in zip(a[offset:], b)):
            return True
    return False

def degnucs_index(a,b):
    ''' Where does b match within a?
    assert failure if not degnucs_contain(a,b)
    '''
    #
    assert len(b)<=len(a), 'degnucs_index: b is not contained in a'
    for offset in range(len(a)-len(b)):
        if all(nucleotide_symbols_match(x,y)
               for x,y in zip(a[offset:], b)):
            return offset
    assert False, 'degnucs_index: b is not contained in a'
    return 0

def reverse_complement( seq ):
    newseq = ''
    L = len(seq)
    for pos in range( L-1, -1, -1 ):
        newseq += base_partner[ seq[ pos ] ]
    assert len( newseq ) == L
    return newseq

