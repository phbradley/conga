#from basic import *
from . import logo_tools
from .genetic_code import genetic_code


bases_plus = list(logo_tools.nucleotide_classes_lower_case.keys())

for a in bases_plus:
    for b in bases_plus:
        for c in bases_plus:
            codon = a+b+c
            if codon in genetic_code: continue

            aas = []
            for a1 in logo_tools.nucleotide_classes_lower_case[a]:
                for b1 in logo_tools.nucleotide_classes_lower_case[b]:
                    for c1 in logo_tools.nucleotide_classes_lower_case[c]:
                        aas.append( genetic_code[ a1+b1+c1 ] )
            if min(aas) == max(aas):
                genetic_code[codon] = aas[0]
            else:
                genetic_code[codon] = 'X'



def get_translation( seq, frame='+1' ): ## STUPID-- frame is 1/2/3 with a +/- in front
    assert frame[0] in '+-'
    if frame[0] == '-': seq = logo_tools.reverse_complement( seq )
    offset = abs( int( frame ))-1
    assert offset in range(3)
    seq = seq[offset:].lower()
    naa = len(seq)//3
    protseq = ''
    for i in range(naa):
        codon = seq[3*i:3*i+3]
        if '#' in codon:
            protseq += '#'
        else:
            protseq += genetic_code.get( codon, 'X' )
    return protseq

def get_translation_and_codons( seq, frame='+1' ): ## STUPID-- frame is 1/2/3 with a +/- in front
    assert frame[0] in '+-'
    if frame[0] == '-': seq = logo_tools.reverse_complement( seq )
    offset = abs( int( frame ))-1
    assert offset in range(3)
    seq = seq[offset:].lower()
    naa = len(seq)//3
    protseq = ''
    codons = []
    for i in range(naa):
        codon = seq[3*i:3*i+3]
        codons.append( codon )
        if '#' in codon:
            protseq += '#'
        else:
            protseq += genetic_code.get( codon, 'X' )
    return protseq,codons

