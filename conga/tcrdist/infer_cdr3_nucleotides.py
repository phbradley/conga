import os
from os.path import exists
import pandas as pd
#from . import basic
#from .amino_acids import amino_acids
from .all_genes import all_genes
from .genetic_code import aa2degenerate_codons
from . import tcr_sampler
from . import logo_tools
from . import translation

# need to choose a nucleotide in N insertions
# could get fancier: use nucleotide context, whether v-d d-j or v-j
#  (which determines strand)
#
N_insertion_frequency_order = {'c':0, 'g':1, 't':2, 'a':3}

def get_degenerate_cdr3_nucseqs(cdr3):
    ''' Convert from amino acid cdr3 to list of degenerate coding seqs
    '''

    seqs= ['']
    while cdr3:
        seqs = [x+y for x in seqs for y in aa2degenerate_codons[cdr3[0]]]
        cdr3 = cdr3[1:]

    return seqs

def infer_cdr3_nucleotides(
        organism,
        v,
        j,
        cdr3,
        verbose=False,
        return_stats=False,
        allow_internal_mismatches=False,
        allow_allele_changes=False,
):
    ''' returns cdr3_nucseq and (optionally) the number of insertions

    if allow_allele_changes: returns cdr3_nucseq, jgene, and (optionally) inserts

    Still need to IMPROVE the logic to fixup alleles based on protein seq
    matches?
    '''

    ########################### special handling of common allelic aa variation ########
    if allow_allele_changes: # need to return the j gene
        special_human_jgenes = ['TRAJ24','TRAJ37']
        jg = j.split('*')[0]
        if organism == 'human' and jg in special_human_jgenes:
            best_inserts = 1000
            for allele in ['*01','*02']:
                new_j = jg+allele
                for aim in [True, False]:
                    cdr3_nucseq, total_inserts = infer_cdr3_nucleotides(
                        organism, v, new_j, cdr3, return_stats=True,
                        allow_internal_mismatches=aim,
                        allow_allele_changes=False,
                    )
                    if allow_internal_mismatches:
                        total_inserts += 1
                    #print(j, allow_internal_mismatches, best_inserts, total_inserts,
                    #      cdr3, cdr3_nucseq)
                    if total_inserts < best_inserts:
                        best_inserts = total_inserts
                        new_tcr = (v, new_j, cdr3, cdr3_nucseq)
            _, new_j, _, cdr3_nucseq = new_tcr
        else:
            new_j = j
            cdr3_nucseq, best_inserts = infer_cdr3_nucleotides(
                organism, v, new_j, cdr3, return_stats=True,
                allow_internal_mismatches=allow_internal_mismatches,
                allow_allele_changes=False,
            )

        if return_stats:
            return cdr3_nucseq, new_j, best_inserts
        else:
            return cdr3_nucseq, new_j
    ####################################################################################


    if allow_internal_mismatches:
        mismatch_score = -4 # default for junction analysis
    else:
        mismatch_score = -1000
    cdr3_nucseqs = get_degenerate_cdr3_nucseqs(cdr3)

    sortl = []
    for cdr3_nucseq in cdr3_nucseqs:
        result = tcr_sampler.analyze_junction(
            organism, v, j, cdr3, cdr3_nucseq, mismatch_score=mismatch_score,
            return_cdr3_nucseq_src=True)

        inserts = result[-2]
        assert len(inserts) == 4
        total_insert= sum(inserts[1:])
        if verbose:
            print(cdr3, cdr3_nucseq, total_insert, result[0], result[1])
        # ( v_trim, d0_trim, d1_trim, j_trim ) = trims
        # ( best_d_id, n_vd_insert, n_dj_insert, n_vj_insert ) = inserts

        sortl.append((total_insert, cdr3_nucseq, result))
    sortl.sort()

    # now reconstruct top nucleotide sequence
    v_nucseq = tcr_sampler.get_v_cdr3_nucseq( organism, v )
    j_nucseq = tcr_sampler.get_j_cdr3_nucseq( organism, j )
    total_insert, cdr3_degseq, result = sortl[0]

    cdr3_nucseq = list(cdr3_degseq) # can modify in place

    trims, inserts, cdr3_nucseq_src = result[-3:]

    num_matched_v = cdr3_nucseq_src.count('V')
    num_matched_j = cdr3_nucseq_src.count('J')
    num_matched_d = cdr3_nucseq_src.count('D')

    # V
    for i,(a,b) in enumerate(zip(v_nucseq, cdr3_degseq)):
        if i>=num_matched_v:break
        if logo_tools.nucleotide_symbols_match(a,b):
            cdr3_nucseq[i] = a
        else:
            assert allow_internal_mismatches

    # J
    for i,(a,b) in enumerate(zip(reversed(j_nucseq), reversed(cdr3_degseq))):
        if i>=num_matched_j:break
        if logo_tools.nucleotide_symbols_match(a,b):
            cdr3_nucseq[-1-i] = a
        else:
            assert allow_internal_mismatches

    ( v_trim, d0_trim, d1_trim, j_trim ) = trims
    ( best_d_id, n_vd_insert, n_dj_insert, n_vj_insert ) = inserts

    if num_matched_d: # D
        assert best_d_id>0
        assert all_genes[organism][v].chain == 'B'
        d_nucseq = tcr_sampler.all_trbd_nucseq[organism][best_d_id]
        match_start = cdr3_nucseq_src.index('D')
        d_match_seq = d_nucseq[d0_trim:None if d1_trim==0 else -d1_trim]
        assert len(d_match_seq) == num_matched_d
        for i,(a,b) in enumerate(zip(d_match_seq, cdr3_degseq[match_start:])):
            if logo_tools.nucleotide_symbols_match(a,b):
                cdr3_nucseq[match_start+i] = a
            else:
                assert allow_internal_mismatches

    # there may still be degenerate symbols in the N-inserted regions
    for i in range(len(cdr3_nucseq)):
        n = cdr3_nucseq[i]
        if n not in 'acgt':
            cdr3_nucseq[i] = min(logo_tools.nucleotide_classes_lower_case[n],
                                 key=lambda x:N_insertion_frequency_order[x])

    cdr3_nucseq = ''.join(cdr3_nucseq)
    assert all(x in 'acgt' for x in cdr3_nucseq)
    assert translation.get_translation(cdr3_nucseq) == cdr3

    if return_stats:
        return cdr3_nucseq, total_insert
    else:
        return cdr3_nucseq






