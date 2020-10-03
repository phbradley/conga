import os
from os.path import exists
import pandas as pd
from . import basic
from .amino_acids import amino_acids
from .tcr_distances_blosum import blosum
from . import translation

cdrs_sep = ';'
gap_character = '.'

all_genes = {}

class TCR_Gene:
    def __init__( self, l ):
        # l comes from pandas dataframe.itertuples
        self.id = l.id
        self.organism = l.organism
        self.chain = l.chain
        self.region = l.region
        self.nucseq = l.nucseq
        self.alseq = l.aligned_protseq
        if pd.isna(l.cdrs):
            self.cdrs = []
            self.cdr_columns = []
        else:
            self.cdrs = l.cdrs.split(cdrs_sep)
            ## these are still 1-indexed !!!!!!!!!!!!!!
            self.cdr_columns = [ list(map( int,x.split('-'))) for x in l.cdr_columns.split(cdrs_sep) ]
        frame = l.frame
        #assert frame in ['+1','+2','+3','1','2','3']
        assert frame in [1,2,3] # now parsed by pandas, so string converted to int
        #self.nucseq_offset = int( frame[-1] )-1 ## 0, 1 or 2 (0-indexed for python)
        self.nucseq_offset = frame-1 ## 0, 1 or 2 (0-indexed for python)
        self.protseq = translation.get_translation( self.nucseq, f'+{frame}' )
        assert self.protseq == self.alseq.replace(gap_character,'')
        # sanity check
        if self.cdrs:
            assert self.cdrs == [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]

def trim_allele_to_gene( id ):
    return id[: id.index('*') ] #will fail if id doesn't contain '*'

db_file = os.path.dirname(os.path.realpath(__file__))+'/db/'+basic.db_file
assert exists(db_file)

df = pd.read_csv(db_file, sep='\t')

for l in df.itertuples():
    g = TCR_Gene( l )
    if g.organism not in all_genes:
        all_genes[g.organism] = {} # map from id to TCR_Gene objects
    all_genes[g.organism][g.id] = g


verbose = ( __name__ == '__main__' )

for organism,genes in all_genes.items():

    for ab in 'AB':
        org_merged_loopseqs = {}
        for id,g in genes.items():
            if g.chain == ab and g.region == 'V':
                loopseqs = g.cdrs[:-1] ## exclude CDR3 Nterm
                org_merged_loopseqs[id] = ' '.join( loopseqs )

        all_loopseq_nbrs = {}
        all_loopseq_nbrs_mm1 = {}
        for id1,seq1 in org_merged_loopseqs.items():
            g1 = genes[id1]
            cpos = g1.cdr_columns[-1][0] - 1 #0-indexed
            alseq1 = g1.alseq
            minlen = cpos+1
            assert len(alseq1) >= minlen
            if alseq1[cpos] != 'C' and verbose:
                print('funny cpos:',id1,alseq1,g1.cdrs[-1])

            all_loopseq_nbrs[id1] = []
            all_loopseq_nbrs_mm1[id1] = []
            for id2,seq2 in org_merged_loopseqs.items():
                g2 = genes[id2]
                alseq2 = g2.alseq
                assert len(alseq2) >= minlen
                assert len(seq1) == len(seq2)
                if seq1 == seq2:
                    all_loopseq_nbrs[id1].append( id2 )
                    all_loopseq_nbrs_mm1[id1].append( id2 )
                    continue

                ## count mismatches between these two, maybe count as an "_mm1" nbr
                loop_mismatches = 0
                loop_mismatches_cdrx = 0
                loop_mismatch_seqs =[]
                spaces=0
                for a,b in zip( seq1,seq2):
                    if a==' ':
                        spaces+=1
                        continue
                    if a!= b:
                        if a in '*.' or b in '*.':
                            loop_mismatches += 10
                            break
                        else:
                            if not (a in amino_acids and b in amino_acids):
                                print( id1,id2,a,b)
                            assert a in amino_acids and b in amino_acids
                            if spaces<=1:
                                loop_mismatches += 1
                                loop_mismatch_seqs.append( ( a,b ) )
                            else:
                                assert spaces==2
                                loop_mismatches_cdrx += 1
                            if loop_mismatches>1:
                                break
                if loop_mismatches <=1:
                    all_mismatches = 0
                    for a,b in zip( alseq1[:cpos+2],alseq2[:cpos+2]):
                        if a!= b:
                            if a in '*.' or b in '*.':
                                all_mismatches += 10
                            else:
                                if not (a in amino_acids and b in amino_acids):
                                    print( id1,id2,a,b)
                                assert a in amino_acids and b in amino_acids
                                all_mismatches += 1
                    #dist = tcr_distances.blosum_sequence_distance( seq1, seq2, gap_penalty=10 )
                    if loop_mismatches<=1 and loop_mismatches + loop_mismatches_cdrx <= 2 and all_mismatches<=10:
                        if loop_mismatches == 1:
                            blscore= blosum[(loop_mismatch_seqs[0][0],loop_mismatch_seqs[0][1])]
                        else:
                            blscore = 100
                        if blscore>=1:
                            all_loopseq_nbrs_mm1[id1].append( id2 )
                            if loop_mismatches>0 and verbose:
                                mmstring = ','.join(['%s/%s'%(x[0],x[1]) for x in loop_mismatch_seqs])
                                gene1 = trim_allele_to_gene( id1 )
                                gene2 = trim_allele_to_gene( id2 )
                                if gene1 != gene2 and verbose:
                                    print('v_mismatches:',organism,mmstring,blscore,id1,id2,\
                                        loop_mismatches,loop_mismatches_cdrx,all_mismatches,seq1)
                                    print('v_mismatches:',organism,mmstring,blscore,id1,id2,\
                                        loop_mismatches,loop_mismatches_cdrx,all_mismatches,seq2)


        for id in all_loopseq_nbrs:
            rep = min( all_loopseq_nbrs[id] )
            assert org_merged_loopseqs[id] == org_merged_loopseqs[ rep ]
            genes[id].rep = rep
            if verbose:
                print('vrep %s %15s %15s %s'%(organism, id, rep, org_merged_loopseqs[id]))


        ## merge mm1 nbrs to guarantee transitivity
        while True:
            new_nbrs = False
            for id1 in all_loopseq_nbrs_mm1:
                new_id1_nbrs = False
                for id2 in all_loopseq_nbrs_mm1[id1]:
                    for id3 in all_loopseq_nbrs_mm1[id2]:
                        if id3 not in all_loopseq_nbrs_mm1[id1]:
                            all_loopseq_nbrs_mm1[id1].append( id3 )
                            if verbose:
                                print('new_nbr:',id1,'<--->',id2,'<--->',id3)
                            new_id1_nbrs = True
                            break
                    if new_id1_nbrs:
                        break
                if new_id1_nbrs:
                    new_nbrs = True
            if verbose:
                print('new_nbrs:',ab,organism,new_nbrs)
            if not new_nbrs:
                break

        for id in all_loopseq_nbrs_mm1:
            rep = min( all_loopseq_nbrs_mm1[id] )
            genes[id].mm1_rep = rep
            if verbose:
                print('mm1vrep %s %15s %15s %s'%(organism, id, rep,org_merged_loopseqs[id]))


    ## setup Jseq reps
    for ab in 'AB':
        jloopseqs = {}
        for id,g in genes.items():
            if g.chain == ab and g.region == 'J':
                num = len( g.cdrs[0].replace( gap_character, '' ) )
                jloopseq = g.protseq[:num+3] ## go all the way up to and including the GXG
                jloopseqs[id] = jloopseq
        all_jloopseq_nbrs = {}
        for id1,seq1 in jloopseqs.items():
            all_jloopseq_nbrs[id1] = []
            for id2,seq2 in jloopseqs.items():
                if seq1 == seq2:
                    all_jloopseq_nbrs[id1].append( id2 )
        for id in all_jloopseq_nbrs:
            rep = min( all_jloopseq_nbrs[id] )
            genes[id].rep = rep
            genes[id].mm1_rep = rep # just so we have an mm1_rep field defined...
            assert jloopseqs[id] == jloopseqs[ rep ]
            if verbose:
                print('jrep %s %15s %15s %15s'%(organism, id, rep, jloopseqs[id]))



    ## setup a mapping that we can use for counting when allowing mm1s and also ignoring alleles

    # allele2mm1_rep_gene_for_counting = {}
    # def get_mm1_rep_ignoring_allele( gene, organism ): # helper fxn
    #     rep = get_mm1_rep( gene, organism )
    #     rep = rep[:rep.index('*')]
    #     return rep

    #allele2mm1_rep_gene_for_counting[ organism ] = {}

    if not basic.CLASSIC_COUNTREPS:
        # simpler scheme for choosing the 'count_rep' field
        for id, g in all_genes[organism].items():
            g.count_rep = trim_allele_to_gene(id)
    else:
        for chain in 'AB':
            for vj in 'VJ':
                allele_gs = [ (id,g) for (id,g) in all_genes[organism].items() if g.chain==chain and g.region==vj]

                gene2rep = {}
                gene2alleles = {}
                rep_gene2alleles = {}

                for allele,g in allele_gs:
                    #assert allele[2] == chain
                    gene = trim_allele_to_gene( allele )
                    rep_gene = trim_allele_to_gene( g.mm1_rep )
                    if rep_gene not in rep_gene2alleles:
                        rep_gene2alleles[ rep_gene ] = []
                    rep_gene2alleles[ rep_gene ].append( allele )

                    if gene not in gene2rep:
                        gene2rep[gene] = set()
                        gene2alleles[gene] = []
                    gene2rep[ gene ].add( rep_gene )
                    gene2alleles[gene].append( allele )

                merge_rep_genes = {}
                for gene,reps in gene2rep.items():
                    if len(reps)>1:
                        if verbose:
                            print('multireps:',organism, gene, reps)
                            for allele in gene2alleles[gene]:
                                print(' '.join(all_genes[organism][allele].cdrs), allele, \
                                    all_genes[organism][allele].rep, \
                                    all_genes[organism][allele].mm1_rep)
                        assert vj=='V'

                        ## we are going to merge these reps
                        ## which one should we choose?
                        l = [ (len(rep_gene2alleles[rep]), rep ) for rep in reps ]
                        l.sort()
                        l.reverse()
                        assert l[0][0] > l[1][0]
                        toprep = l[0][1]
                        for (count,rep) in l:
                            if rep in merge_rep_genes:
                                # ACK need to think more about this, should probably just kill this logic!
                                assert rep == toprep and merge_rep_genes[rep] == rep
                            merge_rep_genes[ rep ] = toprep


                for allele,g in allele_gs:
                    count_rep = trim_allele_to_gene( g.mm1_rep ) #get_mm1_rep_ignoring_allele( allele, organism )
                    if count_rep in merge_rep_genes:
                        count_rep = merge_rep_genes[ count_rep ]
                    g.count_rep = count_rep #allele2mm1_rep_gene_for_counting[ organism ][ allele] = count_rep
                    if verbose:
                        print('countrep:',organism, allele, count_rep)


if __name__ == '__main__':
    for org, genes in all_genes.items():
        for id, g in genes.items():
            print( org, id, g.cdrs )

