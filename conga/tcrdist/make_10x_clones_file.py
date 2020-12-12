from .basic import *
from .all_genes import all_genes
#import parse_tsv
from collections import Counter
from itertools import chain
import sys
import pandas as pd
from ..util import organism2vdj_type, IG_VDJ_TYPE

MIN_CDR3_LEN = 6 # otherwise tcrdist barfs; actually we could also set this to 5 and be OK

def show(tcr):
    "For debugging"
    if type( tcr[0] ) is str:
        return ' '.join(tcr[:3])
    else:
        return ' '.join( show(x) for x in tcr )

def fixup_gene_name( gene, gene_suffix, expected_gene_names ):

    if gene in expected_gene_names:
        return gene # ALL DONE

    if '*' not in gene:
        gene += gene_suffix

    if gene in expected_gene_names:
        return gene # ALL DONE

    vj = gene[3]
    assert vj in 'VJ'

    if vj=='V' and 'DV' in gene:
        # sometimes they just delete the '/'
        # sometimes they convert it into a '-'
        new_gene = gene[:gene.index('DV')]+'/'+gene[gene.index('DV'):]
        if new_gene in expected_gene_names:
            gene = new_gene
        else:
            new_gene = gene[:gene.index('DV')-1]+'/'+gene[gene.index('DV'):]
            if new_gene in expected_gene_names:
                gene = new_gene
            else:
                print('trouble parsing V gene with DV in it:', gene)

    return gene # may still not be in expected_gene_names, will check for that later

def get_ab_from_10x_chain(chain, organism):
    ''' Returns None if the chain is not valid for this 'organism'
    '''
    if organism in ['human', 'mouse']:
        if chain in ['TRA','TRB']:
            return chain[2]
        else:
            return None
    elif organism in ['human_gd','mouse_gd']:
        if chain in ['TRA','TRG','TRD']:
            return 'A' if chain=='TRG' else 'B'
        else:
            return None
    elif organism in ['human_ig','mouse_ig']:
        if chain in ['IGH', 'IGK', 'IGL']:
            return 'B' if chain=='IGH' else 'A'
        else:
            return None
    else:
        print('unrecognized organism in get_ab_from_10x_chain:', organism)
        sys.exit()


def read_tcr_data(
        organism,
        contig_annotations_csvfile,
        consensus_annotations_csvfile,
        allow_unknown_genes = False,
        verbose = False,
        prefix_clone_ids_with_tcr_type = False,
):
    """ Parse tcr data, only taking 'productive' tcrs

    Returns:

    clonotype2tcrs, clonotype2barcodes

    """

    expected_gene_names = set(all_genes[organism].keys())

    #from cdr3s_human import all_align_fasta

    gene_suffix = '*01' # may not be used


    if prefix_clone_ids_with_tcr_type:
        if organism2vdj_type[organism] == IG_VDJ_TYPE:
            clone_id_prefix = 'bcr_'
        else:
            clone_id_prefix = 'tcr_'
    else:
        clone_id_prefix = ''

    # read the contig annotations-- map from clonotypes to barcodes
    # barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id
    # AAAGATGGTCTTCTCG-1,True,AAAGATGGTCTTCTCG-1_contig_1,True,695,TRB,TRBV5-1*01,TRBD2*02,TRBJ2-3*01,TRBC2*01,True,True,CASSPLAGYAADTQYF,TGCGCCAGCAGCCCCCTAGCGGGATACGCAGCAGATACGCAGTATTTT,9427,9,clonotype14,clonotype14_consensus_1
    assert exists( contig_annotations_csvfile )

    #_, lines = parse_csv_file(contig_annotations_csvfile)
    df = pd.read_csv(contig_annotations_csvfile)
    df['productive'] = df['productive'].astype(str) #sometimes it already is if there are 'Nones' in there...
    clonotype2barcodes = {}
    clonotype2tcrs_backup = {} ## in case we dont have a consensus_annotations_csvfile
    for l in df.itertuples():
        bc = l.barcode
        clonotype = clone_id_prefix + l.raw_clonotype_id
        # annoying: pandas sometimes converts to True/False booleans and sometimes not.
        assert l.productive in [ 'None', 'False', 'True']
        if clonotype =='None':
            continue
        if clonotype not in clonotype2barcodes:
            clonotype2barcodes[clonotype] = []
        if bc in clonotype2barcodes[clonotype]:
            pass
            #print 'repeat barcode'
        else:
            clonotype2barcodes[clonotype].append( bc )

        ## experimenting here ########################################3
        if l.productive != 'True':
            continue
        if l.cdr3.lower() == 'none' or l.cdr3_nt.lower() == 'none':
            continue

        chain = l.chain
        ab = get_ab_from_10x_chain(chain, organism)
        if ab is None:
            continue
        if clonotype not in clonotype2tcrs_backup:
            clonotype2tcrs_backup[ clonotype ] = {'A':Counter(), 'B':Counter() }
        # stolen from below
        vg = fixup_gene_name(l.v_gene, gene_suffix, expected_gene_names)
        jg = fixup_gene_name(l.j_gene, gene_suffix, expected_gene_names)

        if vg not in expected_gene_names:
            print('unrecognized V gene:', organism, vg)
            if not allow_unknown_genes:
                continue
        if jg not in expected_gene_names:
            print('unrecognized J gene:', organism, jg)
            if not allow_unknown_genes:
                continue
        #assert vg in all_align_fasta[organism]
        #assert jg in all_align_fasta[organism]

        tcr_chain = ( vg, jg, l.cdr3, l.cdr3_nt.lower() )

        clonotype2tcrs_backup[clonotype][ab][tcr_chain] += int(l.umis)

    for id in clonotype2tcrs_backup:
        for ab in 'AB':
            for t1,count1 in clonotype2tcrs_backup[id][ab].items():
                for t2, count2 in clonotype2tcrs_backup[id][ab].items():
                    if t2<=t1:continue
                    if t1[3] == t2[3]:
                        print('repeat??', count1, count2, t1, t2)



    if consensus_annotations_csvfile is None:
        clonotype2tcrs = clonotype2tcrs_backup
    else:

        ## now read details on the individual chains for each clonotype
        # ==> tcr/human/JCC176_TX2_TCR_consensus_annotations.csv <==
        # clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis
        # clonotype100,clonotype100_consensus_1,550,TRB,TRBV24-1*01,TRBD1*01,TRBJ2-7*01,TRBC2*01,True,True,CATSDPGQGGYEQYF,TGTGCCACCAGTGACCCCGGACAGGGAGGATACGAGCAGTACTTC,8957,9

        assert exists(consensus_annotations_csvfile)
        df = pd.read_csv( consensus_annotations_csvfile )
        df['productive'] = df['productive'].astype(str) #sometimes it already is if there are 'Nones' in there...


        ## first get clonotypes with one alpha and one beta
        clonotype2tcrs = {}

        for l in df.itertuples():
            assert l.productive in [ 'None', 'False', 'True']
            if l.productive == 'True':
                id = l.clonotype_id
                if id not in clonotype2tcrs:
                    # dictionaries mapping from tcr to umi-count
                    clonotype2tcrs[id] = { 'A':Counter(), 'B':Counter() } #, 'G':[], 'D': [] }
                    assert id in clonotype2barcodes

                ch = l.chain
                ab = get_ab_from_10x_chain(ch, organism)
                if ab is None:
                    print('skipline:', consensus_annotations_csvfile, ch, l.v_gene, l.j_gene)
                    continue
                vg = fixup_gene_name(l.v_gene, gene_suffix, expected_gene_names)
                jg = fixup_gene_name(l.j_gene, gene_suffix, expected_gene_names)

                if vg not in expected_gene_names:
                    print('unrecognized V gene:', organism, vg)
                    if not allow_unknown_genes:
                        continue
                if jg not in expected_gene_names:
                    print('unrecognized J gene:', organism, jg)
                    if not allow_unknown_genes:
                        continue
                #assert vg in all_align_fasta[organism]
                #assert jg in all_align_fasta[organism]
                tcr_chain = ( vg, jg, l.cdr3, l.cdr3_nt.lower() )

                if tcr_chain not in clonotype2tcrs[id][ab]:
                    umis = int( l.umis )
                    clonotype2tcrs[id][ab][ tcr_chain ] = umis
                    old_umis = clonotype2tcrs_backup[id][ab][tcr_chain]
                    if umis != old_umis:
                        print('diff_umis:',umis, old_umis, id,ab,tcr_chain)
                else:
                    print('repeat?',id,ab,tcr_chain)
            # else:
            #     if l.productive not in [ 'None','False' ]:
            #         print('unproductive?',l.productive)


        if verbose:
            idl1 = sorted( clonotype2tcrs_backup.keys())
            idl2 = sorted( clonotype2tcrs.keys())
            print('same ids:', len(idl1), len(idl2), idl1==idl2)
            for id in clonotype2tcrs_backup:
                if id in clonotype2tcrs:
                    for ab in 'AB':
                        tl1 = sorted(clonotype2tcrs_backup[id][ab].keys())
                        tl2 = sorted(clonotype2tcrs[id][ab].keys())
                        if tl1 != tl2:
                            print('diffids:',id,ab,tl1,tl2)

    return clonotype2tcrs, clonotype2barcodes

def read_tcr_data_batch(
        organism,
        metadata_file,
        allow_unknown_genes = False,
        verbose = False,
        prefix_clone_ids_with_tcr_type = False,
):
    """ Parse tcr data, only taking 'productive' tcrs

    Returns:

    clonotype2tcrs, clonotype2barcodes

    """
    assert exists( metadata_file )

    md = pd.read_csv(metadata_file, dtype=str)#"string")

    if prefix_clone_ids_with_tcr_type:
        if organism2vdj_type[organism] == IG_VDJ_TYPE:
            clone_id_prefix = 'bcr_'
        else:
            clone_id_prefix = 'tcr_'
    else:
        clone_id_prefix = ''

    # read in contig files and update suffix to match GEX matrix
    contig_list = []
    for x in range(len(md['file'])):
        dfx = pd.read_csv( md.loc[ x , 'file'] ) #read contig_df

        suffix = md.loc[ x , 'suffix'] # new suffix
        barcodes = dfx['barcode'].str.split('-').str.get(0)

        #hacky
        dfx['barcode'] = barcodes + '-' + suffix
        dfx['contig_id'] = barcodes + '-' + suffix + '_' + \
            dfx['contig_id'].str.split('_').str.get(1) + \
            '_' + dfx['contig_id'].str.split('_').str.get(2) # currently unused, but can't hurt

        # giving each library a tag here really boosted the number of clones I got back
        dfx['raw_clonotype_id'] = clone_id_prefix + dfx['raw_clonotype_id'] + '_' + suffix
        dfx['raw_consensus_id'] = clone_id_prefix + dfx['raw_consensus_id'] + '_' + suffix # currently unused, can't hurt

        contig_list.append(dfx)

    df = pd.concat(contig_list)

    expected_gene_names = set(all_genes[organism].keys())

    #from cdr3s_human import all_align_fasta

    gene_suffix = '*01' # may not be used


    # read the contig annotations-- map from clonotypes to barcodes
    # barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id
    # AAAGATGGTCTTCTCG-1,True,AAAGATGGTCTTCTCG-1_contig_1,True,695,TRB,TRBV5-1*01,TRBD2*02,TRBJ2-3*01,TRBC2*01,True,True,CASSPLAGYAADTQYF,TGCGCCAGCAGCCCCCTAGCGGGATACGCAGCAGATACGCAGTATTTT,9427,9,clonotype14,clonotype14_consensus_1

    #_, lines = parse_csv_file(contig_annotations_csvfile)
    #df = pd.read_csv(contig_annotations_csvfile)
    df['productive'] = df['productive'].astype(str) #sometimes it already is if there are 'Nones' in there...
    clonotype2tcrs = {}
    clonotype2barcodes = {}
    for l in df.itertuples():
        # the fields we use:   barcode  raw_clonotype_id  productive  cdr3  cdr3_nt  chain  v_gene  j_gene  umis
        bc = l.barcode
        clonotype = l.raw_clonotype_id
        # annoying: pandas sometimes converts to True/False booleans and sometimes not.
        assert l.productive in [ 'None', 'False', 'True']
        if clonotype =='None':
            continue
        if clonotype not in clonotype2barcodes:
            clonotype2barcodes[clonotype] = []
        if bc in clonotype2barcodes[clonotype]:
            pass
            #print 'repeat barcode'
        else:
            clonotype2barcodes[clonotype].append( bc )

        ## experimenting here ########################################3
        if l.productive != 'True':
            continue
        if l.cdr3.lower() == 'none' or l.cdr3_nt.lower() == 'none':
            continue

        chain = l.chain
        ab = get_ab_from_10x_chain(chain, organism)
        if ab is None:
            continue
        if clonotype not in clonotype2tcrs:
            clonotype2tcrs[ clonotype ] = {'A':Counter(), 'B':Counter() }
        # stolen from below
        vg = fixup_gene_name(l.v_gene, gene_suffix, expected_gene_names)
        jg = fixup_gene_name(l.j_gene, gene_suffix, expected_gene_names)

        if vg not in expected_gene_names:
            print('unrecognized V gene:', organism, vg)
            if not allow_unknown_genes:
                continue
        if jg not in expected_gene_names:
            print('unrecognized J gene:', organism, jg)
            if not allow_unknown_genes:
                continue
        #assert vg in all_align_fasta[organism]
        #assert jg in all_align_fasta[organism]

        tcr_chain = ( vg, jg, l.cdr3, l.cdr3_nt.lower() )

        clonotype2tcrs[clonotype][ab][tcr_chain] += int(l.umis)

    for id in clonotype2tcrs:
        for ab in 'AB':
            for t1,count1 in clonotype2tcrs[id][ab].items():
                for t2, count2 in clonotype2tcrs[id][ab].items():
                    if t2<=t1:continue
                    if t1[3] == t2[3]:
                        print('repeat??', count1, count2, t1, t2)


    return clonotype2tcrs, clonotype2barcodes

def _make_clones_file( organism, outfile, clonotype2tcrs, clonotype2barcodes, verbose=False ):
    ''' Make a clones file with information parsed from the 10X csv files

    organism is one of ['mouse','human']

    outfile is the name of the clones file to be created

    '''

    #tmpfile = outfile+'.tmp' # a temporary intermediate file
    bc_mapfile = outfile+'.barcode_mapping.tsv'
    outmap = open(bc_mapfile,'w')
    outmap.write('clone_id\tbarcodes\n')

    outfields = 'clone_id subject clone_size va_gene ja_gene va2_gene ja2_gene vb_gene jb_gene cdr3a cdr3a_nucseq cdr3a2 cdr3a2_nucseq cdr3b cdr3b_nucseq'\
        .split()
    extra_fields = 'alpha_umi alpha2_umi beta_umi num_alphas num_betas'.split()
    outfields += extra_fields

    # we used to make a temporary file and then run tcr-dist/file_converter.py on it
    # for this slim version, just take the temporary file
    #out = open(tmpfile,'w')
    out = open(outfile,'w')
    out.write('\t'.join( outfields )+'\n' )

    for clonotype in sorted(clonotype2tcrs.keys()):
        tcrs = clonotype2tcrs[clonotype]
        if len(tcrs['A']) >= 1 and len(tcrs['B']) >= 1:
            atcrs = tcrs['A'].most_common()
            btcrs = tcrs['B'].most_common()
            if len(atcrs)>1:
                if verbose:
                    print('multiple alphas, picking top umi:',' '.join( str(x) for _,x in atcrs ))
                atcr2, atcr2_umi = atcrs[1]
            else:
                atcr2, atcr2_umi = ('', '', '', ''), 0
            if len(btcrs)>1 and verbose:
                print('multiple  betas, picking top umi:',' '.join( str(x) for _,x in btcrs ))
            atcr, atcr_umi = atcrs[0]
            btcr, btcr_umi = btcrs[0]
            outl = {}
            outl['clone_id'] = clonotype
            outl['subject'] = 'UNK_S'
            outl['clone_size'] = len(clonotype2barcodes[clonotype])
            outl['va_gene']      = atcr[0]
            outl['ja_gene']      = atcr[1]
            outl['cdr3a']        = atcr[2]
            outl['cdr3a_nucseq'] = atcr[3]
            outl['alpha_umi']    = str(atcr_umi)
            outl['va2_gene']      = atcr2[0]
            outl['ja2_gene']      = atcr2[1]
            outl['cdr3a2']        = atcr2[2]
            outl['cdr3a2_nucseq'] = atcr2[3]
            outl['alpha2_umi']    = str(atcr2_umi)
            outl['num_alphas']   = str(len(atcrs))
            outl['vb_gene']      = btcr[0]
            outl['jb_gene']      = btcr[1]
            outl['cdr3b']        = btcr[2]
            outl['cdr3b_nucseq'] = btcr[3]
            outl['beta_umi']     = str(btcr_umi)
            outl['num_betas']    = str(len(btcrs))
            if len(outl['cdr3a']) < MIN_CDR3_LEN or  len(outl['cdr3b']) < MIN_CDR3_LEN:
                if verbose:
                    print('Warning: skipping clonotype with short cdr3s: {} {} {}'\
                          .format(clonotype, outl['cdr3a'], outl['cdr3b'] ))
                continue
            out.write( '\t'.join(str(outl[x]) for x in outfields)+'\n')
            outmap.write('{}\t{}\n'.format(clonotype,','.join(clonotype2barcodes[clonotype])))
    out.close()
    outmap.close()

    # for the time being, go with the clones file we just made even though it doesn't have all the usual stupid info


def setup_filtered_clonotype_dicts(
        clonotype2tcrs,
        clonotype2barcodes,
        min_repeat_count_fraction = 0.33,
        verbose = False
):
    ''' returns new_clonotype2tcrs, new_clonotype2barcodes
    '''

    # get count of how many cells support each pairing
    #
    pairing_counts, pairing_counts_by_umi = Counter(), Counter()
    for cid, tcrs in clonotype2tcrs.items():
        size = len(clonotype2barcodes[cid])
        for t1, umi1 in tcrs['A'].items():
            for t2, umi2 in tcrs['B'].items():
                pairing_counts[ (t1,t2) ] += size
                pairing_counts_by_umi[ (t1,t2) ] += min(umi1, umi2)

    # this is no longer going to be a bijection!
    chain_partner = {}

    valid_ab_pairings = set()

    for ( (t1,t2), count ) in list(pairing_counts.most_common()): # make a copy for sanity
        if t1 in chain_partner or t2 in chain_partner:
            t1_p2 = chain_partner.get(t1,None)
            t2_p2 = chain_partner.get(t2,None)

            oldcount1 = pairing_counts[ (t1,t1_p2) ]
            oldcount2 = pairing_counts[ (t2_p2,t2) ]

            if count >= min_repeat_count_fraction * max(oldcount1,oldcount2):
                # take it anyway -- second alpha or genuinely shared alpha?

                if verbose:
                    print('take_rep_partners:', count, oldcount1, oldcount2, t1, t2, t1_p2, t2_p2)
                # dont overwrite the old pairing... might not do either of these!
                if t1 not in chain_partner:
                    chain_partner[t1] = t2
                if t2 not in chain_partner:
                    chain_partner[t2] = t1

                valid_ab_pairings.add( (t1, t2 ) )
            else:
                if verbose:
                    print('skip_rep_partners:', count, oldcount1, oldcount2, t1, t2, t1_p2, t2_p2)

        else: # neither chain already in chain_partner
            #
            # NOTE: removed the code checking for TIES!!!!!!!!!!!
            if verbose:
                print('norep_partners:', count, t1, t2)
            chain_partner[t1] = t2
            chain_partner[t2] = t1

            valid_ab_pairings.add( ( t1, t2 ) )


    # now let's revisit the clonotypes
    pairs_tuple2clonotypes = {}
    ab_counts = Counter() # for diagnostics
    for (clone_size, cid) in reversed( sorted( (len(y), x) for x,y in clonotype2barcodes.items() ) ):
        if cid not in clonotype2tcrs:
            #print('WHOAH missing tcrs for clonotype', clone_size, cid, clonotype2barcodes[cid])
            continue
        tcrs = clonotype2tcrs[cid]
        was_good_clone = len(tcrs['A']) >= 1 and len(tcrs['B']) >= 1
        pairs = []
        for atcr in tcrs['A']:
            for btcr in tcrs['B']:
                if ( atcr,btcr ) in valid_ab_pairings:
                    pairs.append( (atcr,btcr) )
        if pairs:
            pairs_tuple = None
            if len(pairs)>1:
                alphas = set( x[0] for x in pairs )
                betas = set( x[1] for x in pairs )
                ab_counts[ (len(alphas),len(betas) ) ] += clone_size

                if len(alphas) == 2 and len(betas) == 1: ## only allow double-alphas
                    assert len(pairs) == 2
                    pairs_tuple = tuple(x[1] for x in reversed(sorted( [ (pairing_counts[x],x) for x in pairs ] )))
                    assert len( pairs_tuple)==2
                    # confirm ordered by counts
                    assert pairing_counts[pairs_tuple[0]] >= pairing_counts[pairs_tuple[1]]
            else:
                ab_counts[ (1,1) ] += clone_size
                pairs_tuple = tuple(pairs)
                assert len(pairs_tuple) == 1

            if pairs_tuple:
                pairs_tuple2clonotypes.setdefault( pairs_tuple, [] ).append( cid )
            else:
                if verbose:
                    print('SKIPCLONE:', was_good_clone, cid, clone_size, pairs, 'bad_pairs')
        else:
            if verbose:
                print('SKIPCLONE:', was_good_clone, cid, clone_size, 'no_valid_pairs')


    ## reorder pairs_tuples in the case of ties, using umis
    reorder = []
    for pt, clonotypes in pairs_tuple2clonotypes.items():
        assert len(pt) in [1,2]
        if len(pt) == 2:
            assert pt[0][1] == pt[1][1] # same beta chain
            at1, bt = pt[0]
            at2, _  = pt[1]
            count1, count2 = pairing_counts[(at1,bt)], pairing_counts[(at2,bt)]
            assert count1 >= count2
            if count1 == count2:
                # no way to figure out which is the primary alpha!
                # look at umis?
                c1 = sum( clonotype2tcrs[x]['A'][at1] for x in clonotypes )
                c2 = sum( clonotype2tcrs[x]['A'][at2] for x in clonotypes )
                if verbose:
                    print('alphatie:', count1, count2, c1, c2, show(at1), show(at2), show(bt))
                if c2 > c1:
                    reorder.append( pt )

    for pt in reorder:
        #print('reorder:', show(pt))
        assert len(pt) == 2
        rpt = (pt[1], pt[0])
        assert pt in pairs_tuple2clonotypes and rpt not in pairs_tuple2clonotypes
        pairs_tuple2clonotypes[rpt] = pairs_tuple2clonotypes[pt][:]
        del pairs_tuple2clonotypes[pt]


    ## look for len1 pairs_tuples that overlap with two different len2 pairs_tuples
    merge_into_pairs = []
    for pt1 in pairs_tuple2clonotypes:
        if len(pt1) == 1:
            overlaps = []
            for pt2 in pairs_tuple2clonotypes:
                if len(pt2) == 2 and pt2[0] == pt1[0]:
                    overlaps.append( pt2 )
                elif len(pt2) == 2 and pt2[1] == pt1[0]:
                    # this 'elif' was added 2020-12-12; before that we were missing some
                    # overlaps...
                    # print('MISSED OVERLAP:', pairing_counts[pt2[0]], pairing_counts[pt2[1]],
                    #       pairing_counts_by_umi[pt2[0]], pairing_counts_by_umi[pt2[1]],
                    #       show(pt1), show(pt2))
                    overlaps.append( pt2 )
            if len(overlaps)>1:
                if verbose:
                    print('badoverlaps:', len(overlaps), show(pt1), show(overlaps))
            elif len(overlaps)==1:
                pt2 = overlaps[0]
                merge_into_pairs.append( (pt1,pt2) )

    for pt1, pt2 in merge_into_pairs:
        assert len(pt1) == 1 and len(pt2) == 2
        #print('mergeinto:', show(pt1), show(pt2))
        pairs_tuple2clonotypes[pt2].extend( pairs_tuple2clonotypes[pt1] )
        del pairs_tuple2clonotypes[pt1]


    ## look for pairs_tuples that will give the same clones file line
    if verbose:
        for pt1, clonotypes in pairs_tuple2clonotypes.items():
            for pt2 in pairs_tuple2clonotypes:
                if pt1 < pt2 and pt1[0] == pt2[0]:
                    print('overlap:', len(pt1), len(pt2), pt1, pt2)


    ## now setup new clonotype2tcrs, clonotype2barcodes mappings
    new_clonotype2tcrs = {}
    new_clonotype2barcodes = {}

    for pairs_tuple, clonotypes in pairs_tuple2clonotypes.items():
        c0 = clonotypes[0]
        if len(clonotypes)>1:
            if verbose:
                print('merging:', ' '.join(clonotypes))
        tcrs = {'A':Counter(), 'B':Counter()}
        for (atcr,btcr) in pairs_tuple:
            tcrs['A'][atcr] += pairing_counts[(atcr, btcr)]
            tcrs['B'][btcr] += pairing_counts[(atcr, btcr)]
        if len(pairs_tuple)==2:
            a1, a2 = pairs_tuple[0][0], pairs_tuple[1][0]
            if tcrs['A'][a1] == tcrs['A'][a2]:
                tcrs['A'][a1] += 1 # ensure the desired order
            else:
                assert tcrs['A'][a1] > tcrs['A'][a2]
        assert len(tcrs['A']) in [1,2]
        assert len(tcrs['B']) == 1
        assert c0 not in new_clonotype2tcrs
        new_clonotype2tcrs[c0] = tcrs
        new_clonotype2barcodes[c0] = list(chain( *(clonotype2barcodes[x] for x in clonotypes)))
        # print 'new:', new_clonotype2barcodes[c0]
        # for c in clonotypes:
        #     print 'old:', clonotype2barcodes[c]
        # print len(new_clonotype2barcodes[c0]), sum( len(clonotype2barcodes[x]) for x in clonotypes)
        assert len(new_clonotype2barcodes[c0]) == sum( len(clonotype2barcodes[x]) for x in clonotypes)


    print('ab_counts:', ab_counts.most_common())
    old_good_clonotypes = [ x for x,y in clonotype2tcrs.items() if len(y['A']) >= 1 and len(y['B']) >= 1 ]
    old_num_barcodes = sum( len(clonotype2barcodes[x]) for x in old_good_clonotypes )
    new_num_barcodes = sum( len(x) for x in list(new_clonotype2barcodes.values()) )

    print('old_num_barcodes:', old_num_barcodes, 'new_num_barcodes:', new_num_barcodes)

    return new_clonotype2tcrs, new_clonotype2barcodes


def make_10x_clones_file(
        filtered_contig_annotations_csvfile,
        organism,
        clones_file, # the OUTPUT file, the one we're making
        stringent = True, # dont believe the 10x clonotypes; reduce 'duplicated' and 'fake' clones
        consensus_annotations_csvfile = None,
):

    clonotype2tcrs, clonotype2barcodes = read_tcr_data( organism, filtered_contig_annotations_csvfile,
                                                        consensus_annotations_csvfile )

    if stringent:
        clonotype2tcrs, clonotype2barcodes = setup_filtered_clonotype_dicts( clonotype2tcrs, clonotype2barcodes )


    _make_clones_file( organism, clones_file, clonotype2tcrs, clonotype2barcodes )

def make_10x_clones_file_batch(
        metadata_file,
        organism,
        clones_file, # the OUTPUT file, the one we're making
        stringent = True, # dont believe the 10x clonotypes; reduce 'duplicated' and 'fake' clones
):

    clonotype2tcrs, clonotype2barcodes = read_tcr_data_batch( organism, metadata_file )

    if stringent:
        clonotype2tcrs, clonotype2barcodes = setup_filtered_clonotype_dicts( clonotype2tcrs, clonotype2barcodes )


    _make_clones_file( organism, clones_file, clonotype2tcrs, clonotype2barcodes )
