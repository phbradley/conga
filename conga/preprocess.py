import scanpy as sc
import pandas as pd
from os.path import exists
from sklearn.metrics import pairwise_distances
from sklearn.utils import sparsefuncs
import numpy as np
import scipy
from anndata import AnnData
from sys import exit
from . import tcr_scores


def get_feature_types_varname( adata ):
    ''' SPECIAL HACK for AGBT dataset
    Figure out the correct "feature_types" varname, e.g. feature_types-0 or feature_types-0-0
    '''
    for name in adata.var: # the columns of the var dataframe
        if name.startswith('feature_types'):
            return name
    print('unable to find feature_types varname')
    print(adata.var_names)
    return None

def check_if_raw_matrix_is_logged( adata ):
    return adata.uns.get( 'raw_matrix_is_logged', False )

def set_raw_matrix_is_logged_to_true( adata ):
    adata.uns[ 'raw_matrix_is_logged' ] = True

def normalize_and_log_the_raw_matrix( adata, counts_per_cell_after = 1e4 ):
    '''
    '''
    if check_if_raw_matrix_is_logged( adata ):
        print('normalize_and_log_the_raw_matrix:: matrix is already logged')
        return adata

    print('Normalize and logging matrix...')

    ft_varname = get_feature_types_varname(adata)
    if ft_varname:
        ngenes = sum( adata.raw.var[ft_varname] != 'Antibody Capture' )
    else:
        ngenes = adata.raw.shape[1]

    X_gex = adata.raw.X[:,:ngenes]
    X_ab  = adata.raw.X[:,ngenes:]

    counts_per_cell = np.sum( X_gex, axis=1 ).A1 # A1 since X_gex is sparse
    assert np.min( counts_per_cell ) > 0
    if np.median( counts_per_cell ) < 100:
        print('WARNING normalize_and_log_the_raw_matrix: low median counts_per_cell.', np.median(counts_per_cell),'\n',
              'has the matrix already been log1p-ed???')
        exit()

    counts_per_cell /= counts_per_cell_after

    sparsefuncs.inplace_row_scale(X_gex, 1/counts_per_cell)

    new_counts_per_cell = np.sum( X_gex, axis=1 ).A1 # A1 since X_gex is sparse
    assert min(new_counts_per_cell) > counts_per_cell_after-1 and max(new_counts_per_cell) < counts_per_cell_after+1

    new_X = scipy.sparse.hstack( [X_gex, X_ab], format="csr" )
    np.log1p( new_X.data, out = new_X.data )

    adata_new = AnnData( X = new_X, obs = adata.obs, var = adata.raw.var )

    adata.raw = adata_new

    set_raw_matrix_is_logged_to_true( adata )

    print(adata)
    return adata

tcr_keys = 'va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq'.split() # ugly

def store_tcrs_in_adata(adata, tcrs):
    ''' returns NOTHING, modifies adata
    tcrs is a list of (atcr,btcr) tuples, where atcr = (va,ja,cdr3a,cdr3a_nucseq) ...
    '''
    assert len(tcrs) == adata.shape[0]

    global tcr_keys
    tcr_indices = ((x,y) for x in range(2) for y in range(4)) # some better itertools soln...

    for tag, (i,j) in zip( tcr_keys, tcr_indices):
        adata.obs[tag] = [x[i][j] for x in tcrs ]
    return


def retrieve_tcrs_from_adata(adata):
    global tcr_keys
    arrays = [ adata.obs[x] for x in tcr_keys ]
    tcrs = []
    for va, ja, cdr3a, cdr3a_nucseq, vb, jb, cdr3b, cdr3b_nucseq in zip( *arrays):
        tcrs.append( ( ( va, ja, cdr3a, cdr3a_nucseq), (vb, jb, cdr3b, cdr3b_nucseq) ) )

    return tcrs

def setup_X_igex( adata ):
    ''' Side effect: will log1p-and-normalize the raw matrix if that's not already done
    '''
    # tmp hacking
    all_genes_file = '/home/pbradley/gitrepos/conga/conga/data/igenes_all_v1.txt'
    kir_genes = 'KIR2DL1 KIR2DL3 KIR2DL4 KIR3DL1 KIR3DL2 KIR3DL3 KIR3DX1'.split()
    extra_genes = [ 'CD4','CD8A','CD8B','CCR7', 'SELL','CCL5','KLRF1','KLRG1','IKZF2','PDCD1','GNG4','MME',
                    'ZNF683','KLRD1','NCR3','KIR3DL1','NCAM1','ITGAM','KLRC2','KLRC3', #'MME',
                    'GNLY','IFNG', 'GZMH','GZMA','CCR6', 'TRAV1-2','KLRB1','ZBTB16', 'GATA3',
                    'TRBC1','EPHB6','SLC4A10','DAD1','ITGB1', 'KLRC1','CD45RA_TotalSeqC','CD45RO_TotalSeqC',
                    'ZNF683','KLRD1','NCR3','KIR3DL1','NCAM1','ITGAM','KLRC2','KLRC3','GNLY' ]
    all_genes = sorted( set( [x[:-1] for x in open(all_genes_file,'r')] + extra_genes + kir_genes ) )
    normalize_and_log_the_raw_matrix(adata) # just in case
    var_names = list( adata.raw.var_names )
    good_genes = [ x for x in all_genes if x in var_names ]
    indices = [ var_names.index(x) for x in good_genes ]
    #print('igex_genes: {}'.format(' '.join(good_genes)))
    X_igex = adata.raw[:,indices].X.toarray()
    return X_igex, good_genes

def read_dataset(
        gex_data,
        gex_data_type,
        clones_file,
        make_var_names_unique = True,
        keep_cells_without_tcrs = False
):
    ''' returns adata

    stores the tcr-dist kPCA info in adata.obsm under the key 'X_pca_tcr'
    stores the tcr info in adata.obs under multiple keys (see store_tcrs_in_adata(...) function)
    '''

    include_tcr_nucseq = True

    print('reading:', gex_data, 'of type', gex_data_type)
    if gex_data_type == 'h5ad':
        adata = sc.read_h5ad( gex_data )

    elif gex_data_type == '10x_mtx':
        adata = sc.read_10x_mtx( gex_data )

    elif gex_data_type == '10x_h5':
        adata = sc.read_10x_h5( gex_data, gex_only=True )

    else:
        print('unrecognized gex_data_type:', gex_data_type, "should be one of ['h5ad', '10x_mtx', '10x_h5']")
        exit()

    if adata.isview: # this is so weird
        adata = adata.copy()

    assert not adata.isview

    if make_var_names_unique:
        adata.var_names_make_unique() # added

    barcodes = set( adata.obs.index )
    print('total barcodes:',len(barcodes),adata.shape)

    # read the kpcs, etc
    bcmap_file = clones_file+'.barcode_mapping.tsv'
    kpcs_file = clones_file[:-4]+'_AB.dist_50_kpcs'
    assert exists(clones_file)
    assert exists(bcmap_file)
    assert exists(kpcs_file)

    print('reading:',clones_file)
    tmpdata = open(clones_file,'r')
    clones_file_header = tmpdata.readline()
    tmpdata.close()
    #clones_file_header = popen('head -n1 '+clones_file).readlines()[0]
    lines = pd.read_csv( clones_file,sep='\t')
    id2tcr = {}
    tcr2id = {}
    tcr2clones_line = {}
    for l in lines.itertuples(index=False):
        if include_tcr_nucseq:
            atcr = ( l.va_gene, l.ja_gene, l.cdr3a, l.cdr3a_nucseq )
            btcr = ( l.vb_gene, l.jb_gene, l.cdr3b, l.cdr3b_nucseq )
        else:
            atcr = ( l.va_gene, l.ja_gene, l.cdr3a )#, l.cdr3a_nucseq )
            btcr = ( l.vb_gene, l.jb_gene, l.cdr3b )#, l.cdr3b_nucseq )
        tcr = ( atcr, btcr )
        #tcr = ( atcr, btcr, l.clone_id ) # hack to guarantee uniqueness!!!
        id2tcr[ l.clone_id ] = tcr
        tcr2id[ tcr ] = l.clone_id
        tcr2clones_line[ tcr ] = '\t'.join(str(x) for x in l)+'\n'


    print('reading:',kpcs_file)
    tcr2kpcs = {}
    lencheck= None
    for line in open(kpcs_file,'r'):
        l = line.split()
        assert l[0] == 'pc_comps:'
        id = l[1]
        kpcs = [float(x) for x in l[2:] if x!='nan' ]
        if lencheck is None:
            lencheck = len(kpcs)
        else:
            assert lencheck == len(kpcs)
        tcr2kpcs[ id2tcr[id] ] = kpcs


    # read the barcode/clonotype mapping info
    barcode2kpcs = {}
    barcode2tcr = {}

    for line in open(bcmap_file,'r'):
        l = line[:-1].split('\t')
        if l[0] == 'clone_id': continue # header line
        if not l[1]: continue
        barcodes = l[1].split(',')
        clone_id = l[0]
        if clone_id not in id2tcr: continue # maybe short cdr3?
        tcr = id2tcr[clone_id]
        kpcs = tcr2kpcs[tcr]
        for bc in barcodes:
            barcode2tcr[ bc ] = tcr
            barcode2kpcs[ bc ] = kpcs
            assert bc in barcodes # the barcodes list before preprocessing...

    mask = [ x in barcode2tcr for x in adata.obs.index ]

    print('slicing:',adata.shape, len(mask), np.sum(mask))
    assert not adata.isview
    #adata = adata[mask,:]
    adata = adata[mask,:].copy()
    assert not adata.isview
    print('sliced:',adata.shape)
    #stdout.flush()


    # stash the kPCA info in adata.obsm
    X_kpca = np.array( [ barcode2kpcs[x] for x in adata.obs.index ] )
    adata.obsm['X_pca_tcr'] = X_kpca

    tcrs = [ barcode2tcr[x] for x in adata.obs.index ]
    store_tcrs_in_adata( adata, tcrs )

    adata.obs['mhci_score']      = [ tcr_scores.mhci_score_tcr(x) for x in tcrs ]
    adata.obs['cd8_score']       = [ tcr_scores.cd8_score_tcr(x) for x in tcrs ]
    adata.obs['alphadist_score'] = [ tcr_scores.alphadist_score_tcr(x) for x in tcrs ]
    adata.obs['cdr3len_score']   = [ tcr_scores.cdr3len_score_tcr(x) for x in tcrs ]

    return adata



#adapt if want to compare more than one data frame
#return or pass in adata?
def filter_normalize_and_hvg(
        adata,
        min_genes=200,
        min_cells=3,
        n_genes=2000,
        percent_mito=0.1,
        hvg_min_mean = 0.0125,
        hvg_max_mean=3,
        antibody = False,
        hvg_min_disp=0.5,
        exclude_TR_genes = True,
):
    '''Filters cells and genes to find highly variable genes'''

    ## notes:
    ## * AGTB specific things that we removed
    ##   - exclude "antibody" features before normalizing
    ##   - "qcheck" thingy


    #filters out cells with less than given number of genes expressed (200)
    #filters out genes present in less than given number of cells (3)
    #adds n_genes to obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)


    #find mitochondrial genes
    #add percent_mito and n_counts to obs
    mito_genesA = ( adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-') )
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genesA].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1

    #filter n_genes and percent_mito based on param
    print('filtered out {} cells with more than {} genes'\
          .format( np.sum( adata.obs['n_genes'] >= n_genes ), n_genes ) )
    adata = adata[adata.obs['n_genes'] < n_genes, :]
    print('filtered out {} cells with more than {} percent mito'\
          .format( np.sum( adata.obs['percent_mito'] >= percent_mito ), percent_mito ) )
    adata = adata[adata.obs['percent_mito'] < percent_mito, :]

    adata.raw = adata

    feature_types_colname = get_feature_types_varname( adata )
    if feature_types_colname:
        mask = (adata.var[feature_types_colname]=='Antibody Capture' )
        print('num antibody features:', np.sum(mask))
        if np.sum(mask):
            assert not antibody # sanity check
            oldshape = adata.shape
            oldrawshape = adata.raw.shape
            adata = adata[:,~mask]
            newshape = adata.shape
            newrawshape = adata.raw.shape
            removed_at_end = ( np.sum( mask[:newshape[1]] )==0 )
            print('Removed {} antibody features from adata, using colname {}'\
                  .format( np.sum(mask), feature_types_colname ), oldshape, oldrawshape, newshape, newrawshape,
                  removed_at_end )
            assert removed_at_end # want to make this assumption somewhere else

    #normalize and log data
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

    #find and filter by highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=hvg_min_mean, max_mean=hvg_max_mean, min_disp=hvg_min_disp)

    # exclude TCR genes
    if exclude_TR_genes:
        is_tr = []
        for gene in adata.var.index:
            is_tr.append( gene.lower().startswith('trav') or gene.lower().startswith('trbv') or \
                            gene.lower().startswith('traj') or gene.lower().startswith('trbj') )
        print('excluding {} TR genes'.format(sum(is_tr)))
        adata.var['highly_variable'][ is_tr ] = False
        assert sum( adata.var['highly_variable'][ is_tr ] )==0 # sanity check


    print('total of',np.sum(adata.var['highly_variable']),'variable genes',adata.shape)
    oldshape = adata.shape
    oldrawshape = adata.raw.shape
    adata = adata[:, adata.var['highly_variable']]
    newshape = adata.shape
    newrawshape = adata.raw.shape
    print('after slice:', len(adata.var_names), oldshape, oldrawshape, newshape, newrawshape )

    return adata


def cluster_and_tsne_and_umap(adata, louvain_resolution= None, compute_pca_gex= True):
    '''calculates neighbors, tsne, louvain for both GEX and TCR

    stores in  adata:

    obsm: X_pca_gex
    obsm: X_tsne_gex
    obsm: X_tsne_tcr
    obsm: X_gex_2d (same as) X_umap_gex
    obsm: X_tcr_2d (same as) X_umap_tcr
    obs: clusters_gex (int version of) louvain_gex
    obs: clusters_tcr (int version of) louvain_tcr

    '''
    ncells = adata.shape[0]
    n_components = min(ncells, 50)

    if compute_pca_gex:
        sc.tl.pca(adata, n_comps=n_components) # re-run now that we have reduced to a single cell per clone
        adata.obsm['X_pca_gex'] = adata.obsm['X_pca']
    assert 'X_pca_gex' in adata.obsm_keys()
    assert 'X_pca_tcr' in adata.obsm_keys()


    for tag in ['gex','tcr']:
        adata.obsm['X_pca'] = adata.obsm['X_pca_'+tag]

        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40, n_components))
        sc.tl.tsne(adata)
        adata.obsm['X_tsne_'+tag] = adata.obsm['X_tsne']
        sc.tl.umap(adata)
        adata.obsm['X_umap_'+tag] = adata.obsm['X_umap']
        sc.tl.louvain(adata, key_added='louvain_'+tag, resolution=louvain_resolution)


    ## set better obsm keys and data types
    adata.obs['clusters_gex'] = np.copy(adata.obs['louvain_gex']).astype(int)
    adata.obs['clusters_tcr'] = np.copy(adata.obs['louvain_tcr']).astype(int)

    adata.obsm['X_gex_2d'] = adata.obsm['X_umap_gex']
    adata.obsm['X_tcr_2d'] = adata.obsm['X_umap_tcr']


    return adata

def filter_and_scale( adata ):

    ## now process as before
    adata = filter_normalize_and_hvg(adata, exclude_TR_genes= True, percent_mito=0.1)

    ## should consider adding cell cycle here:
    sc.pp.regress_out(adata, ['n_counts','percent_mito'])

    sc.pp.scale(adata, max_value=10)

    return adata


def reduce_to_single_cell_per_clone(
        adata,
        n_pcs=50
):
    ''' returns adata

    needs the X array to be ready for pca

    the raw array is normalized and logged if not already done

    stashes info in adata:
    obs: clone_sizes
    obsm: X_igex
    obsm: X_pca  (by sc.tl.pca() using the X array; but note that this gets recomputed in cluster_and_tsne_and_umap())
    uns: X_igex_genes
    '''

    # compute pcs
    print('compute pca to find rep cell for each clone', adata.shape)
    sc.tl.pca(adata, n_comps=min(adata.shape[0], n_pcs))

    # for each clone, try to identify the most 'representative' cell based on gex
    tcrs_with_duplicates = retrieve_tcrs_from_adata(adata)
    tcrs = sorted( set( tcrs_with_duplicates))

    num_clones = len(tcrs)
    print('num_clones:',num_clones)

    tcr2clone_id = { y:x for x,y in enumerate(tcrs) }

    #
    # get the interesting genes matrix before we subset
    # this will normalize and log the raw array if not done
    X_igex, good_genes = setup_X_igex(adata)

    clone_ids = np.array( [ tcr2clone_id[x] for x in tcrs_with_duplicates ] )

    X_pca = adata.obsm['X_pca']
    clone_sizes = []
    rep_cell_indices = [] # parallel
    new_X_igex = []

    ## for each clone (tcr) we pick a single representative cell, stored in rep_cell_indices
    ## rep_cell_indices is parallel with and aligned to the tcrs list
    for c in range(num_clones):
        clone_cells = np.array( [ x for x,y in enumerate(clone_ids) if y==c] )
        clone_size = clone_cells.shape[0]
        clone_sizes.append( clone_size )
        if clone_size==1:
            rep_cell_index = clone_cells[0]
        else:
            X_pca_clone = X_pca[ clone_cells, : ]
            #print('X_pca_clone:',X_pca_clone.shape)
            assert X_pca_clone.shape[0] == clone_size
            D_gex_clone = pairwise_distances( X_pca_clone, X_pca_clone )
            assert D_gex_clone.shape  == ( clone_size,clone_size )
            avgdist = D_gex_clone.sum(axis=1)
            #print('min_avgdist', np.min(avgdist)/clone_size, clone_size, np.argmin(avgdist))
            rep_cell_index = clone_cells[ np.argmin( avgdist ) ]
        rep_cell_indices.append(rep_cell_index)
        new_X_igex.append( np.sum( X_igex[clone_cells,:], axis=0 ) / clone_size )

    rep_cell_indices = np.array(rep_cell_indices)
    new_X_igex = np.array(new_X_igex)
    assert new_X_igex.shape == ( num_clones, len(good_genes))


    print('before subset barcodes:', adata.shape)
    adata = adata[ rep_cell_indices, : ].copy() ## seems like we need to copy here, something to do with adata 'views'
    print('after subset barcodes:', adata.shape)
    adata.obs['clone_sizes'] = np.array( clone_sizes )
    adata.obsm['X_igex'] = new_X_igex
    adata.uns['X_igex_genes'] = good_genes

    tcrs_redo = retrieve_tcrs_from_adata(adata)
    assert tcrs_redo == tcrs # sanity check

    adata = normalize_and_log_the_raw_matrix( adata ) # prob not necessary

    return adata


