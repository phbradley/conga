import scanpy as sc
import random
import pandas as pd
from os.path import exists
from pathlib import PurePath
from collections import Counter
from sklearn.metrics import pairwise_distances
from sklearn.utils import sparsefuncs
from sklearn.decomposition import KernelPCA
import numpy as np
import scipy
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform, cdist
from scipy.sparse import issparse#, csr_matrix
from anndata import AnnData
import sys
import os
from sys import exit
from . import tcr_scoring
from . import util
from . import pmhc_scoring
from . import plotting
from .tcrdist.tcr_distances import TcrDistCalculator

# silly hack
all_sexlinked_genes = frozenset('XIST DDX3Y EIF1AY KDM5D LINC00278 NLGN4Y RPS4Y1 TTTY14 TTTY15 USP9Y UTY ZFY'.split())


def check_if_raw_matrix_is_logged( adata ):
    return adata.uns.get( 'raw_matrix_is_logged', False )

def set_raw_matrix_is_logged_to_true( adata ):
    adata.uns[ 'raw_matrix_is_logged' ] = True


def add_mait_info_to_adata_obs( adata, key_added = 'is_mait' ):
    ''' Note that for mouse we are actually doing iNKT cells since those are the common ones
    Also this is a total approximation based on the TCR alpha chain
    '''
    if key_added in adata.obs:
        #print('adata.obs already has mait info')
        return
    tcrs = retrieve_tcrs_from_adata(adata)
    organism = 'human' if 'organism' not in adata.uns_keys() else adata.uns['organism']
    if 'human' in organism:
        is_mait = [ tcr_scoring.is_human_mait_alpha_chain(x[0]) for x in tcrs ]
    else:
        assert 'mouse' in organism
        is_mait = [ tcr_scoring.is_mouse_inkt_alpha_chain(x[0]) for x in tcrs ]
    adata.obs['is_mait'] = is_mait


def normalize_and_log_the_raw_matrix( adata, counts_per_cell_after = 1e4 ):
    '''
    '''
    if check_if_raw_matrix_is_logged( adata ):
        print('normalize_and_log_the_raw_matrix:: matrix is already logged')
        return adata

    print('Normalize and logging matrix...')

    ft_varname = pmhc_scoring.get_feature_types_varname(adata)
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

    #print(adata)
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
    all_genes_file = PurePath.joinpath( util.path_to_data , 'igenes_all_v1.txt')
    kir_genes = 'KIR2DL1 KIR2DL3 KIR2DL4 KIR3DL1 KIR3DL2 KIR3DL3 KIR3DX1'.split()
    extra_genes = [ 'CD4','CD8A','CD8B','CCR7', 'SELL','CCL5','KLRF1','KLRG1','IKZF2','PDCD1','GNG4','MME',
                    'ZNF683','KLRD1','NCR3','KIR3DL1','NCAM1','ITGAM','KLRC2','KLRC3', #'MME',
                    'GNLY','IFNG', 'GZMH','GZMA','CCR6', 'TRAV1-2','KLRB1','ZBTB16', 'GATA3',
                    'TRBC1','EPHB6','SLC4A10','DAD1','ITGB1', 'KLRC1','CD45RA_TotalSeqC','CD45RO_TotalSeqC',
                    'ZNF683','KLRD1','NCR3','KIR3DL1','NCAM1','ITGAM','KLRC2','KLRC3','GNLY',
                    'CCR7_TotalSeqC','CD3_TotalSeqC','IgG1_TotalSeqC','PD-1_TotalSeqC','CD5' ]

    all_genes = sorted( set( [x[:-1] for x in open(all_genes_file,'r')] + extra_genes + kir_genes ) )

    organism = 'human'
    if 'organism' in adata.uns_keys():
        organism = adata.uns['organism']

    if 'mouse' in organism: # this is a temporary hack -- could actually load a mapping between mouse and human genes
        all_genes = [x.capitalize() for x in all_genes]

    for g in plotting.default_logo_genes[organism] + plotting.default_gex_header_genes[organism]:
        if g not in all_genes:
            all_genes.append(g)

    normalize_and_log_the_raw_matrix(adata) # just in case
    var_names = list( adata.raw.var_names )
    good_genes = [ x for x in all_genes if x in var_names ]
    print('found {} of {} good_genes in var_names  organism={}'.format(len(good_genes), len(all_genes), organism))
    assert good_genes
    indices = [ var_names.index(x) for x in good_genes ]
    X_igex = adata.raw[:,indices].X.toarray()
    return X_igex, good_genes


def read_adata(
        gex_data, # filename
        gex_data_type # string describing file type
):
    ''' Split this out so that other code can use it. Read GEX data
    '''
    print('reading:', gex_data, 'of type', gex_data_type)
    if gex_data_type == 'h5ad':
        adata = sc.read_h5ad( gex_data )

    elif gex_data_type == '10x_mtx':
        adata = sc.read_10x_mtx( gex_data ) # consider-- do we need gex_only here? is it also an option for this fxn?

    elif gex_data_type == '10x_h5':
        adata = sc.read_10x_h5( gex_data, gex_only=True )

    elif gex_data_type == 'loom':
        adata = sc.read_loom( gex_data )

    else:
        print('unrecognized gex_data_type:', gex_data_type, "should be one of ['h5ad', '10x_mtx', '10x_h5', 'loom']")
        exit()

    if adata.isview: # this is so weird
        adata = adata.copy()
    return adata

def read_dataset(
        gex_data,
        gex_data_type,
        clones_file,
        make_var_names_unique = True,
        keep_cells_without_tcrs = False,
        kpca_file = None,
        allow_missing_kpca_file = False,
):
    ''' returns adata

    stores the tcr-dist kPCA info in adata.obsm under the key 'X_pca_tcr'
    stores the tcr info in adata.obs under multiple keys (see store_tcrs_in_adata(...) function)
    '''

    include_tcr_nucseq = True

    adata = read_adata(gex_data, gex_data_type)

    if make_var_names_unique:
        adata.var_names_make_unique() # added

    barcodes = set( adata.obs.index )
    print('total barcodes:',len(barcodes),adata.shape)

    # read the kpcs, etc
    bcmap_file = clones_file+'.barcode_mapping.tsv'
    if kpca_file is None:
        kpca_file = clones_file[:-4]+'_AB.dist_50_kpcs'
    assert exists(clones_file)
    assert exists(bcmap_file)
    if not allow_missing_kpca_file:
        assert exists(kpca_file)

    print('reading:',clones_file)
    tmpdata = open(clones_file,'r')
    clones_file_header = tmpdata.readline()
    tmpdata.close()
    #clones_file_header = popen('head -n1 '+clones_file).readlines()[0]
    lines = pd.read_csv( clones_file,sep='\t')
    id2tcr = {}
    tcr2id = {}
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

    tcr2kpcs = {}
    if not exists(kpca_file):
        if not allow_missing_kpca_file:
            print('ERROR: missing kpca_file:', kpca_file)
            sys.exit(1)
        else:
            missing_kpca_file = True
            print('WARNING: missing kpca_file:', kpca_file)
            print('WARNING: X_tcr_pca will be empty')
    else:
        missing_kpca_file = False
        print('reading:',kpca_file)
        lencheck= None
        for line in open(kpca_file,'r'):
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
        if not missing_kpca_file:
            kpcs = tcr2kpcs[tcr]
        for bc in barcodes:
            barcode2tcr[ bc ] = tcr
            if not missing_kpca_file:
                barcode2kpcs[ bc ] = kpcs
            assert bc in barcodes # the barcodes list before preprocessing...

    mask = [ x in barcode2tcr for x in adata.obs.index ]

    print(f'Reducing to the {np.sum(mask)} barcodes (out of {adata.shape[0]}) with paired TCR sequence data')
    assert not adata.isview
    #adata = adata[mask,:]
    adata = adata[mask,:].copy()
    assert not adata.isview


    if not missing_kpca_file: # stash the kPCA info in adata.obsm
        X_kpca = np.array( [ barcode2kpcs[x] for x in adata.obs.index ] )
        adata.obsm['X_pca_tcr'] = X_kpca

    tcrs = [ barcode2tcr[x] for x in adata.obs.index ]
    store_tcrs_in_adata( adata, tcrs )

    return adata



#adapt if want to compare more than one data frame
#return or pass in adata?
def filter_normalize_and_hvg(
        adata,
        min_genes=None,
        n_genes=None,
        percent_mito=None,
        min_cells=3,
        hvg_min_mean = 0.0125,
        hvg_max_mean=3,
        antibody = False,
        hvg_min_disp=0.5,
        exclude_TR_genes = True,
        also_exclude_TR_constant_region_genes = True, ## this is NEW and not the default for the old TCR analysis!!!!
        exclude_sexlinked = False,
):
    '''Filters cells and genes to find highly variable genes'''

    organism = adata.uns['organism']

    if min_genes is None:
        min_genes=200
        print('min_genes not set. Using default ' + str(min_genes) )

    if n_genes is None:
        n_genes=2000
        print('n_genes not set. Using default ' + str(n_genes) )

    if percent_mito is None:
        percent_mito=0.1
        print('percent_mito not set. Using default ' + str(percent_mito) )

    ## notes:
    ## * AGTB specific things that we removed
    ##   - exclude "antibody" features before normalizing
    ##   - "qcheck" thingy

    global all_sexlinked_genes

    #filters out cells with less than given number of genes expressed (200)
    #filters out genes present in less than given number of cells (3)
    #adds n_genes to obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    assert not adata.isview
    sc.pp.filter_genes(adata, min_cells=min_cells)
    assert not adata.isview


    #find mitochondrial genes
    #add percent_mito and n_counts to obs
    mito_genesA = ( adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-') )
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genesA].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1

    #filter n_genes and percent_mito based on param
    print('filtered out {} cells with more than {} genes'\
          .format( np.sum( adata.obs['n_genes'] >= n_genes ), n_genes ) )
    assert not adata.isview
    adata = adata[adata.obs['n_genes'] < n_genes, :].copy()
    assert not adata.isview
    print('filtered out {} cells with more than {} percent mito'\
          .format( np.sum( adata.obs['percent_mito'] >= percent_mito ), percent_mito ) )
    adata = adata[adata.obs['percent_mito'] < percent_mito, :].copy()
    assert not adata.isview

    adata.raw = adata

    feature_types_colname = pmhc_scoring.get_feature_types_varname( adata )
    if feature_types_colname:
        mask = (adata.var[feature_types_colname]=='Antibody Capture' )
        print('num antibody features:', np.sum(mask))
        if np.sum(mask):
            assert not antibody # sanity check
            oldshape = adata.shape
            oldrawshape = adata.raw.shape
            adata = adata[:,~mask].copy()
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
    hvg_mask = np.array(adata.var['highly_variable'])

    # exclude TCR genes
    if exclude_TR_genes:
        is_tr = []
        for gene in adata.var.index:
            is_tr.append( util.is_vdj_gene(gene, organism,
                                           include_constant_regions=also_exclude_TR_constant_region_genes))
        print('excluding {} TR genes ({} variable)'.format(sum(is_tr), sum(hvg_mask[is_tr])))
        hvg_mask[is_tr] = False
        assert sum( hvg_mask[is_tr] )==0 # sanity check

    if exclude_sexlinked:
        is_sexlinked = [ x in all_sexlinked_genes for x in adata.var.index ]
        print('excluding {} sexlinked genes'.format(sum(is_sexlinked)))
        hvg_mask[is_sexlinked] = False

    print('total of', np.sum(hvg_mask), 'variable genes', adata.shape)
    oldshape = adata.shape
    oldrawshape = adata.raw.shape
    adata = adata[:, hvg_mask].copy()
    newshape = adata.shape
    newrawshape = adata.raw.shape
    #print('after slice:', len(adata.var_names), oldshape, oldrawshape, newshape, newrawshape )

    return adata


def cluster_and_tsne_and_umap(
        adata,
        louvain_resolution= None,
        compute_pca_gex= True,
        skip_tsne=True, # yes this is silly
        skip_tcr=False,
        clustering_method=None
):
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
    assert clustering_method in [None, 'louvain', 'leiden']

    ncells = adata.shape[0]
    n_components = min(ncells-1, 50)

    if compute_pca_gex:
        # switch to arpack for better reproducibility
        # re-run now that we have reduced to a single cell per clone
        sc.tl.pca(adata, svd_solver='arpack', n_comps=n_components)
        adata.obsm['X_pca_gex'] = adata.obsm['X_pca']
    assert 'X_pca_gex' in adata.obsm_keys()
    if not skip_tcr:
        assert 'X_pca_tcr' in adata.obsm_keys()


    for tag in ['gex','tcr']:
        if skip_tcr and tag=='tcr': # hack for analyzing GEX before reducing to clones
            continue

        adata.obsm['X_pca'] = adata.obsm['X_pca_'+tag]
        if tag == 'gex': # tmp hack...
            # raw to louvain uses 50 for pca, but only 40 for neighbors
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40,n_components)) #
        else:
            # analyze_sorted_generic uses 50
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_components) # had a 40 in there...
        if not skip_tsne:
            sc.tl.tsne(adata, n_pcs=n_components)
            adata.obsm['X_tsne_'+tag] = adata.obsm['X_tsne']
        sc.tl.umap(adata)
        adata.obsm['X_umap_'+tag] = adata.obsm['X_umap']
        resolution = 1.0 if louvain_resolution is None else louvain_resolution
        if clustering_method=='louvain':
            cluster_key_added = 'louvain_'+tag
            sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
            print('ran louvain clustering:', cluster_key_added)
        elif clustering_method=='leiden':
            cluster_key_added = 'leiden_'+tag
            sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key_added)
            print('ran leiden clustering:', cluster_key_added)
        else: # try both, first louvain, then leiden
            try:
                cluster_key_added = 'louvain_'+tag
                sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
                print('ran louvain clustering:', cluster_key_added)
            except ImportError: # hacky
                cluster_key_added = 'leiden_'+tag
                sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key_added)
                print('ran leiden clustering:', cluster_key_added)

        ## set better obsm keys and data types
        # generic names-- these will be the ones we use in other stuff
        adata.obs['clusters_'+tag] = np.copy(adata.obs[cluster_key_added]).astype(int)
        adata.obsm['X_{}_2d'.format(tag)] = adata.obsm['X_umap_'+tag]


    del adata.obsm['X_umap']
    del adata.obsm['X_pca']

    return adata

def filter_and_scale(
        adata,
        min_genes = None,
        n_genes= None,
        percent_mito= None
):
    ## now process as before
    adata = filter_normalize_and_hvg(
        adata, min_genes=min_genes, n_genes=n_genes, percent_mito=percent_mito,
        exclude_TR_genes= True, exclude_sexlinked=True)

    ## should consider adding cell cycle here:
    sc.pp.regress_out(adata, ['n_counts','percent_mito'])

    sc.pp.scale(adata, max_value=10)

    #stash as a layer for gex analysis

    adata.layers['scaled'] = sc.pp.scale(adata, max_value=10, copy=True).X

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
    uns: X_igex_genes
    '''

    # compute pcs
    print('compute pca to find rep cell for each clone', adata.shape)
    # switch to arpack for better reproducibility
    sc.tl.pca(adata, svd_solver='arpack', n_comps=min(adata.shape[0]-1, n_pcs))

    # for each clone, try to identify the most 'representative' cell based on gex
    tcrs_with_duplicates = retrieve_tcrs_from_adata(adata)
    tcrs = sorted( set( tcrs_with_duplicates))

    num_clones = len(tcrs)
    print('num_clones:',num_clones)

    tcr2clone_id = { y:x for x,y in enumerate(tcrs) }

    #
    # get the interesting genes matrix before we subset
    # this will normalize and log the raw array if not done!!
    X_igex, good_genes = setup_X_igex(adata)


    if 'pmhc_var_names' in adata.uns_keys():
        pmhc_var_names = adata.uns['pmhc_var_names']
        X_pmhc = pmhc_scoring._get_X_pmhc(adata, pmhc_var_names)
        new_X_pmhc = []
    else:
        pmhc_var_names = None

    clone_ids = np.array( [ tcr2clone_id[x] for x in tcrs_with_duplicates ] )

    X_pca = adata.obsm['X_pca']
    clone_sizes = []
    rep_cell_indices = [] # parallel
    new_X_igex = []

    ## for each clone (tcr) we pick a single representative cell, stored in rep_cell_indices
    ## rep_cell_indices is parallel with and aligned to the tcrs list
    for c in range(num_clones):
        if c%1000==0:
            print('choose representative cell for clone:', c, num_clones, adata.shape)
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
        if pmhc_var_names:
            new_X_pmhc.append( np.sum( X_pmhc[clone_cells,:], axis=0 ) / clone_size )

    rep_cell_indices = np.array(rep_cell_indices)
    new_X_igex = np.array(new_X_igex)
    assert new_X_igex.shape == ( num_clones, len(good_genes))

    print(f'reduce from {adata.shape[0]} cells to {len(rep_cell_indices)} cells (one per clonotype)')
    adata = adata[ rep_cell_indices, : ].copy() ## seems like we need to copy here, something to do with adata 'views'
    adata.obs['clone_sizes'] = np.array( clone_sizes )
    adata.obsm['X_igex'] = new_X_igex
    adata.uns['X_igex_genes'] = good_genes
    if pmhc_var_names:
        new_X_pmhc = np.array(new_X_pmhc)
        assert new_X_pmhc.shape == ( num_clones, len(pmhc_var_names))
        adata.obsm['X_pmhc'] = new_X_pmhc

    tcrs_redo = retrieve_tcrs_from_adata(adata)
    assert tcrs_redo == tcrs # sanity check

    adata = normalize_and_log_the_raw_matrix( adata ) # prob not necessary; already done for X_igex ?

    del adata.obsm['X_pca'] # delete the old pcas
    return adata


def write_proj_info( adata, outfile ):
    clusters_gex = np.array(adata.obs['clusters_gex'])
    clusters_tcr = np.array(adata.obs['clusters_tcr'])
    X_gex_2d = adata.obsm['X_gex_2d']
    X_tcr_2d = adata.obsm['X_tcr_2d']
    tcrs = retrieve_tcrs_from_adata(adata)

    out = open(outfile,'w')
    for ii in range(adata.shape[0]):
        out.write('ii: {} X_gex_2d: {:.3f} {:.3f} X_tcr_2d: {:.3f} {:.3f} cluster_gex: {} cluster_tcr: {} tcr: {} {}'\
                  .format(ii, X_gex_2d[ii,0], X_gex_2d[ii,1], X_tcr_2d[ii,0], X_tcr_2d[ii,1],
                          clusters_gex[ii], clusters_tcr[ii],
                          ' '.join(tcrs[ii][0]), ' '.join(tcrs[ii][1])))
    out.close()



def filter_cells_by_ribo_norm(adata):
    ''' returns  new filtered adata
    will normalize_and_log_the_raw_matrix if not already done
    '''
    from scipy.stats import gaussian_kde

    normalize_and_log_the_raw_matrix(adata)

    X_norm = adata.raw.X

    X_norm_ribo = []
    var_names = [x.upper() for x in adata.raw.var_names ] ## NOTE UPPER-- help for mouse
    for ig,g in enumerate(var_names):
        if g.startswith('RP') and len(g)>=4 and g[2] in 'SL' and g[3].isdigit():
            X_norm_ribo.append( X_norm[:,ig].toarray() )

    X_norm_ribo = np.sum(np.hstack( X_norm_ribo ), axis=1 )/len(X_norm_ribo)

    dens = gaussian_kde(X_norm_ribo)
    # look for two maxima and a minima in between
    # or look for the closest minimum to 2.0
    step = 0.01
    x=2.0
    while True:
        y0 = dens(x)[0]
        yl = dens(x-step)[0]
        yr = dens(x+step)[0]
        print('filter_cells_by_ribo_norm:: minfind:',x,y0,yl,yr)
        if yl < yr:
            if yl < y0:
                x = x-step
            else:
                # found a local min
                ribo_norm_threshold=x
                break
        else:
            if yr < y0:
                x = x+step
            else:
                ribo_norm_threshold=x
                break
    print('filter_cells_by_ribo_norm:: ribo_norm localmin:',ribo_norm_threshold)

    mask = (X_norm_ribo>ribo_norm_threshold)
    print('filter_cells_by_ribo_norm::', adata.shape, 'downto:', np.sum(mask))
    sys.stdout.flush()
    adata = adata[mask,:].copy()
    assert np.sum(mask) == adata.shape[0]
    return adata


def setup_tcr_groups( adata ):
    tcrs = retrieve_tcrs_from_adata(adata)

    atcrs = sorted( set( x[0] for x in tcrs ) )
    btcrs = sorted( set( x[1] for x in tcrs ) )

    atcr2agroup = dict( (y,x) for x,y in enumerate(atcrs))
    btcr2bgroup = dict( (y,x) for x,y in enumerate(btcrs))

    agroups = np.array( [ atcr2agroup[x[0]] for x in tcrs] )
    bgroups = np.array( [ btcr2bgroup[x[1]] for x in tcrs] )

    return agroups, bgroups


def _calc_nndists( D, nbrs ):
    batch_size, num_nbrs = nbrs.shape
    assert D.shape[0] == batch_size
    sample_range = np.arange(batch_size)[:, np.newaxis]
    nbrs_sorted = nbrs[sample_range, np.argsort(D[sample_range, nbrs])]
    D_nbrs_sorted = D[sample_range, nbrs_sorted]
    wts = np.linspace(1.0, 1.0/num_nbrs, num_nbrs)
    wts /= np.sum(wts)
    nndists = np.sum( D_nbrs_sorted * wts[np.newaxis,:], axis=1)
    assert nndists.shape==(batch_size,)
    return nndists

def _calc_nndists_old( D, nbrs ):
    num_clones, num_nbrs = nbrs.shape
    assert D.shape ==(num_clones, num_clones)
    sample_range = np.arange(num_clones)[:, np.newaxis]
    nbrs_sorted = nbrs[sample_range, np.argsort(D[sample_range, nbrs])]
    D_nbrs_sorted = D[sample_range, nbrs_sorted]
    wts = np.linspace(1.0, 1.0/num_nbrs, num_nbrs)
    wts /= np.sum(wts)
    nndists = np.sum( D_nbrs_sorted * wts[np.newaxis,:], axis=1)
    assert nndists.shape==(num_clones,)
    return nndists



def calc_nbrs_batched(
        adata,
        nbr_fracs,
        obsm_tag_gex = 'X_pca_gex',
        obsm_tag_tcr = 'X_pca_tcr', # set to None to skip tcr calc
        also_calc_nndists = False,
        nbr_frac_for_nndists = None,
        target_N_for_batching = 8192,
):
    ''' returns dict mapping from nbr_frac to [nbrs_gex, nbrs_tcr]

    nbrs exclude self and any clones in same atcr group or btcr group
    '''
    if also_calc_nndists:
        assert nbr_frac_for_nndists in nbr_fracs

    N = adata.shape[0]
    # try for a dataset size of target_N, ie batch_size*N = target_N**2, batch_size = target_N**2 / N
    batch_size = max(10, int(target_N_for_batching**2/N))
    num_batches = (N-1)//batch_size + 1

    agroups, bgroups = setup_tcr_groups(adata)

    all_nbrs = {}
    for nbr_frac in nbr_fracs:
        all_nbrs[nbr_frac] = [ [], [] ]

    nndists = [ [], [] ]

    for bb in range(num_batches):
        b_start = bb*batch_size
        b_stop = min(N, (bb+1)*batch_size)
        batch_indices = np.arange(b_start, b_stop)

        for itag, (tag, obsm_tag) in enumerate( [['gex', obsm_tag_gex], ['tcr', obsm_tag_tcr]] ):
            if obsm_tag is None:
                print('skipping nbr calc for', tag, obsm_tag)
                continue
            print(f'compute D {tag} batch= {bb} num_batches= {num_batches} N= {N} batch_size= {batch_size}')
            X = adata.obsm[obsm_tag]
            D = cdist( X[b_start:b_stop, :], X )

            for ii in batch_indices:
                D[ii-b_start, (agroups==agroups[ii]) ] = 1e3
                D[ii-b_start, (bgroups==bgroups[ii]) ] = 1e3

            for nbr_frac in nbr_fracs:
                num_neighbors = max(1, int(nbr_frac*adata.shape[0]))
                print(f'argpartitions: {tag} batch= {bb} nbr_frac= {nbr_frac}')
                nbrs = np.argpartition( D, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
                assert nbrs.shape == (D.shape[0], num_neighbors)
                all_nbrs[nbr_frac][itag].append(nbrs)

                if also_calc_nndists and nbr_frac == nbr_frac_for_nndists:
                    print('calculate nndists:', tag, nbr_frac)
                    nndists[itag].append( _calc_nndists( D, nbrs ))
                    print('DONE calculating nndists:', nbr_frac)

    for nbr_frac in nbr_fracs:
        nbrslist_gex, nbrslist_tcr = all_nbrs[nbr_frac]
        all_nbrs[nbr_frac] = [ np.vstack(nbrslist_gex), None if obsm_tag_tcr is None else np.vstack(nbrslist_tcr) ]

    if also_calc_nndists:
        nndists_gex = np.hstack(nndists[0])
        nndists_tcr = None if obsm_tag_tcr is None else np.hstack(nndists[1])

        return all_nbrs, nndists_gex, nndists_tcr
    else:
        return all_nbrs


def calc_nbrs(
        adata,
        nbr_fracs,
        obsm_tag_gex = 'X_pca_gex',
        obsm_tag_tcr = 'X_pca_tcr', # set to None to skip calc
        also_calc_nndists = False,
        nbr_frac_for_nndists = None,
        target_N_for_batching = 8192,
):
    ''' returns dict mapping from nbr_frac to [nbrs_gex, nbrs_tcr]

    nbrs exclude self and any clones in same atcr group or btcr group
    '''
    if adata.shape[0] > 1.25*target_N_for_batching: ## EARLY RETURN
        return calc_nbrs_batched(adata, nbr_fracs, obsm_tag_gex, obsm_tag_tcr, also_calc_nndists, nbr_frac_for_nndists,
                                 target_N_for_batching)

    if also_calc_nndists:
        assert nbr_frac_for_nndists in nbr_fracs

    all_nbrs = {}
    for nbr_frac in nbr_fracs:
        all_nbrs[nbr_frac] = [ None, None ]
    nndists = [ None, None ]

    agroups, bgroups = setup_tcr_groups(adata)
    for itag, (tag, obsm_tag) in enumerate([['gex', obsm_tag_gex], ['tcr', obsm_tag_tcr]]):
        if obsm_tag is None:
            print('skipping', tag, 'nbr calc:', obsm_tag)
            continue

        print('compute D', tag, adata.shape[0])
        D = pairwise_distances( adata.obsm[obsm_tag], metric='euclidean' )
        for ii,(a,b) in enumerate(zip(agroups, bgroups)):
            D[ii, (agroups==a) ] = 1e3
            D[ii, (bgroups==b) ] = 1e3

        for nbr_frac in nbr_fracs:
            num_neighbors = max(1, int(nbr_frac*adata.shape[0]))
            print('argpartitions:', nbr_frac, adata.shape[0], tag)
            nbrs = np.argpartition( D, num_neighbors-1 )[:,:num_neighbors] # will NOT include self in there
            assert nbrs.shape == (adata.shape[0], num_neighbors)
            all_nbrs[nbr_frac][itag] = nbrs

            if also_calc_nndists and nbr_frac == nbr_frac_for_nndists:
                print('calculate nndists:', tag, nbr_frac)
                nndists[itag] = _calc_nndists( D, nbrs )
                print('DONE calculating nndists:', tag, nbr_frac)

    if also_calc_nndists:
        return all_nbrs, nndists[0], nndists[1]
    else:
        return all_nbrs


def recalculate_tcrdist_nbrs(
        adata,
        nbr_fracs,
        nbr_frac_for_nndists = None, # if not None, calculate nndists at this nbr fraction
):
    ''' returns dict mapping from nbr_frac to nbrs_tcr

    nbrs exclude self and any clones in same atcr group or btcr group
    '''
    agroups, bgroups = setup_tcr_groups(adata)

    tcrs = retrieve_tcrs_from_adata(adata)

    tcrdist = TcrDistCalculator(adata.uns['organism'])

    num_clones = adata.shape[0]

    all_nbrs = {}
    for nbr_frac in nbr_fracs:
        all_nbrs[nbr_frac] = []

    nndists = []

    for ii in range(num_clones):
        if ii%1000==0:
            print('recalculate_tcrdist_nbrs:', ii, num_clones)
            sys.stdout.flush()
        ii_tcr = tcrs[ii]
        dists = np.array([ tcrdist(ii_tcr, x) for x in tcrs])
        dists[ agroups==agroups[ii] ] = 1e3
        dists[ bgroups==bgroups[ii] ] = 1e3
        for nbr_frac in nbr_fracs: # could do this more efficiently by going in decreasing order, saving partitions...
            num_neighbors = max(1, int(nbr_frac*num_clones))
            ii_nbrs = np.argpartition(dists, num_neighbors-1 )[:num_neighbors]
            all_nbrs[nbr_frac].append(ii_nbrs)
        if nbr_frac_for_nndists is not None:
            num_neighbors = max(1, int(nbr_frac_for_nndists*num_clones))
            lowdists = np.sort( np.partition(dists, num_neighbors-1 )[:num_neighbors] )
            wts = np.linspace(1.0, 1.0/num_nbrs, num_nbrs)
            nndists.append(np.sum( lowdists * wts )/np.sum(wts))

    for nbr_frac in nbr_fracs:
        num_neighbors = max(1, int(nbr_frac*num_clones))
        all_nbrs[nbr_frac] = np.vstack(all_nbrs[nbr_frac])
        assert all_nbrs[nbr_frac].shape == (num_clones, num_neighbors)

    if nbr_frac_for_nndists is None:
        return all_nbrs
    else:
        assert len(nndists) == num_clones
        return all_nbrs, np.array(nndists)


def calculate_tcrdist_nbrs_cpp(
        adata,
        nbr_fracs,
        nbr_frac_for_nndists = None,
        tmpfile_prefix = None
):
    ''' returns dict mapping from nbr_frac to nbrs_tcr

    nbrs exclude self and any clones in same atcr group or btcr group
    '''
    if tmpfile_prefix is None:
        tmpfile_prefix = PurePath('./tmp_nbrs{}'.format(random.randrange(1,10000)))

    print('calculate_tcrdist_nbrs_cpp:', adata.shape, nbr_fracs, tmpfile_prefix)

    agroups, bgroups = setup_tcr_groups(adata)

    agroups_filename = str(tmpfile_prefix) +'_agroups.txt'
    bgroups_filename = str(tmpfile_prefix) +'_bgroups.txt'
    np.savetxt(agroups_filename, agroups, fmt='%d')
    np.savetxt(bgroups_filename, bgroups, fmt='%d')

    tcrs_filename = str(tmpfile_prefix) +'_tcrs.tsv'
    adata.obs['va cdr3a vb cdr3b'.split()].to_csv(tcrs_filename, sep='\t', index=False)

    if os.name == 'posix':
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors')
    else:
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors.exe') 

    if not exists(exe):
        print('need to compile c++ exe:', exe)
        exit(1)

    db_filename = PurePath.joinpath(util.path_to_tcrdist_cpp_db, 'tcrdist_info_{}.txt'.format(adata.uns['organism']))

    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)

    max_nbr_frac = max(nbr_fracs)
    num_nbrs = max(1, int(max_nbr_frac*adata.shape[0]))

    outprefix = str(tmpfile_prefix) +'_calc_tcrdist'

    cmd = '{} -f {} -n {} -d {} -o {} -a {} -b {}'\
    .format(exe, tcrs_filename, num_nbrs, db_filename, outprefix, agroups_filename, bgroups_filename)

    if os.name == 'posix':
        util.run_command(cmd, verbose=True)
    else:
        util.run_command_windows(cmd, verbose=True)

    knn_indices_filename = outprefix+'_knn_indices.txt'
    knn_distances_filename = outprefix+'_knn_distances.txt'

    if not exists(knn_indices_filename) or not exists(knn_distances_filename):
        print('find_neighbors failed:', exists(knn_indices_filename), exists(knn_distances_filename))
        exit(1)

    knn_indices = np.loadtxt(knn_indices_filename, dtype=int)
    knn_distances = np.loadtxt(knn_distances_filename, dtype=float)

    all_nbrs = {}
    all_nbrs[max_nbr_frac] = knn_indices
    for nbr_frac in nbr_fracs:
        if nbr_frac == max_nbr_frac:
            continue
        num_nbrs = max(1, int(nbr_frac*adata.shape[0]))
        assert num_nbrs <= knn_indices.shape[1]
        inds = np.argpartition( knn_distances, num_nbrs-1)[:,:num_nbrs]
        all_nbrs[nbr_frac] = knn_indices[np.arange(adata.shape[0])[:,np.newaxis], inds]

    for filename in [tcrs_filename, agroups_filename, bgroups_filename, knn_indices_filename, knn_distances_filename]:
        os.remove(filename)

    if nbr_frac_for_nndists is None:
        return all_nbrs
    else:
        num_nbrs = max(1, int(nbr_frac_for_nndists*adata.shape[0]))
        assert num_nbrs <= knn_indices.shape[1]
        dists = np.sort( np.partition( knn_distances, num_nbrs-1)[:,:num_nbrs], axis=1)
        wts = np.linspace(1.0, 1.0/num_nbrs, num_nbrs)
        wts /= np.sum(wts)
        nndists = np.sum( dists * wts[np.newaxis,:], axis=1)
        assert nndists.shape==(adata.shape[0],)
        return all_nbrs, nndists


def get_vfam(vgene):
    #assert vgene.startswith('TR') and vgene[3]=='V'
    assert vgene[3]=='V'
    pos = 4
    while pos<len(vgene) and vgene[pos].isdigit():
        pos += 1
    vno = int(vgene[4:pos])
    return '{}V{:d}'.format(vgene[2], vno)


def setup_tcr_cluster_names(adata):

    organism = adata.uns['organism']

    #clusters_gex = adata.obs['clusters_gex']
    clusters_tcr = np.array(adata.obs['clusters_tcr'])

    tcrs = retrieve_tcrs_from_adata(adata)

    num_clusters = np.max(clusters_tcr)+1
    names = []
    for c in range(num_clusters):
        cluster_size = np.sum(clusters_tcr==c)
        ctcrs = [ x for x,y in zip(tcrs,clusters_tcr) if y==c]
        counts = Counter( [ get_vfam(x[0][0]) for x in ctcrs]) +Counter([get_vfam(x[1][0]) for x in ctcrs])
        #print( c, cluster_size, counts.most_common(3))
        top_vfam, top_count = counts.most_common(1)[0]
        eps=1e-3
        if top_count+eps >= 0.75*cluster_size:
            names.append('{}_{}'.format(c, top_vfam))
        elif top_count+eps >= 0.5*cluster_size:
            names.append('{}_{}'.format(c, top_vfam.lower()))
        else:
            if organism=='human': # special hack
                c2 = counts['AV14']+counts['AV38']
                c3 = counts['AV14']+counts['AV38']+counts['AV19']
                if c2 >= 0.75*cluster_size:
                    names.append('{}_AV14+'.format(c, top_vfam))
                elif c3 >= 0.75*cluster_size:
                    names.append('{}_AV14++'.format(c, top_vfam))
                elif c2 >= 0.5*cluster_size:
                    names.append('{}_av14+'.format(c, top_vfam))
                elif c3 >= 0.5*cluster_size:
                    names.append('{}_av14++'.format(c, top_vfam))
                else:
                    names.append(str(c))
            else:
                names.append(str(c))
    print('setup_tcr_cluster_names:', names)
    adata.uns['clusters_tcr_names'] = names


def make_tcrdist_kernel_pcs_file_from_clones_file(
        clones_file,
        organism,
        n_components_in=50,
        kernel=None, # either None (-->default) or 'gaussian'
        gaussian_kernel_sdev=100.0, #unused unless kernel=='gaussian'
        verbose = False,
        outfile = None,
        input_distfile = None,
        output_distfile = None,
        force_Dmax = None,
):
    if outfile is None: # this is the name expected by read_dataset above (with n_components_in==50)
        outfile = '{}_AB.dist_{}_kpcs'.format(clones_file[:-4], n_components_in)

    tcrdist_calculator = TcrDistCalculator(organism)

    df = pd.read_csv(clones_file, sep='\t')

    # in conga we usually also have cdr3_nucseq but we don't need it for tcrdist; we also don't need the jgene but hey
    tcrs = [ ( ( l.va_gene, l.ja_gene, l.cdr3a ), ( l.vb_gene, l.jb_gene, l.cdr3b ) ) for l in df.itertuples() ]
    ids = [ l.clone_id for l in df.itertuples() ]


    if input_distfile is None: ## tcr distances
        print(f'compute tcrdist distance matrix for {len(tcrs)} clonotypes')
        D = np.array( [ tcrdist_calculator(x,y) for x in tcrs for y in tcrs ] ).reshape( (len(tcrs), len(tcrs)) )
    else:
        print(f'reload tcrdist distance matrix for {len(tcrs)} clonotypes')
        D = np.loadtxt(input_distfile)

    if output_distfile is not None:
        np.savetxt( distfile, D.astype(float), fmt='%.1f')

    n_components = min( n_components_in, D.shape[0] )

    print(f'running KernelPCA with {kernel} kernel distance matrix shape= {D.shape} D.max()= {D.max()} force_Dmax= {force_Dmax}')

    pca = KernelPCA(kernel='precomputed', n_components=n_components)

    if kernel is None:
        if force_Dmax is None:
            force_Dmax = D.max()
        gram = np.maximum(0.0, 1 - ( D / force_Dmax ))
    elif kernel == 'gaussian':
        gram = np.exp(-0.5 * (D/gaussian_kernel_sdev)**2 )
    else:
        print('conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file:: unrecognized kernel:', kernel)
        sys.exit(1)

    xy = pca.fit_transform(gram)

    if verbose: #show the eigenvalues
        for ii in range(n_components):
            print( 'eigenvalue: {:3d} {:.3f}'.format( ii, pca.lambdas_[ii]))

    # this is the kpca_file that conga.preprocess.read_dataset is expecting:
    #kpca_file = clones_file[:-4]+'_AB.dist_50_kpcs'
    print( 'writing TCRdist kernel PCs to outfile:', outfile)
    out = open(outfile,'w')

    for ii in range(D.shape[0]):
        out.write('pc_comps: {} {}\n'\
                  .format( ids[ii], ' '.join( '{:.6f}'.format(xy[ii,j]) for j in range(n_components) ) ) )
    out.close()
    return


def condense_clones_file_and_barcode_mapping_file_by_tcrdist(
        old_clones_file,
        new_clones_file,
        tcrdist_threshold,
        organism,
        output_distfile=None,
        force_tcrdist_cpp=False
):

    df = pd.read_csv(old_clones_file, sep='\t')
    N = df.shape[0]

    all_barcodes = pd.read_csv(old_clones_file+'.barcode_mapping.tsv', sep='\t')
    all_barcodes.set_index('clone_id', inplace=True)
    all_barcodes = all_barcodes['barcodes']
    assert type(all_barcodes) is pd.Series


    if output_distfile is None and (force_tcrdist_cpp or (util.tcrdist_cpp_available() and N>5000)):
        tcrdist_threshold = int(tcrdist_threshold+0.001) # cpp tcrdist threshold is integer

        if os.name == 'posix':
            exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors')
        else:
            exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors.exe') 

        db_filename = PurePath.joinpath( util.path_to_tcrdist_cpp_db, 'tcrdist_info_{}.txt'.format(organism) ) 

        outprefix = old_clones_file + '_calc_tcrdist'

        cmd = '{} -f {} -t {} -d {} -o {}'.format(exe, old_clones_file, tcrdist_threshold, db_filename, outprefix)

        if os.name == 'posix':
            util.run_command(cmd, verbose=True)
        else:
            util.run_command_windows(cmd, verbose=True)

        nbr_indices_filename = '{}_nbr{}_indices.txt'.format(outprefix, tcrdist_threshold)
        nbr_distances_filename = '{}_nbr{}_distances.txt'.format(outprefix, tcrdist_threshold)

        if not exists(nbr_indices_filename) or not exists(nbr_distances_filename):
            print('find_neighbors failed:', exists(nbr_indices_filename), exists(nbr_distances_filename))
            exit(1)

        all_nbrs = []
        all_distances = []
        all_smallest_nbr = []
        for line1, line2 in zip(open(nbr_indices_filename,'r'), open(nbr_distances_filename,'r')):
            l1 = line1.split()
            l2 = line2.split()
            assert len(l1) == len(l2)
            ii = len(all_nbrs)
            all_nbrs.append([ii]+[int(x) for x in l1])
            all_distances.append([0.]+[float(x) for x in l2])
            all_smallest_nbr.append(min(all_nbrs[-1]))
        assert len(all_nbrs) == N

        # now do single linkage clustering
        all_smallest_nbr = np.array(all_smallest_nbr)
        while True:
            updated = False
            for ii in range(N):
                nbr = all_smallest_nbr[ii]
                new_nbr = min(nbr, np.min(all_smallest_nbr[all_nbrs[ii]]))
                if nbr != new_nbr:
                    print('update:', nbr, new_nbr)
                    all_smallest_nbr[ii] = new_nbr
                    updated = True
            if not updated:
                break
        # define clusters, choose cluster centers
        clusters = np.array([-1]*N)
        cluster_centers = []
        clusters_set = []
        for ii, nbr in enumerate(all_smallest_nbr):
            if ii==nbr:
                cmask = all_smallest_nbr==nbr
                csize = np.sum(cmask)
                cluster_number = len(clusters_set)
                clusters_set.append(cluster_number)
                clusters[cmask] = cluster_number
                # choose center
                min_avgdist = 1e6
                center = None
                for jj in np.nonzero(cmask)[0]:
                    assert all_smallest_nbr[jj] == nbr
                    avgdist = (np.sum(all_distances[jj]) + (tcrdist_threshold+1.)*(csize-len(all_distances[jj])))/csize
                    if avgdist<min_avgdist:
                        min_avgdist = avgdist
                        center = jj
                assert center is not None
                cluster_centers.append(center) # will be parallel with clusters_set
        assert not np.any(clusters==-1)

    else: # use python tcrdist, compute full distance matrix
        # in conga we usually also have cdr3_nucseq but we don't need it for tcrdist
        tcrs = [ ( ( l.va_gene, l.ja_gene, l.cdr3a ), ( l.vb_gene, l.jb_gene, l.cdr3b ) ) for l in df.itertuples() ]
        tcrdist_calculator = TcrDistCalculator(organism)
        print(f'compute tcrdist distance matrix for {len(tcrs)} clonotypes')
        D = np.array( [ tcrdist_calculator(x,y) for x in tcrs for y in tcrs ] ).reshape( (len(tcrs), len(tcrs)) )

        DT = squareform(D, force='tovector')

        # single linkage clustering of the distance matrix: any clonotypes with dist<tcrdist_threshold
        #  should end up in the same cluster
        Z = hierarchy.single(DT)

        clusters = hierarchy.fcluster(Z, t=tcrdist_threshold, criterion='distance')
        clusters_set = sorted(set(clusters))

        cluster_centers = []
        for c in clusters_set:
            # choose a representative clone based on distance
            cmask = clusters==c
            members = np.nonzero(cmask)[0]
            assert len(members) == np.sum(cmask)
            if len(members) == 1:
                center = members[0]
            else:
                cdist = D[cmask,:][:,cmask]
                dists = np.sum(cdist,axis=1)/(len(members)-1)
                icenter = np.argmin(dists)
                center = members[icenter]
                print('center_avgdist: {:3d} {:7.2f} avg {:7.2f}'\
                      .format(len(members), dists[icenter], np.mean(dists)))
            cluster_centers.append(center)

    print('num_clusters:', len(clusters_set))

    new_clones_dfl = []
    new_bcmap_dfl = []

    for c, center in zip(clusters_set, cluster_centers):
        # choose a representative clone based on distance
        cmask = clusters==c
        members = np.nonzero(cmask)[0]
        cdf = df[cmask]
        center_df = pd.Series( df.iloc[center] )
        clone_size = sum(x.clone_size for _,x in cdf.iterrows())
        center_df.clone_size = clone_size
        new_clones_dfl.append( center_df)
        cbarcodes = []
        for _,row in cdf.iterrows():
            cbarcodes.extend( all_barcodes[row.clone_id].split(',') )
        assert len(cbarcodes) == clone_size
        new_bcmap_dfl.append( dict(clone_id=center_df.clone_id, barcodes=','.join(cbarcodes)))

    new_clones_df = pd.DataFrame(new_clones_dfl)
    new_bcmap_df = pd.DataFrame(new_bcmap_dfl)['clone_id barcodes'.split()] # ensure order?

    new_clones_df.to_csv(new_clones_file, sep='\t', index=False)
    new_bcmap_df.to_csv(new_clones_file+'.barcode_mapping.tsv', sep='\t', index=False)

    if output_distfile is not None:
        new_D = D[cluster_centers,:][:,cluster_centers]
        np.savetxt( output_distfile, new_D.astype(float), fmt='%.1f')


def calc_tcrdist_nbrs_umap_clusters_cpp(
        adata,
        num_nbrs,
        outfile_prefix,
        umap_key_added = 'X_tcrdist_2d',
        cluster_key_added = 'clusters_tcrdist',
        louvain_resolution = None,
        n_components_umap = 2,
):
    tcrs_filename = outfile_prefix+'_tcrs.tsv'
    adata.obs['va cdr3a vb cdr3b'.split()].to_csv(tcrs_filename, sep='\t', index=False)

    if os.name == 'posix':
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors')
    else:
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors.exe') 

    if not exists(exe):
        print('need to compile c++ exe:', exe)
        exit(1)

    db_filename = PurePath.joinpath(util.path_to_tcrdist_cpp_db,\
        'tcrdist_info_{}.txt'.format(adata.uns['organism']) )

    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)

    outprefix = outfile_prefix+'_calc_tcrdist'

    cmd = '{} -f {} -n {} -d {} -o {}'.format(exe, tcrs_filename, num_nbrs, db_filename, outprefix)

    if os.name == 'posix':
        util.run_command(cmd, verbose=True)
    else:
        util.run_command_windows(cmd, verbose=True)

    knn_indices_filename = outprefix+'_knn_indices.txt'
    knn_distances_filename = outprefix+'_knn_distances.txt'

    if not exists(knn_indices_filename) or not exists(knn_distances_filename):
        print('find_neighbors failed:', exists(knn_indices_filename), exists(knn_distances_filename))
        exit(1)

    knn_indices = np.loadtxt(knn_indices_filename, dtype=int)
    knn_distances = np.loadtxt(knn_distances_filename, dtype=float)

    # distances = sc.neighbors.get_sparse_matrix_from_indices_distances_numpy(
    #     knn_indices, knn_distances, adata.shape[0], num_nbrs)

    distances, connectivities = sc.neighbors.compute_connectivities_umap(
        knn_indices, knn_distances, adata.shape[0], num_nbrs)

    if issparse(connectivities):
        from scipy.sparse.csgraph import connected_components
        connected_components = connected_components(connectivities)
        number_connected_components = connected_components[0]
        print('number_connected_components:', number_connected_components)

    ################
    # stash the stuff in adata, stolen from scanpy/neighbors/__init__.py
    #
    adata.uns['neighbors'] = {}
    adata.uns['neighbors']['params'] = {'n_neighbors': num_nbrs, 'method': 'umap'}
    adata.uns['neighbors']['params']['metric'] = 'tcrdist'#metric
    # if metric_kwds:
    #     adata.uns['neighbors']['params']['metric_kwds'] = metric_kwds
    # if use_rep is not None:
    #     adata.uns['neighbors']['params']['use_rep'] = use_rep
    # if n_pcs is not None:
    #     adata.uns['neighbors']['params']['n_pcs'] = n_pcs
    adata.uns['neighbors']['distances'] = distances
    adata.uns['neighbors']['connectivities'] = connectivities

    # as far as I can tell, these are used if there are multiple components in the graph...
    print('putting random pca vectors into adata!!!')
    fake_pca = np.random.randn(adata.shape[0], 10)
    adata.obsm['X_pca'] = fake_pca

    print('running umap', adata.shape)
    sc.tl.umap(adata, n_components=n_components_umap)
    print('DONE running umap')
    adata.obsm[umap_key_added] = adata.obsm['X_umap']

    print('running louvain', adata.shape)
    resolution = 1.0 if louvain_resolution is None else louvain_resolution
    sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key_added)
    adata.obs[cluster_key_added] = np.copy(adata.obs[cluster_key_added]).astype(int)
    print('DONE running louvain', cluster_key_added)

    del adata.obsm['X_pca'] # delete the fake pcas

def calc_tcrdist_matrix_cpp(
        tcrs,
        organism,
        tmpfile_prefix = None,
):
    if tmpfile_prefix is None:
        tmpfile_prefix = PurePath('./tmp_tcrdists{}'.format(random.randrange(1,10000)))

    tcrs_filename = str(tmpfile_prefix) +'_tcrs.tsv'

    df = pd.DataFrame(dict(va=[x[0][0] for x in tcrs], cdr3a=[x[0][2] for x in tcrs],
                           vb=[x[1][0] for x in tcrs], cdr3b=[x[1][2] for x in tcrs]))
    df.to_csv(tcrs_filename, sep='\t', index=False)

    if os.name == 'posix':
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors')
    else:
        exe = PurePath.joinpath( util.path_to_tcrdist_cpp_bin , 'find_neighbors.exe') 

    if not exists(exe):
        print('need to compile c++ exe:', exe)
        exit(1)

    db_filename = PurePath.joinpath(util.path_to_tcrdist_cpp_db, 'tcrdist_info_{}.txt'.format( organism ) )

    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)

    cmd = '{} -f {} --only_tcrdists -d {} -o {}'.format(exe, tcrs_filename, db_filename, tmpfile_prefix)

    if os.name == 'posix':
        util.run_command(cmd, verbose=True)
    else:
        util.run_command_windows(cmd, verbose=True)

    tcrdist_matrix_filename = str(tmpfile_prefix) +'_tcrdists.txt'

    if not exists(tcrdist_matrix_filename):
        print('find_neighbors failed, missing', tcrdist_matrix_filename)
        exit(1)

    D = np.loadtxt(tcrdist_matrix_filename).astype(float)

    for filename in [tcrs_filename, tcrdist_matrix_filename]:
        os.remove(filename)

    return D

"""
New experimental functions for enhancing interative flexibility
"""
def Prep_for_CoNGA(
        adata,
        clones_file,
        output_prefix = None,
        second_clones_file = None,
        write_full_clone_df = True,
):
    ''' prepares adata object generated thru std scanpy workflow for CoNGA analysis.
    
    Uses clones DataFrame generated by make_10x_clone_file or make_10x_clone_file_batch as input.
    Clonotype information gets stored in adata.obs under multiple keys 
    TCRdist kPCA info in adata.obsm under the key 'X_pca_tcr' or 'X_pca_bcr'
    returns adata 
    
    '''
    
    combo_pipe = False
    single_pipe = False
    organism = adata.uns['organism']

    if second_clones_file is None:
        single_pipe = True
    else:
        combo_pipe = True

        
    if single_pipe:
        
        bcmap_file = clones_file +'.barcode_mapping.tsv'
        
        assert exists(clones_file)
        assert exists(bcmap_file)
        
        #these steps might be worth moving up to make_10x_clone_file
        clone_df = pd.read_csv( clones_file, sep='\t')
        barcode_map = pd.read_csv( bcmap_file ,sep='\t')
        
        tag = clone_df.iloc[0,0].split('_')[0]
        
        if tag == 'bcr': 
            organism = organism + '_ig'
        
        #melt the barcode map
        wide_bc = barcode_map['barcodes'].str.split(',', expand=True).copy()
        wide_bc['clone_id'] = barcode_map['clone_id']
    
        long_bc = ( wide_bc.melt( id_vars=['clone_id'] , var_name='cell', value_name='barcode')
                   .dropna()
                   .drop(columns=['cell'])
                   .copy()
                   )
        
        #merge barcodes and clone_df
        clone_df_full = pd.merge( long_bc , clone_df, how='left', on= 'clone_id')
        clone_df_full = clone_df_full.set_index('barcode')
        clone_df_full['rep_type'] = tag
        clone_df_full = clone_df_full.rename(columns={ "va_gene": "va", "ja_gene": "ja", "vb_gene": "vb", "jb_gene": "jb"})
        mapped_cells = len( clone_df_full.index )
        uclones = clone_df_full['clone_id'].nunique()
        
        print(f'{mapped_cells} cells with {tag} info')
        print(f'{uclones} unique {tag} clonotypes')
        
        if write_full_clone_df:
            full_file = clones_file[:-4]+'_full.tsv'
            clone_df_full.to_csv( full_file , sep = '\t')
        
        #add clone info into obs
        adata.obs = adata.obs.join(clone_df_full)
        
        #add kpca to adata.obsm
        
        # add option to tune kernel pcs parameters in the future
        kpcs = make_tcrdist_kernel_pcs_file_from_clones_file_V2( clones_file, organism, output_prefix)
        
        kpcs_df = pd.DataFrame( data= kpcs ) 
        kpcs_df['clone_id'] = clone_df['clone_id']
        clones_tmp = adata.obs['clone_id']
        kpcs_full = pd.merge( clones_tmp, kpcs_df, how='left',on= 'clone_id')
        kpcs_full['clone_id'] = adata.obs.index
    
        X_kpca = kpcs_full.to_numpy()
        adata.obsm[f'X_pca_{tag}'] = X_kpca
        
        return adata
    
    elif combo_pipe:
        
        bcmap_file_1 = clones_file +'.barcode_mapping.tsv'
        bcmap_file_2 = second_clones_file +'.barcode_mapping.tsv'
        
        assert exists(clones_file)
        assert exists(bcmap_file_1)
        assert exists(second_clones_file)
        assert exists(bcmap_file_2)
        
        list_in = [[clones_file, bcmap_file_1] , [second_clones_file, bcmap_file_2] ]
        combo_clones = pd.DataFrame()
        
        for row in list_in:
            clones_file = row[0]
            bcmap_file = row[1]
            
            #these steps might be worth moving up to make_10x_clone_file
            clone_df = pd.read_csv( clones_file, sep='\t')
            barcode_map = pd.read_csv( bcmap_file ,sep='\t')
            tag = clone_df.iloc[0,0].split('_')[0]
            
            organism = adata.uns['organism']
        
            if tag == 'bcr': 
                organism = organism + '_ig'
        
            #melt the barcode map
            wide_bc = barcode_map['barcodes'].str.split(',', expand=True).copy()
            wide_bc['clone_id'] = barcode_map['clone_id']
        
            long_bc = ( wide_bc.melt( id_vars=['clone_id'] , var_name='cell', value_name='barcode')
                       .dropna()
                       .drop(columns=['cell'])
                       .copy()
                       )
        
            #merge barcodes and clone_df
            clone_df_full = pd.merge( long_bc , clone_df, how='left', on= 'clone_id')
            clone_df_full = clone_df_full.set_index('barcode')
            clone_df_full['rep_type'] = tag
            clone_df_full = clone_df_full.rename(columns={ "va_gene": "va", "ja_gene": "ja", "vb_gene": "vb", "jb_gene": "jb"})
            mapped_cells = len( clone_df_full.index )
            uclones = clone_df_full['clone_id'].nunique()
            print(f'{mapped_cells} cells with repertoire info')
            print(f'{uclones} unique clonotypes')
        
            write_full_clone_df = True
            if write_full_clone_df:
                full_file = clones_file[:-4]+'_full.tsv'
                clone_df_full.to_csv( full_file , sep = '\t')
                
            combo_clones = pd.concat([combo_clones , clone_df_full])
        
        
        #filter out collisions
        combo_clones['collisions'] = combo_clones.index.value_counts()
        doubles = combo_clones[ combo_clones.collisions > 1]
        combo_clones_f = combo_clones[ combo_clones.collisions == 1].copy()
        combo_clones_f = combo_clones_f.drop( columns = ['collisions'] ) 
        
        #add to obs
        adata.obs = adata.obs.join(combo_clones_f) 
        
        #filter out
        print( f'{int(len(doubles) /2) } cells containing TCR and BCR information filtered out')
        adata = adata[ ~ adata.obs.index.isin( doubles.index.to_list() ) , :]
        
        
        # run tcrdist and perform kpca on clones tables without doublets
        new_tcr = (adata.obs[ adata.obs.rep_type == 'tcr'].copy()
                   .drop_duplicates(subset = 'clone_id')
                   )
        
        new_bcr = (adata.obs[ adata.obs.rep_type == 'bcr'].copy()
                   .drop_duplicates(subset = 'clone_id')
                   )
        
        
        tick = 0 
        for df in new_tcr, new_bcr:
            
            organism = adata.uns['organism']
        
            if tick == 1: 
                tag = 'bcr'
                organism = organism + '_ig'
            else:
                tag = 'tcr'
        
            #calc TCRdist 
            kpcs = make_tcrdist_kernel_pcs_file_from_clones_file_V2( df, organism)
            kpcs_df = pd.DataFrame( data= kpcs, index = df.clone_id) 
            kpcs_df = kpcs_df.reset_index()
            
            clones_tmp =  ( adata.obs['clone_id'].copy()
                           .reset_index()
                           )
            
            kpcs_full = pd.merge( clones_tmp, kpcs_df, 'left', on= 'clone_id')
            kpcs_full = kpcs_full.set_index('index')
            
            X_kpca = kpcs_full.to_numpy()
            adata.obsm[f'X_pca_{tag}'] = X_kpca
            
            tick += 1

        
        return adata


def make_tcrdist_kernel_pcs_file_from_clones_file_V2(
        df,
        organism,
        output_prefix = None,
        n_components_in=50,
        kernel=None, # either None (-->default) or 'gaussian'
        gaussian_kernel_sdev=100.0, #unused unless kernel=='gaussian'
        verbose = False,
        force_Dmax = None,
):
    

    # in conga we usually also have cdr3_nucseq but we don't need it for tcrdist; we also don't need the jgene but hey
    tcrs = [ ( ( l.va, l.ja, l.cdr3a ), ( l.vb, l.jb, l.cdr3b ) ) for l in df.itertuples() ]
    #ids = [ l.clone_id for l in df.itertuples() ]
    
    print(f'compute tcrdist distance matrix for {len(tcrs)} clonotypes')

    D = calc_tcrdist_matrix_cpp(tcrs, organism, output_prefix )

    n_components = min( n_components_in, D.shape[0] )

    print(f'running KernelPCA with {kernel} kernel distance matrix shape= {D.shape} D.max()= {D.max()} force_Dmax= {force_Dmax}')

    pca = KernelPCA(kernel='precomputed', n_components=n_components)

    if kernel is None:
        if force_Dmax is None:
            force_Dmax = D.max()
        gram = np.maximum(0.0, 1 - ( D / force_Dmax ))
    elif kernel == 'gaussian':
        gram = np.exp(-0.5 * (D/gaussian_kernel_sdev)**2 )
    else:
        print('conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file:: unrecognized kernel:', kernel)
        sys.exit(1)

    xy = pca.fit_transform(gram)

    if verbose: #show the eigenvalues
        for ii in range(n_components):
            print( 'eigenvalue: {:3d} {:.3f}'.format( ii, pca.lambdas_[ii]))

    # expt with saving in numpy format
    #need to make more specific name
    #array_out =  f'{tag}.dist_{n_components}_kpcs.npy' 
    #np.save( array_out,  xy)

    return xy
