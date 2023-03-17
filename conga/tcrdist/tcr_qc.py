# Used for cleaning up for metaconga. could be dropped later.
from conga.tcrdist.all_genes import all_genes
import re

def tcr_qc(adata, species):
    assert species in ['mouse', 'human', 'mouse_gd', 'human_gd', 'human_ig']
    gene_set = list(all_genes[species])

    #lower case the nucseq  just in case. does it care?
    for nuc in ['cdr3a_nucseq', 'cdr3a_nucseq']:
        adata.obs[nuc] = adata.obs[nuc].str.lower()
    
    #repair missing allele info and filter missing genes
    for gene in ['va','ja','vb','jb']:
        adata.obs[gene] = adata.obs[gene].apply( lambda x: x if re.search(r'[*]0.$',x) is not None else x + '*01' )
        good_genes = adata.obs[gene].isin(gene_set)
        n_out = adata.obs.shape[0] - good_genes.value_counts()[0]
        print(f'Removing {n_out} chains with missing {gene}')
        adata = adata[good_genes, :].copy()

    #filter out non-functional, bad endings, and short/long
    for chain in ['cdr3a','cdr3b']:
        good_chains = (adata.obs[chain].str.contains('F|W$') & 
                       ~adata.obs[chain].str.contains('[*]|_') & 
                       adata.obs[chain].str.len().between(6,22) ) 
        n_outc = adata.obs.shape[0] - good_chains.value_counts()[0]
        print(f'Removing {n_outc} bad {chain}')
        adata = adata[good_chains, :].copy()



    return adata