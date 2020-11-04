import argparse
from os.path import exists
import sys
import os
import numpy as np

parser = argparse.ArgumentParser(description="Merge multiple datasets and generate a single clones_file and gex_data file for subsequent CoNGA analysis. Each individual dataset should be previously setup to run through conga. Use the --samples argument to provide a tsv file with three rows that give, for each dataset, the location of the clones file and gex datafile and the format of the gex data" )

parser.add_argument('--samples', required=True, help='tsvfile with 3 columns, "clones_file" "gex_data" "gex_data_type" corresponding to the arguments for run_conga.py')
parser.add_argument('--output_clones_file', required=True)
parser.add_argument('--output_gex_data', required=True, help='The format of this output file will be "h5ad" ie scanpy h5')
parser.add_argument('--output_distfile', help='(optional) Save the tcrdist distance matrix in numpy savetxt format to this file')
parser.add_argument('--organism', choices=['mouse', 'human', 'mouse_gd', 'human_gd', 'human_ig'], required=True)
parser.add_argument('--condense_clonotypes_by_tcrdist', action='store_true', help='Merge clonotypes with TCRdist distances less than the threshold specified by --tcrdist_threshold_for_condensing. This can be useful for BCR data to merge families of clonally-related cells')
parser.add_argument('--tcrdist_threshold_for_condensing', type=float, default=50.)
parser.add_argument('--no_tcrdists', action='store_true', help='Don\'t compute tcrdists or kernel PCs; instead generate a random matrix of kernel PCs. This might be useful for preprocessing very big GEX datasets to isolate subsets of interest')
parser.add_argument('--no_kpcs', action='store_true')
parser.add_argument('--force_tcrdist_cpp', action='store_true')
parser.add_argument('--batch_keys', type=str, nargs='*')


args = parser.parse_args()

if args.no_tcrdists:
    assert not args.condense_clonotypes_by_tcrdist
    assert not args.output_distfile

# put this after arg parsing because it's so dang slow
sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) ) # so we can import conga
import conga
import conga.preprocess
import pandas as pd
import scanpy as sc


df = pd.read_csv(args.samples, sep='\t')

all_data = []

add_batch_suffix = (df.shape[0]>1)

for ii, row in df.iterrows():
    bcmap_file = row.clones_file+'.barcode_mapping.tsv'
    assert exists(bcmap_file)

    clones_df = pd.read_csv(row.clones_file, sep='\t')
    bcmap_df = pd.read_csv(bcmap_file, sep='\t')
    adata = conga.preprocess.read_adata(row.gex_data, row.gex_data_type, gex_only=False)

    if args.batch_keys is not None:
        for k in args.batch_keys:
            assert k in df.columns
            val = int(row[k]) # conga expects integers
            adata.obs[k] = np.array([val]*adata.shape[0])
            print(f'saving batch_key {k} value {val} for sample row {row.clones_file}')

    if add_batch_suffix:
        suffix = '-'+str(ii)
        clones_df.clone_id += suffix
        bcmap_df.clone_id += suffix
        new_barcodes = []
        for bcs in bcmap_df.barcodes:
            barcodes = bcs.split(',')
            new_barcodes.append(','.join([x+suffix for x in barcodes]))
        bcmap_df.barcodes = new_barcodes

    adata.obs['batch_gex_data'] = row.gex_data
    adata.obs['batch_clones_file'] = row.clones_file

    all_data.append( [clones_df, bcmap_df, adata ] )
    print(adata.shape, row.gex_data)
    sys.stdout.flush()

new_clones_df = pd.concat([x[0] for x in all_data])
new_bcmap_df = pd.concat([x[1] for x in all_data])
if len(all_data)==1:
    new_adata = all_data[0][2] # no batch obs key or batch suffix!
else:
    new_adata = all_data[0][2].concatenate(*[x[2] for x in all_data[1:]])

# this all assumes that when scanpy concatenates it adds '-N' to the Nth datasets barcodes
if args.condense_clonotypes_by_tcrdist:
    tmpfile = args.output_clones_file+'.uncondensed.tsv'
    new_clones_df.to_csv(tmpfile, sep='\t', index=False)
    new_bcmap_df.to_csv(tmpfile+'.barcode_mapping.tsv', sep='\t', index=False)
    if args.no_kpcs:
        input_distfile = None
    else:
        input_distfile = args.output_clones_file[:-4]+'_AB.dist'
    conga.preprocess.condense_clones_file_and_barcode_mapping_file_by_tcrdist(
        tmpfile, args.output_clones_file, args.tcrdist_threshold_for_condensing, args.organism,
        output_distfile=input_distfile, force_tcrdist_cpp=args.force_tcrdist_cpp)
else:
    print('writing', new_clones_df.shape[0], 'clonotypes to merged clones_file', args.output_clones_file)
    new_clones_df.to_csv(args.output_clones_file, sep='\t', index=False)
    new_bcmap_df.to_csv(args.output_clones_file+'.barcode_mapping.tsv', sep='\t', index=False)
    input_distfile=None

print('writing anndata object of shape:', new_adata.shape, 'to file:', args.output_gex_data)
new_adata.write_h5ad(args.output_gex_data)

# now we have to compute the merged tcrdist distances and kpcs
if args.no_tcrdists:
    # the filename and output code is stolen from
    #  conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file
    num_clones = new_clones_df.shape[0]
    clone_ids = list(new_clones_df.clone_id)
    num_components = 50
    outfile = '{}_AB.dist_{}_kpcs'.format(args.output_clones_file[:-4], num_components)
    kpcs = np.random.rand(num_clones, min(num_clones, num_components))

    print('writing fake (random) kernel PCs of shape', kpcs.shape, 'to', outfile)
    out = open(outfile,'w')
    for ii in range(num_clones):
        out.write('pc_comps: {} {}\n'\
                  .format( clone_ids[ii], ' '.join( '{:.6f}'.format(x) for x in kpcs[ii,:])))
    out.close()
elif args.no_kpcs:
    pass
else: # the usual route
    conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file(
        args.output_clones_file, args.organism, input_distfile=input_distfile,
        output_distfile=args.output_distfile )

if args.output_distfile is None and input_distfile is not None:
    os.remove(input_distfile)

print('DONE')
