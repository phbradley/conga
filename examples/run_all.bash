# Run this after running setup_all.bash and changing these next few lines if necessary
#
# this will launch all three calculations concurrently on your machine. If you want them to run
# sequentially, remove the '&' symbols at the end of the lines that start "$cmd"
#
# change this to point to the correct python executable
PYTHON=/home/pbradley/anaconda2/envs/scanpy_new/bin/python

# change this to point to the place where you put the conga repository
CONGA_REPO=/home/pbradley/gitrepos/conga/

# change this to point to wherever you downloaded the 10X example datasets
# this should be where the clones files and associated data was put by the setup_all.bash script
DATADIR="./conga_example_datasets/"

################################################################################

SCRIPT="${CONGA_REPO}scripts/run_conga.py"

# Mouse PBMC
OUTPREFIX="tcr_mm_pbmc"
CLONES="${DATADIR}vdj_v1_mm_balbc_pbmc_t_filtered_contig_annotations_tcrdist_clones.tsv"
GEX="${DATADIR}vdj_v1_mm_balbc_pbmc_5gex_filtered_gene_bc_matrices_h5.h5"

cmd="nice $PYTHON $SCRIPT --all --organism mouse --clones_file $CLONES  --gex_data $GEX --gex_data_type 10x_h5 --outfile_prefix ${OUTPREFIX}"
echo ; echo $cmd
$cmd > ${OUTPREFIX}.log 2> ${OUTPREFIX}.err &


# Human PBMC
OUTPREFIX="tcr_hs_pbmc"
CLONES="${DATADIR}vdj_v1_hs_pbmc3_t_filtered_contig_annotations_tcrdist_clones.tsv"
GEX="${DATADIR}vdj_v1_hs_pbmc3_5gex_filtered_gene_bc_matrices_h5.h5"

cmd="nice $PYTHON $SCRIPT --all --organism human --clones_file $CLONES  --gex_data $GEX --gex_data_type 10x_h5 --outfile_prefix ${OUTPREFIX}"
echo ; echo $cmd
$cmd > ${OUTPREFIX}.log 2> ${OUTPREFIX}.err &


# Human melanoma B cells, merged clonotype definitions based on distance
OUTPREFIX="bcr_hs_melanoma"
CLONES="${DATADIR}sc5p_v1p1_hs_melanoma_10k_b_filtered_contig_annotations_tcrdist_clones_condensed.tsv"
GEX="${DATADIR}sc5p_v1p1_hs_melanoma_10k_filtered_feature_bc_matrix.h5"

cmd="nice $PYTHON $SCRIPT --all --organism human_ig --clones_file $CLONES --gex_data $GEX --gex_data_type 10x_h5 --outfile_prefix $OUTPREFIX"
echo ; echo $cmd
$cmd > ${OUTPREFIX}.log 2> ${OUTPREFIX}.err &


