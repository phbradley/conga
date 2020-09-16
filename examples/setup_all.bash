# change this to point to the correct python executable
PYTHON=/home/pbradley/anaconda2/envs/scanpy_new/bin/python

# change this to point to the place where you put the conga repository
CONGA_REPO=/home/pbradley/gitrepos/conga/

# change this to point to wherever you downloaded the 10X example datasets
DATADIR="./conga_example_datasets/"

########################################

SCRIPT="${CONGA_REPO}scripts/setup_10x_for_conga.py"

# Mouse PBMC
CONTIGS="${DATADIR}vdj_v1_mm_balbc_pbmc_t_filtered_contig_annotations.csv"

cmd="nice $PYTHON $SCRIPT --organism mouse --filtered_contig_annotations_csvfile $CONTIGS"
echo "Running:" $cmd
$cmd > setup1.log 2> setup1.err &

# Human PBMC
CONTIGS="${DATADIR}vdj_v1_hs_pbmc3_t_filtered_contig_annotations.csv"

cmd="nice $PYTHON $SCRIPT --organism human --filtered_contig_annotations_csvfile $CONTIGS"
echo "Running:" $cmd
$cmd > setup2.log 2> setup2.err &

# Human melanoma B cells -- use stricter tcrdist based clonotypes
CONTIGS="${DATADIR}sc5p_v1p1_hs_melanoma_10k_b_filtered_contig_annotations.csv"

cmd="nice $PYTHON $SCRIPT --organism human_ig --filtered_contig_annotations_csvfile $CONTIGS --condense_clonotypes_by_tcrdist --tcrdist_threshold_for_condensing 50"
echo "Running:" $cmd
$cmd > setup4.log 2> setup4.err &

