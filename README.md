# Clonotype Neighbor Graph Analysis (CoNGA) -- pre-beta

This repository contains the `conga` python package and associated scripts
and workflows. `conga` was developed to detect correlation between
T cell gene expression profile and TCR sequence in single-cell datasets.
`conga` is in active development right now so the interface may change in
the next few months. Questions and requests can be directed to `pbradley` at `fredhutch` dot `org`.

# Running

Running `conga` on a single-cell dataset is a two- (or more) step process (see examples in `workflows/`).
Python scripts are provided in the `scripts/` directory but analysis steps can also be accessed interactively
in jupyter notebooks or in your own python scripts through the interface in the `conga` python package.

1. **SETUP**: The TCR data is converted to a form that can be read by `conga` and then
a matrix of `TCRdist` distances is computed. KernelPCA is applied to this distance
matrix to generate a PC matrix that can be used in clustering and dimensionality reduction. This
is accomplished with the python script `scripts/setup_10x_for_conga.py` for 10x datasets. For example:

```
python conga/scripts/setup_10x_for_conga.py --filtered_contig_annotations_csvfile vdj_v1_hs_pbmc_t_filtered_contig_annotations.csv --organism human
```

2. **ANALYZE**: The `scripts/run_conga.py` script has an implementation of the main pipeline and can be run
as follows:

```
python conga/scripts/run_conga.py --gex_data data/vdj_v1_hs_pbmc_5gex_filtered_gene_bc_matrices_h5.h5  --gex_data_type 10x_h5 --clones_file vdj_v1_hs_pbmc_t_filtered_contig_annotations_tcrdist_clones.tsv --find_nbrhood_overlaps --organism human  --outfile_prefix tmp_hs_pbmc
```


# Installation



# Dependencies



