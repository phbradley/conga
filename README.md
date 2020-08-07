# Clonotype Neighbor Graph Analysis (CoNGA) -- pre-beta

This repository contains the `conga` python package and associated scripts
and workflows. `conga` was developed to detect correlation between
T cell gene expression profile and TCR sequence in single-cell datasets.
`conga` is in active development right now so the interface may change in
the next few months. Questions and requests can be directed to `pbradley` at `fredhutch` dot `org`.

# Running

Running `conga` on a single-cell dataset is a two- (or more) step process, as outlined below.
Python scripts are provided in the `scripts/` directory but analysis steps can also be accessed interactively
in jupyter notebooks (for example, [a simple pipeline](simple_conga_pipeline.ipynb) and
[Seurat to conga](Seurat_to_Conga.ipynb) in the top directory of this repo)
or in your own python scripts through the interface in the `conga` python package.
The examples below and in the jupyter notebooks feature publicly available data from 10X Genomics,
which can be downloaded
using these links: [GEX data](https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_pbmc_5gex) and
[TCR data](https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_pbmc_t).


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
python conga/scripts/run_conga.py --graph_vs_graph --gex_data data/vdj_v1_hs_pbmc_5gex_filtered_gene_bc_matrices_h5.h5 --gex_data_type 10x_h5 --clones_file vdj_v1_hs_pbmc_t_filtered_contig_annotations_tcrdist_clones.tsv --organism human --outfile_prefix tmp_hs_pbmc
```

3. **RE-ANALYZE**: Step 2 will generate a processed `.h5ad` file that contains all the gene expression
and TCR sequence information along with the results of clustering and dimensionality reduction. It can then
be much faster to perform subsequent re-analysis or downstream analysis by "restarting" from those files.

```
python conga/scripts/run_conga.py --restart tmp_hs_pbmc_final.h5ad --graph_vs_tcr_features --graph_vs_gex_features --outfile_prefix tmp_hs_pbmc_restart
```


# Installation

`conga` relies heavily on the wonderful `scanpy` python package for single-cell analysis. See the `scanpy`
instructions for installation: <https://scanpy.readthedocs.io/en/stable/installation.html>.
We highly recommend using anaconda/miniconda for managing python environments. The calculations in the
`conga` manuscript were conducted with the following package versions:

```
scanpy==1.4.3 anndata==0.6.18 umap==0.3.9 numpy==1.16.2 scipy==1.2.1 pandas==0.24.1 scikit-learn==0.20.2 statsmodels==0.9.0 python-igraph==0.7.1 louvain==0.6.1
```

which might possibly be installed with the following `conda` command:
```
conda create -n conga_classic_env ipython python=3.6 scanpy=1.4.3 umap-learn=0.3.9
```


We've also been able to re-run everything, albeit with some numerical changes, with a current (2020-05-25) scanpy
installation and these package versions:
```
scanpy==1.5.1 anndata==0.7.3 umap==0.4.3 numpy==1.17.5 scipy==1.4.1 pandas==1.0.3 scikit-learn==0.23.1 statsmodels==0.11.1 python-igraph==0.8.2 louvain==0.6.1 leidenalg==0.8.0
```

Which was installed with the following `conda` commands (following the `scanpy` docs):
```
conda create -n conga_new_env ipython python=3.6
conda activate conga_new_env   (or source activate conga_new_env)
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg
conda install -c conda-forge louvain
pip install scanpy
```

(And consider also adding `conda install -c conda-forge notebook` for Jupyter notebook stuff.)

Preliminary results suggest that, at least with default clustering parameters, the older `louvain`
clustering algorithm seems to give slightly 'better' results than the newer `leiden` algorithm,
ie finds a few more GEX/TCR associations, probably because there seem to be fewer, larger clusters.
If the `louvain` package is installed `conga` will use that. 


# svg to png
The `conga` image-making pipeline requires an svg to png conversion. There seem to be a variety of
options for doing this, with the best choice being somewhat platform dependent. We've had good luck with
ImageMagick `convert` (on linux) and Inkscape (on mac). The conversion is handled in the file
`conga/convert_svg_to_png.py`, so you can modify that file if things are not working and you have
a tool installed; `conga` may not be looking in the right place. 

If you are having trouble and are using anaconda/miniconda, you could try
`conda install -c conda-forge imagemagick` in the relevant conda environment.
