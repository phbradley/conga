################################################################################
##
## these are strings used as keys when we store results various places
## like in adata.uns['conga_results']
##
## (maybe) useful if we need to repeat the same string in multiple places
##
## also they are strings that appear in the default figure filenames
##  typically: <outfile_prefix>_<tag>.png

# table tags:
GRAPH_VS_GRAPH_STATS = 'graph_vs_graph_stats'
GRAPH_VS_GRAPH = 'graph_vs_graph'
TCR_DB_MATCH = 'tcr_db_match'
TCR_CLUMPING = 'tcr_clumping'

TCR_GRAPH_VS_GEX_FEATURES = 'tcr_graph_vs_gex_features'
TCR_GENES_VS_GEX_FEATURES = 'tcr_genes_vs_gex_features'
GEX_GRAPH_VS_TCR_FEATURES = 'gex_graph_vs_tcr_features'

HOTSPOT_FEATURES = 'hotspot_features'

# figure tags:
GRAPH_VS_GRAPH_LOGOS = 'graph_vs_graph_logos'
TCR_GRAPH_VS_GEX_FEATURES_PLOT = 'tcr_graph_vs_gex_features_plot'
TCR_GRAPH_VS_GEX_FEATURES_PANELS = 'tcr_graph_vs_gex_features_panels'
#TCR_GRAPH_VS_GEX_FEATURES_CLUSTERMAP = 'tcr_graph_vs_gex_features_clustermap'
TCR_GENES_VS_GEX_FEATURES_PANELS = 'tcr_genes_vs_gex_features_panels'
TCR_DB_MATCH_PLOT = 'tcr_db_match_plot'
GEX_GRAPH_VS_TCR_FEATURES_PLOT = 'gex_graph_vs_tcr_features_plot'
GEX_GRAPH_VS_TCR_FEATURES_PANELS = 'gex_graph_vs_tcr_features_panels'
GRAPH_VS_FEATURES_GEX_CLUSTERMAP = 'graph_vs_features_gex_clustermap'
GRAPH_VS_FEATURES_TCR_CLUSTERMAP = 'graph_vs_features_tcr_clustermap'
TCR_CLUMPING_LOGOS = 'tcr_clumping_logos'
GRAPH_VS_SUMMARY = 'graph_vs_summary'
HOTSPOT_GEX_UMAP  = 'hotspot_gex_umap'
HOTSPOT_TCR_UMAP  = 'hotspot_tcr_umap'
HOTSPOT_GEX_CLUSTERMAP  = 'hotspot_gex_clustermap'
HOTSPOT_TCR_CLUSTERMAP  = 'hotspot_tcr_clustermap'
GEX_CLUSTERS_TCRDIST_TREES = 'gex_clusters_tcrdist_trees'
TCR_CLUSTERS_TCRDIST_TREES = 'tcr_clusters_tcrdist_trees'
CONGA_THRESHOLD_TCRDIST_TREE = 'conga_threshold_tcrdist_tree'
BATCH_UMAPS = 'batch_umaps'

# this gets added to table_tag or figure_tag to maek the key for storing help
#   messages in adata.uns['conga_results']
#
HELP_SUFFIX = '_help'
