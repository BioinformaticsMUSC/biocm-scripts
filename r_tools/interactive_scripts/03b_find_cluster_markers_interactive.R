#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script finds markers for each cluster using Presto.\n")

opt = list(seurat_file_path = "none",
           project_name = "none",
           seurat_save_name = "none",
           output_path = "none",
           integrated = "none",
           umap_resolution = "none")

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

clviz_path = paste0(opt$output_path, '/viz_cluster')

if(!dir.exists(clviz_path)){dir.create(clviz_path,recursive = T)}


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(presto)
  library(AUCell)
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to SCE...')
seurat <- readRDS(opt$seurat_file_path)


##############################
### RUN PRESTO ON CLUSTERS ###
##############################

presto_markers <- wilcoxauc(seurat, "seurat_clusters")

write.csv(presto_markers, file = paste0(opt$output_path, '/presto_markers_clusters.csv'))

top <- top_markers(presto_markers, n = 10)

write.csv(top, file = paste0(opt$output_path, '/top10_presto_markers_clusters.csv'))

seurat@misc[['markers']][['presto']]$all <- presto_markers
seurat@misc[['markers']][['presto']]$top_20 <- top

seurat@misc[['gene_ranks']][['aucell']]$all <- AUCell_buildRankings(seurat@assays[['SCT']]@counts)

saveRDS(seurat, file=paste0(opt$seurat_output_path, '/03b_presto_markers_clusters_', opt$seurat_save_name))

