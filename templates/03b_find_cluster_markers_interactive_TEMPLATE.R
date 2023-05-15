#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script finds markers for each cluster using Presto.\n")

opt = list(seurat_file_path = "{{seurat_file_path}}",
           project_name = "{{project_name}}",
           output_path = "{{output_path}}")

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)


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
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to SCE...')
seurat <- readRDS(opt$seurat_file_path)

seurat$seurat_clusters <- seurat$integrated_snn_res.0.6


##############################
### RUN PRESTO ON CLUSTERS ###
##############################

presto_markers <- wilcoxauc(seurat, "seurat_clusters")

write.csv(presto_markers, file = paste0(opt$output_path, '/presto_markers_clusters.csv'))

top <- top_markers(presto_markers, n = 20)

write.csv(top, file = paste0(opt$output_path, '/top20_presto_markers_clusters.csv'))



