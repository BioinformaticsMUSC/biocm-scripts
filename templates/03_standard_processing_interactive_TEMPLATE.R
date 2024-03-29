#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script performs standard processing (after SCTransform) to Seurat objects.\n")

opt = list(seurat_file_path = "{{seurat_file_path}}",
           project_name = "{{project_name}}",
           seurat_save_name = "{{seurat_save_name}}",
           seurat_output_path = "{{seurat_output_path}}",
           output_path = "{{output_path}}",
           integrated = "{{integrated}}",
           umap_resolution = "{{umap_resolution}}")

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

clviz_path = paste0(opt$output_path, '/viz_cluster')

if(!dir.exists(clviz_path)){dir.create(clviz_path,recursive = T)}

if(opt$umap_resolution == 'none'){
  opt$umap_resolution <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2)
}

##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(clustree)
  library(harmony)
}))


####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
seurat <- readRDS(opt$seurat_file_path)

if (opt$integrated == 'true'){
  DefaultAssay(seurat) <- 'integrated'
} else {
  DefaultAssay(seurat) <- 'RNA'
}

cat('\nRunning PCA...')
seurat <- RunPCA(object = seurat, 
                            features=NULL, 
                            weight.by.var = TRUE, 
                            ndims.print = 1:5, 
                            nfeatures.print = 30, 
                            npcs = 50, 
                            reduction.name = "pca")
cat('done!\n')

cat('\nRunning Harmony for batch correction...')
seurat <- RunHarmony(seurat, 'sample')
cat('done!\n')

cat('\nFinding Neighbors...')
seurat <- FindNeighbors(object = seurat, 
                                   reduction = "harmony", 
                                   dims = 1:50, 
                                   nn.eps = 0.5)
cat('done!\n')
cat('\nFinding clusters using resolution(s):',opt$umap_resolution)
seurat <- FindClusters(object = seurat, 
                                  resolution = opt$umap_resolution,
                                  algorithm = 1,
                                  n.iter = 1000)
cat('done!\n')

clustree(seurat@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
ggsave(file = paste0(clviz_path, "/07_clustree.pdf"), height = 12, width = 8)
cat("\n\nClustree created - use this to choose a resolution if needed. Trying 0.6 now...\n")

##########################################
### PAUSE HERE TO CHOOSE RESOLUTION ######
##########################################

##Change this to select the resolution
Idents(seurat) <- "integrated_snn_res.0.6"


cat('\nRunning UMAP...')
seurat <- RunUMAP(object = seurat, 
                             reduction = "harmony", 
                             dims = 1:50)
cat('done!\n')

cat('\n\nCreating Dim Plots ...\n')
p1 <- DimPlot(object = seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
ggsave(file = paste0(clviz_path, "/08_DIMPLOT_labels.pdf"))
p2 <- DimPlot(object = seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="sample")
ggsave(file = paste0(clviz_path, "/09_DIMPLOT_sample.pdf"))

cat('\n\nsaving seurat object ...\n')
saveRDS(seurat, file = paste0(opt$seurat_output_path, '/', opt$seurat_save_name))


cat('\ndone with script 03_standard_processing.R\n\n')
