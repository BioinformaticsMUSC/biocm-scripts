#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script performs standard processing (after SCTransform) to Seurat objects.\n")
option_list = list(
  make_option(c("-i", "--seurat_file_path"), type = "character",   metavar="character",   default='none',  help="Path for the seurat object."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default='seurat_processed.rds',  help="nCount_RNA threshold."),
  make_option(c("-t", "--integrated"),      type = "character",   metavar="character",   default='false',  help="Whether the seurat object is integrated (true) or not (false)."),
  make_option(c("-r", "--umap_resolution"),      type = "double",   metavar="double",   default='none',  help="UMAP resolution.")
)

opt = parse_args(OptionParser(option_list=option_list))
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
                            npcs = 30, 
                            reduction.name = "pca")
cat('done!\n')
cat('\nFinding Neighbors...')
seurat <- FindNeighbors(object = seurat, 
                                   reduction = "pca", 
                                   dims = 1:30, 
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


Idents(seurat) <- "integrated_snn_res.0.6"
cat('\nRunning UMAP...')
seurat <- RunUMAP(object = seurat, 
                             reduction = "pca", 
                             dims = 1:30)
cat('done!\n')


# Select the RNA counts slot to be the default assay
DefaultAssay(seurat) <- "RNA"
cat('\nRunning Normalization...')
seurat <- NormalizeData(object = seurat, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)

#seurat@meta.data <- seurat@meta.data %>%
#  rownames_to_column("TMP") %>%
#  select(TMP,orig.ident,Genotype,percent.mt,nCount_SCT,nFeature_SCT,integrated_snn_res.0.6) %>%
#  column_to_rownames("TMP")

cat('\n\nsaving seurat object ...\n')
saveRDS(seurat, file = paste0(opt$output_path, '/', opt$seurat_save_name))

cat('\n\nCreating Dim Plots ...\n')
p1 <- DimPlot(object = seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
ggsave(file = paste0(clviz_path, "/08_DIMPLOT_labels.png"))
p2 <- DimPlot(object = seurat, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="sample")
ggsave(file = paste0(clviz_path, "/09_DIMPLOT_sample.png"))
cat('\ndone with script 03_standard_processing.R\n\n')
