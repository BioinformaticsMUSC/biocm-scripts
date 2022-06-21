#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script finds markers for each cluster using Presto.\n")
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
option_list = list(
  make_option(c("-i", "--seurat_file_path"), type = "character",   metavar="character",   default='none',  help="Path for the seurat object."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default='seurat_processed.rds',  help="nCount_RNA threshold."),
  make_option(c("-t", "--integrated"),      type = "character",   metavar="character",   default='false',  help="Whether the seurat object is integrated (true) or not (false)."),
  make_option(c("-r", "--umap_resolution"),      type = "double",   metavar="double",   default=0.5,  help="UMAP resolution.")
)

opt = parse_args(OptionParser(option_list=option_list))
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



