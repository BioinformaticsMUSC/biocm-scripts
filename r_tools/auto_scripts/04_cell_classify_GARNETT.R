#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script predicts cell types using Garnett.\n")
option_list = list(
  make_option(c("-i", "--seurat_file_path"), type = "character",   metavar="character",   default='none',  help="Path for the seurat object."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default='seurat_processed_celltypes.rds',  help="nCount_RNA threshold."),
  make_option(c("-m", "--model_file_path"),      type = "character",   metavar="character",   default='none',  help="File path for pretrained classifier"),
  make_option(c("-f", "--marker_file_path"),      type = "character",   metavar="character",   default='none',  help="File path for marker file."),
  make_option(c("-d", "--species"),      type = "character",   metavar="character",   default="none",  help="human or mouse.")
)

opt = parse_args(OptionParser(option_list=option_list))
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
  library(clustree)
  library(SingleCellExperiment)
  library(garnett)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
}))

if (opt$species == 'human') {
  db <- org.Hs.eg.db
} else if (opt$species == 'mouse') {
  db <- org.Mm.eg.db
} else {
  stop('Species not found - please specify human or mouse')
}


####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to CellDataSet...')
seurat <- readRDS(opt$seurat_file_path)

#convert to SCE for sclDblFinder
sce <- as.SingleCellExperiment(seurat)

data <- as(as.matrix(seurat@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seurat@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

cat('done!\n')

##############################
### READ AND CHECK MARKERS ###
##############################
cat('\nChecking markers...')
monocle_cds <- estimateSizeFactors(monocle_cds)

marker_check <- check_markers(monocle_cds, opt$marker_file_path,
                              db = db,
                              cds_gene_id_type = 'SYMBOL',
                              marker_file_gene_id_type = 'SYMBOL')

plot_markers(marker_check)
ggsave(filename = paste0(opt$output_path, '/marker_check.png'))
cat('done!\n')

######################
### CLASSIFY CELLS ###
######################
cat('\nLoading Classifier...\n')
pretrained_cls <- readRDS(opt$model_file_path)

cat('\nClassifying cells...')
monocle_classified <- classify_cells(monocle_cds, 
                                     pretrained_cls,
                                     db = db,
                                     rank_prob_ratio = 1.5,
                                     cluster_extend = TRUE,
                                     cds_gene_id_type = "SYMBOL")

seurat[['garnett_celltype']] <- monocle_classified@phenoData@data$cell_type
cat('done!\n')

DimPlot(seurat, group.by = 'garnett_celltype')
ggsave(filename = paste0(opt$output_path, '/DIMPLOT_celltypes.png'))

Labels <- table(seurat@meta.data$garnett_celltype, seurat@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  filter(Var1 != "Unknown") %>% 
  group_by(Var2) %>%
  filter(Freq == max(Freq)) %>%
  as.data.frame()

colnames(Labels) <- c("Cell","seurat_clusters","Freq")
print(Labels)

write.table(Labels, "Integrated_Labels_clusters.txt", sep="\t",quote=F,row.names=F)

cat('\nSaving seurat object with garnett cell annotations...')
saveRDS(seurat, file = paste0(opt$output_path, '/', opt$seurat_save_name))


