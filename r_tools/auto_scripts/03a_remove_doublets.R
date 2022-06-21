#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script removes doublets from the seurat object.\n")
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
  library(clustree)
  library(scDblFinder)
  library(SingleCellExperiment)
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to SCE...')
seurat <- readRDS(opt$seurat_file_path)

#convert to SCE for sclDblFinder
sce <- as.SingleCellExperiment(seurat)
cat('done!\n')

############################
### RUN SCLDOUBLETFINDER ###
############################
cat("\nRunning SCDblFinder...")
sce <- scDblFinder(sce, samples = 'sample')

class_counts <- table(sce@colData@listData[["scDblFinder.class"]])
write.csv(class_counts, file = paste0(opt$output_path, '/singlet-doublet_counts.csv'))
print(class_counts)

seurat[['scDblFinder.class']] <- sce@colData@listData[["scDblFinder.class"]]
cat('done!\n')

#plot singlets/doublets
DimPlot(seurat, group.by = "scDblFinder.class")
ggsave(file = paste0(clviz_path, "/10_DIMPLOT_doublets.png"))


#######################
### FILTER DOUBLETS ###
#######################

#filter out doublets
cat('\nFiltering doublets and processing ...')
seurat_no_dbl <- subset(seurat, subset = scDblFinder.class == 'singlet')
cat('done!\n')
#########################
### STANDARD WORKFLOW ###
#########################
if (opt$integrated == 'true'){
  DefaultAssay(seurat_no_dbl) <- 'integrated'
} else {
  DefaultAssay(seurat_no_dbl) <- 'RNA'
}

#seurat_no_dbl
cat('\nRunning Standard workflow now without doublets\n\n')
cat('\nRunning PCA...')
seurat_no_dbl <- RunPCA(object = seurat_no_dbl,
                        features=NULL, 
                        weight.by.var = TRUE, 
                        ndims.print = 1:5, 
                        nfeatures.print = 30, 
                        npcs = 30, 
                        reduction.name = "pca")
cat('done!\n')
cat('\nFinding Neighbors...')
seurat_no_dbl <- FindNeighbors(object = seurat_no_dbl, 
                               reduction = "pca", 
                               dims = 1:30, 
                               nn.eps = 0.5)
cat('done!\n')
cat('\nFinding clusters using resolution(s):',opt$umap_resolution)
seurat_no_dbl <- FindClusters(object = seurat_no_dbl,
                              resolution = opt$umap_resolution, 
                              algorithm = 1,
                              n.iter = 1000)
cat('done!\n')
cat('\nRunning UMAP...')
seurat_no_dbl <- RunUMAP(object = seurat_no_dbl, 
                         reduction = "pca", 
                         dims = 1:30)
cat('done!\n')


cat('\n\nCreating Dim Plots ...\n')
p1 <- DimPlot(object = seurat_no_dbl, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
ggsave(file = paste0(clviz_path, "/11_DIMPLOT_labels_filt.png"))
p2 <- DimPlot(object = seurat_no_dbl, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="sample")
ggsave(file = paste0(clviz_path, "/12_DIMPLOT_sample_filt.png"))

cat('\nSaving filtered Seurat object ...')
saveRDS(seurat_no_dbl, file = paste0(opt$output_path, '/', opt$seurat_save_name))
cat('done!\n')
