#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script removes doublets from the seurat object.\n")

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
ggsave(file = paste0(clviz_path, "/10_DIMPLOT_doublets.pdf"))


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
ggsave(file = paste0(clviz_path, "/11_DIMPLOT_labels_filt.pdf"))
p2 <- DimPlot(object = seurat_no_dbl, reduction = "umap", label = FALSE, pt.size = 0.5, split.by="sample")
ggsave(file = paste0(clviz_path, "/12_DIMPLOT_sample_filt.pdf"), width = 12)

cat('\nSaving filtered Seurat object ...')
saveRDS(seurat_no_dbl, file = paste0(opt$seurat_output_path, '/', opt$seurat_save_name))
cat('done!\n')
