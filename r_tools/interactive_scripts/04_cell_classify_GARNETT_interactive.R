#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script predicts cell types using Garnett.\n")


opt = list(seurat_file_path='none',
                 project_name='none',
                 output_path='none',
                 seurat_save_name='04_project_garnett.rds',
                 model_file_path="none",
                 marker_file_path="none",
                 species='mouse')

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

if (opt$model_file_path != "none"){
  cls <- readRDS(opt$model_file_path)
} else {
  cls <- train_cell_classifier(cds = monocle_cds,
                               marker_file = opt$marker_file_path,
                               db=org.Mm.eg.db,
                               cds_gene_id_type = "SYMBOL",
                               num_unknown = 50,
                               marker_file_gene_id_type = "SYMBOL")
}




cat('\nClassifying cells...')
monocle_classified <- classify_cells(monocle_cds, 
                                     cls,
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


