#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script loads Seurat Object from 10X data\n")



print(t(t(unlist(opt))))

opt = list(data_dir = "none",
           project_name = "none",
           sample_name = "none",
           h5_file_path = "none",
           output_path = "none",
           min_cells = "none",
           min_feats = "none",
           species = "none")

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
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
if (opt$data_dir == 'none' && opt$h5_file_path == 'none') {
  stop('Data file path not detected.')
} else if (opt$h5_file_path != "none") {
  seurat_data <- Read10X_h5(opt$h5_file_path)
} else {
  seurat_data <- Read10x(opt$data_dir)
}

cat("\nFiltering for Protein Coding Genes\n")

if (opt$species == 'human') {
  load('/Users/bryanwgranger/biocm/biocm-tools/r_tools/protein_coding_genes/HgProteinCodingGenes.rda')
  cat("Prefiltering genes total: ", nrow(seurat_data))
  seurat_data <- seurat_data[rownames(seurat_data) %in% HgProteinCodingGenes,]
  cat("\nPostfiltering genes total: ", nrow(seurat_data), "\n")
} else if (opt$species == 'mouse') {
  load('/Users/bryanwgranger/biocm/biocm-tools/r_tools/protein_coding_genes/MgProteinCodingGenes.rda')
  cat("Prefiltering genes total: ", nrow(seurat_data))
  seurat_data <- seurat_data[rownames(seurat_data) %in% HgProteinCodingGenes,]
  cat("\nPostfiltering genes total: ", nrow(seurat_data), "\n")
} else {
  cat("\nSpecies not recognized - proceeding with all genes...\n")
}

cat("Creating Seurat Object ...")
DATA <- CreateSeuratObject(seurat_data, 
                           project = opt$project_name, 
                           min.cells = opt$min_cells, 
                           min.features = opt$min_feats)
cat("done!\n\n")
print(DATA)

if (opt$sample_name != "none"){
  DATA@meta.data$sample <- opt$sample_name
  cat("\n\nCreating sample column in metadata.")
  file_save_name <- paste0(opt$project_name, "_", opt$sample_name)
} else {
  file_save_name <- opt$project_name
}

cat('\n\nSaving as RDS file ...')
saveRDS(DATA, file=paste0(opt$output_path, '/', file_save_name, '.rds'))
cat(paste0('\n\nSeurat Object saved as RDS file to ', paste0(opt$output_path, '/', file_save_name, '.rds')))

