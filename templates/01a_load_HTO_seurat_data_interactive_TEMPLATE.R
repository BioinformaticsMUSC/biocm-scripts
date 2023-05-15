#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script loads Seurat Object from 10X data\n")

opt = list(data_dir = "{{data_dir}}",
           project_name = "{{project_name}}",
           output_path = "{{output_path}}",
           min_cells = {{min_cells}},
           min_feats = {{min_feats}},
           species = "{{species}}")

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


########################################
### FIND DATA FILES AND SAMPLE NAMES ###
########################################

f <- list.dirs(path = opt$data_dir)

ignore_dir = c('count', 'input_path', opt$project_name, 'vdj')

for (i in 1:length(f)){
  if (!(basename(f[i]) %in% ignore_dir)) {
    sample <- basename(f[i])
    cat('Working on ', sample)
    if (any(str_detect(list.files(f[i]), pattern = "cellbender_matrix_filtered.h5"))) {
      cat('\nReading from cellbender h5 file\n')
      cb_file = list.files(f[i])[str_which(list.files(f[i]), pattern = "cellbender_matrix_filtered.h5")]
      seurat_data <- Read10X_h5(paste0(f[i], "/", cb_file))
    } else if ("filtered_feature_bc_matrix.h5" %in% list.files(f[i])) {
      cat('\nReading from 10X file')
      seurat_data <- Read10X_h5(paste0(f[i], '/filtered_feature_bc_matrix.h5'))
    }

    #hto
    if (typeof(seurat_data) == "list") {
      umis <- seurat_data[['Gene Expression']]
      hto <- seurat_data[['Multiplexing Capture']]

      joint.bcs <- intersect(colnames(umis), colnames(hto))
  
      umis <- umis[,joint.bcs]
      hto <- hto[,joint.bcs]

    } else{
      umis <- seurat_data
    }
  cat("\nFiltering for Protein Coding Genes\n")
  
  if (opt$species == 'human') {
    cat('\nSpecies is human...\n')
    load('/Users/bryanwgranger/biocm/biocm-tools/r_tools/protein_coding_genes/HgProteinCodingGenes.rda')
    cat("Prefiltering genes total: ", nrow(umis))
    umis <- umis[rownames(umis) %in% HgProteinCodingGenes,]
    cat("\nPostfiltering genes total: ", nrow(umis), "\n")
  } else if (opt$species == 'mouse') {
    cat('\nSpecies is mouse...\n')
    load('/Users/bryanwgranger/biocm/biocm-tools/r_tools/protein_coding_genes/MgProteinCodingGenes.rda')
    cat("Prefiltering genes total: ", nrow(umis))
    umis <- umis[rownames(umis) %in% MgProteinCodingGenes,]
    cat("\nPostfiltering genes total: ", nrow(umis), "\n")
  } else {
    cat("\nSpecies not recognized - proceeding with all genes...\n")
  }
  
  cat("Creating Seurat Object ...")
  DATA <- CreateSeuratObject(umis, 
                             project = opt$project_name, 
                             min.cells = opt$min_cells, 
                             min.features = opt$min_feats)
  
  ##Add HTO
  new_bcs <- colnames(DATA)
  hto <- hto[,new_bcs]
  
  DATA[['HTO']] <- CreateAssayObject(counts = hto)
  DATA <- NormalizeData(DATA, assay = 'HTO', normalization.method = 'CLR')
  
  DATA <- RenameCells(DATA, add.cell.id = sample)
  DATA$sample <- sample
  file_save_name <- paste0(opt$project_name, "_", sample)
  cat("done!\n\n")
  print(DATA)
  
  saveRDS(DATA, file=paste0(opt$output_path, '/', file_save_name, '.rds'))
    cat(paste0('\n\nSeurat Object saved as RDS file to ', paste0(opt$output_path, '/', file_save_name, '.rds')))
    
  }
}
