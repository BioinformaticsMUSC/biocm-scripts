#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script begins VDJ analysis using scRepertoire.\n")
option_list = list(
  make_option(c("-v", "--vdj_input_dir"), type = "character",   metavar="character",   default='none',  help="Path for vdj input data."),
  make_option(c("-i", "--seurat_file_path"), type = "character",   metavar="character",   default='none',  help="Path for vdj input data."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default='seurat_vdj.rds',  help="nCount_RNA threshold."),
  make_option(c("-r", "--proportion"),      type = "character",   metavar="character",   default='false',  help="File path for pretrained classifier")
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
  library(patchwork)
  library(dplyr)
  library(scRepertoire)
  library(cowplot)
}))

#opt = data.frame(vdj_input_dir='/Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/input_data/vdj',
#                 seurat_file_path='/Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/seurat_processed_celltypes.rds',
#                 proportion="false")

files <- list.files(opt$vdj_input_dir, 
                    recursive = T, 
                    pattern = 'filtered_contig_annotations')

#files <- list.files(seurat_data_dir, pattern = 'rds')
if (length(files) != 0) {
  cat(length(files), "files found:\n\n")
  for (f in files) {
    cat(f,"\n")
  } 
} else {
  stop('\n\nNo filtered_contig_annotations files were found in the given directory.')
}
##############################
### LOAD/INSTALL LIBRARIES ###
##############################
seurat <- readRDS(opt$seurat_file_path)
cat('\nParsing contig files...\n')

contig_list <- list()
sample_names <- list()
sample_ids <- list()

for (i in 1:length(files)) {
  sample_name <- strsplit(files[i], "/")[[1]][1]
  sample_names[[i]] <- sample_name
  sample_id <- 's'
  sample_ids[[i]] <- sample_id
  contig_csv <- read.csv(file = paste0(opt$vdj_input_dir, "/", files[i]))
  contig_list[[i]] <- contig_csv
}

cat('\nSamples and ids found:\n\nSAMPLE\t\t|\tID\n')
for (j in 1:length(sample_names)){
  cat(sample_names[[j]], '\t\t|\t', sample_ids[[j]], '\n')
}
cat("\nNow creating combined contig list...")
combined <- combineBCR(contig_list, samples = sample_names, ID = sample_ids)
cat("done!\n")

cat("\nCombining with seurat object...")

#create raw_barcode field
seurat$raw_barcode <- rownames(seurat@meta.data)

#modify rownames for seurat obj metadata to match barcodes in combined
##HERE I AM HARD CODING THE ID - IF WORKING WITH MULTIPLE IDS, NEED TO MAKE A COLUMN FOR IT IN SEURAT IN UPSTREAM ANALYSIS
rownames(seurat@meta.data) <- paste0(seurat$sample, '_', sample_ids[[1]], '_', seurat$raw_barcode)
rownames(seurat@meta.data) <- gsub('_1', '', rownames(seurat@meta.data))
rownames(seurat@meta.data) <- gsub('_2', '', rownames(seurat@meta.data))

#MAY WANT TO HAVE CLONETYPES and CLONECALL AS OPTPARSE PARAMETERS
if (opt$proportion == "true") {
  seurat_combined <- combineExpression(combined, seurat, cloneCall = 'gene')
} else {
  seurat_combined <- combineExpression(combined, seurat, cloneCall = 'gene',
                                       proportion = FALSE, 
                                       cloneTypes = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
}
#revert meta.data rownames back so they match seurat_obj columns
rownames(seurat_combined@meta.data) <- seurat_combined$raw_barcode
cat('done!\n\nSaving new seurat object...')
#save seurat object
saveRDS(seurat_combined, file=paste0(opt$output_path, '/', opt$seurat_save_name))
saveRDS(combined, file=paste0(opt$output_path, "/combined_contig_list.rds"))

