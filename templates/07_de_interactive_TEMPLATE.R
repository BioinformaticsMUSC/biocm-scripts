#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis performs differential expression analysis using Libra.\n")

opt = list(seurat_file_path="{{seurat_file_path}}",
           project_name="{{project_name}}",
           output_path="{{output_path}}",
           cell_type_col='none',
           replicate_col='none',
           label_col='none',
           de_family="singlecell",
           de_method="MAST",
           de_type=NULL,
           split_comparisons='none' #add a Seurat meta.data column if you'd like to perform
                                      #additional DE analysis using subsets of seurat
           )
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
  library(Libra)

}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and preparing for Libra...')

seurat <- readRDS(opt$seurat_file_path)

#define function to run libra
run_libra <- function(seurat_obj, cell_type_col, replicate_col, label_col,
                      de_family, de_method, de_type, filename){
  DefaultAssay(seurat_obj) <- "SCT"
  tmp <- seurat_obj
  
  tmp@meta.data$cell_type <- tmp@meta.data[[opt$cell_type_col]]
  tmp@meta.data$replicate <- tmp@meta.data[[opt$replicate_col]]
  tmp@meta.data$label <- tmp@meta.data[[opt$label_col]]
  
  DE_output <- run_de(tmp, de_family = opt$de_family, de_method = opt$de_method,
                       de_type=opt$de_type)
  
  dge <- DE_output %>%
    mutate(Abs = abs(avg_logFC)) %>%
    filter(p_val_adj < 0.05 & Abs > 0.3) %>%
    arrange(desc(Abs))
  
  openxlsx::write.xlsx(dge, 
                       file = paste0(opt$output_path, "/", filename, ".xlsx"), 
                       colNames = TRUE,
                       rowNames = FALSE, 
                       borders = "columns",
                       sheetName="Stats",
                       overwrite=T)
  write.table(dge, paste0(opt$output_path, "/", filename, ".txt"),sep="\t",quote=F)
  
}

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and preparing for Libra...\n')
seurat <- readRDS(opt$seurat_file_path)

run_libra(seurat_obj = seurat,
          cell_type_col = opt$cell_type_col,
          replicate_col = opt$replicate_col,
          label_col = opt$label_col,
          de_family = opt$de_family,
          de_method = opt$de_method,
          de_type = opt$de_type,
          filename = "DE_overall")

cat('\nOverall Differential Expression analysis done!\n')

if (opt$split_comparisons != "none"){
  cat("Now working on subsetted Seurat objects")
  values <- unique(seurat@meta.data[[opt$split_comparisons]])
  for (val in values) {
    tmp_seurat <- seurat
    tmp_seurat$tmp_col <- seurat@meta.data[[split_comparisons]]
    tmp_seurat <- subset(tmp_seurat, tmp_col == val)
    
    run_libra(seurat_obj = tmp_seurat,
              cell_type_col = opt$cell_type_col,
              replicate_col = opt$replicate_col,
              label_col = opt$label_col,
              de_family = opt$de_family,
              de_method = opt$de_method,
              de_type = opt$de_type,
              filename = paste0("DE_", val))
  }
}
cat('done with differential expression')
saveRDS(seurat, file = paste0(opt$output_path, "/", opt$seurat_save_name))