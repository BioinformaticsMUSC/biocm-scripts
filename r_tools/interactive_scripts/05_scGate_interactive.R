#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script uses scGate to check purity of cell types.\n")


opt = list(seurat_file_path='/Users/bryanwgranger/biocm/projects/dubois/seurat_objects/seurat_processed_filt.rds',
           project_name='dubois',
           output_path='/Users/bryanwgranger/biocm/projects/dubois/main_analysis/markers',
           seurat_save_name='seurat_custom_garnett.rds',
           species='mouse'
           )

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

scgate_path <- paste0(opt$output_path, "/scGate")

if(!dir.exists(scgate_path)){dir.create(scgate_path,recursive = T)}
setwd(scgate_path)

##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(scGate)
}))


####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and scGate models...')
seurat <- readRDS(opt$seurat_file_path)

models <- get_scGateDB()

cat('done!\n')

##############################
### CREATE SCGATE FUNCTION ###
##############################

plot_scGate <- function(seurat, scGate_model, save_dir){
  seurat <- scGate(seurat, model=scGate_model, assay = "RNA")
  DimPlot(seurat, cols = c(list(Impure = "gray", Pure = "green"))) + theme(aspect.ratio = 1) + 
    ggtitle(paste0(model_name, " scGate model"))
  ggsave(filename = paste0(save_dir, "/", model_name, '_scgate_dimplot.pdf'))
  
}

######################################
### CREATE SCGATE MODELS and LISTS ###
######################################


bcell_model <- gating_model(name = "Bcell", signature = c('Cd19', 'Ms4a1', 'Cd79a'))
endothelial_model <- gating_model(name = "Endothelial", signature = c('Lyve1', 'Tek', 'Kdr', 'Ramp2', 'Flt1'))
fib_model <- gating_model(name = 'Fibroblast', signature = c('Des', 'Pdgfra', 'Pdgfrb', 'Fap'))
gran_model <- gating_model(name = "Granulocyte", signature = c('Cd11b', 'Ly6g', 'Ly6b', 'Ly6c', 'Cd193', 'Siglec-f'))
epithelial_model <- gating_model(name = "Epithelial", signature = c('Fabp4', 'Guca2a', 'Krt8', 'Gpx2', 'Sult1a1'))
mesenchymal_model <- gating_model(name = 'Mesenchymal', signature = c('Ltbp', 'Tcf21', 'Sod3', 'Bmp4'))

new_models <- list(bcell_model, endothelial_model, fib_model, gran_model, epithelial_model, mesenchymal_model)


#################################
### CREATE SCGATE MODEL PLOTS ###
#################################
for (m in names(models[[opt$species]][['generic']])) {
  cat("\nWorking on model: ", m)
  plot_scGate(seurat, scGate_model = models[[opt$species]][['generic']][[m]], model_name = m, save_dir = scgate_path)
}

for (model in names(new_models)) {
  cat("\nWorking on model: ", model)
  plot_scGate(seurat, scGate_model = new_models[[model]], model_name = model, save_dir = scgate_path)
}

saveRDS(seurat, file = paste0(opt$output_path, '/', opt$seurat_save_name))

