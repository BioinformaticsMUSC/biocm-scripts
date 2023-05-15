#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script predicts cell types using Garnett.\n")


opt = list(seurat_file_path = "{{seurat_file_path}}",
           project_name = "{{project_name}}",
           seurat_save_name = "{{seurat_save_name}}",
           seurat_output_path = "{{seurat_output_path}}",
           output_path = "{{output_path}}",
           model_file_path = "none")

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
  library(SingleCellExperiment)
  library(scPred)
  library(tidyr)
  library(tibble)
  library(ggpubr)
}))


####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and converting to CellDataSet...')
seurat <- readRDS(opt$seurat_file_path)
data <- readRDS(opt$model_file_path)
DefaultAssay(data) <- "integrated"

#create new Seurat obj from reference if needed - maybe if from tabula muris

#reference <- CreateSeuratObject(tiss@data, meta.data = tiss@meta.data)

#downsamples <- subset(seurat, downsample = 300) # I downsample to 300 cells just to make it easier. You can downsample to 1000 or 3000.
seurat@meta.data$Data <- opt$project_name
data@meta.data$Data <- "reference"
#reference@meta.data$Cell <- reference@meta.data$orig.ident

### only if needed
process <- FALSE


if (isTRUE(process)) {
  data <- data %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000, verbose = FALSE) %>%
    ScaleData() %>%
    RunPCA(features=NULL, weight.by.var = TRUE, npcs = 30, reduction.name = "pca") 
}

#check for NA values in cell column and subset
data <- subset(reference, !(data = Cell %>% is.na()))

reference <- getFeatureSpace(data, "Cell")

reference <- trainModel(reference)

query <- scPredict(seurat, reference, max.iter.harmony = 30)

DimPlot(query, group.by = 'scpred_prediction', label = T, repel = T)
ggsave(file = paste0(opt$output_path, '/dimplot_celltype_scPred.pdf'))

query$Cell <- query$scpred_prediction

mat <- crossTab(query, "integrated_snn_res.0.6", "scpred_prediction", output = "prop") %>%
  as.data.frame() %>%
  rownames_to_column("Cell") %>%
  pivot_longer(!Cell, names_to = "Cluster", values_to = "Pred") %>%
  group_by(Cluster) %>%
  filter(!(Cell %in% "unassigned")) %>%
  filter(Pred == max(Pred)) %>%
  arrange(Cluster) %>%
  as.data.frame()

ggbarplot(mat, x = "Cluster", y = "Pred",
          fill = "Cluster",           # change fill color by mpg_level
          color = "black",            # Set bar border colors to white
          palette = mat$Colors,            # jco journal color palett. see ?ggpar
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 45,          # Rotate vertically x axis texts
          ylab = "Deleon, Predicted %",
          rotate = TRUE,
          ggtheme = theme_classic()
) +
  theme(legend.position="none")

write.table(mat, file = paste0(opt$output_path, "/scpred_cluster_table.txt"),
            sep = '\t', quote = F, row.names = F)

saveRDS(query, file = paste0(opt$seurat_save_path, "/", opt$seurat_save_name))
