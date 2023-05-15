#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script creates visualizations.\n")

opt = list(seurat_file_path = "{{seurat_file_path}}",
           project_name = "{{project_name}}",
           output_path = "{{output_path}}",
           seurat_save_name = "{{seurat_save_name}}",
           seurat_output_path = "{{seurat_output_path}}",
           label_marker_file = "none"
          )

print(t(t(unlist(opt))))

biocm_colors <- function(n){
  if (n < 9) {
    
    palettes <- list("1" = c('#003f5c'),
                     "2" = c('#003f5c', '#ffa600'),
                     "3" = c('#003f5c', '#bc5090', '#ffa600'),
                     "4" = c("#003f5c", "#7a5195", "#ef5675", "#ffa600"),
                     "5" = c('#003f5c', '#58508d', '#bc5090', '#ff6361', '#ffa600'),
                     "6" = c("#003f5c", "#444e86", "#955196", "#dd5182", "#ff6e54", "#ffa600"),
                     "7" = c("#003f5c","#374c80","#7a5195","#bc5090","#ef5675","#ff764a","#ffa600"),
                     "8" = c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a", "#ff7c43", "#ffa600"))
    
    return(palettes[n])
  } 
}


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggpubr)
  library(plyr)
  library(broom)
  library(muscat)
  library(Seurat)
  library(leiden)
  library(cowplot)
  library(presto)
  library(ComplexHeatmap)
})

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object...')
seurat <- readRDS(opt$seurat_file_path)

labels <- read.table(opt$label_marker_file, header = T, sep = '\t')
cat('done!\n')

#########################
### UPDATE CELL TYPES ###
#########################
cat("\nUpdating cell types according to label file...")
Idents(seurat) <- 'seurat_clusters'

current.cluster.ids <- labels$seurat_clusters
new.cluster.ids <- as.character(labels$Cell)

seurat@active.ident <- mapvalues(seurat@active.ident, from = current.cluster.ids, to = new.cluster.ids)

seurat@meta.data$Cell <- seurat@active.ident
cat("done!\n")


#############################
### CREATE VISUALIZATIONS ###
#############################

samples <- unique(seurat$sample)
#Overall dimplot by cluster
DimPlot(seurat, group.by = 'seurat_clusters')
ggsave(filename = paste0(opt$output_path, "/OVERALL_dimplot.pdf"), width = 8, height = 8)

#DimPlots split by sample
if (length(unique(seurat$sample)) == 2 | length(unique(seurat$sample)) == 3){
  DimPlot(seurat, group.by = "Cell", split.by = "sample")
  ggsave(filename = paste0(opt$output_path, "/SPLIT_dimplot.pdf"), width = 12, height = 8)
} else if (length(unique(seurat$sample)) == 4) {
  p_list <- list()
  for (i in 1:length(samples)) {
    p_list[i] <- DimPlot(subset(seurat, sample == samples[i])) + ggtitle(samples[i])
  }
  ((p_list[1] + p_list[2]) / (p_list[3] + p_list[4]))
  ggsave(filename = paste0(opt$output_path, "/SPLIT_dimplot.pdf"), width = 8, height = 8)
} else {
  for (s in unique(seurat$sample)){
    DimPlot(subset(seurat, sample == s))
    ggsave(filename = paste0(opt$output_path, "/SPLIT_dimplot_", s, ".pdf"), height = 8, width = 8)
  }
}

#CLEAN UMAP with cell types
umap <- as.data.frame(Embeddings(seurat, reduction = "umap"))
meta <- as.data.frame(seurat@meta.data)

df <- cbind(umap,meta)%>% 
  group_by(Cell) %>% 
  mutate(N = n()) %>%
  ungroup() %>% 
  arrange(sample) %>%
  as.data.frame()

label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  group_by(Cell) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame() 

col <- randomcoloR::distinctColorPalette(length(unique(seurat@meta.data$Cell)))

write.table(col, paste0(opt$output_path, "/Colors_Used_Initial_Data.txt"), sep="\t",quote=F)

ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(aes(colour = Cell),size=0.5) +
  ggrepel::geom_text_repel(data = label_2, aes(label = label),
                           color = "black",
                           #fontface = 'bold',
                           segment.colour = "grey60",
                           box.padding = unit(0.25, "lines"),
                           point.padding = unit(0.5, "lines"),
                           nudge_x = .15,
                           nudge_y = 1,
                           size = 6) + 
  #scale_color_viridis(discrete=TRUE,option="inferno")
  #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
  #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
  scale_colour_manual(values = col)+
  theme_void() +
  theme(legend.position="none") + 
  ggtitle('UMAP plot with cell types')
ggsave(filename = paste0(opt$output_path, '/DIMPLOT_clean_clusters.pdf'))

seurat$cell_color <- plyr::mapvalues(seurat$Cell, from = label_2$Cell, to = col)

cat("\nSaving Seurat object with new cell types...")
saveRDS(seurat, paste0(opt$seurat_output_path, "/", opt$seurat_save_name))
cat('done!\n')

#barplot of sample cell proportions

plot_cell_pct <- function(seurat, split_col) {
  splits <- unique(seurat@meta.data[[split_col]])
  data <- seurat@meta.data %>%
    as.data.frame() %>%
    mutate(new_col = seurat@meta.data[[split_col]])
  
  mat <- data.frame()
  
  for (i in 1:length(splits)){
    tmp <- data %>%
      filter(new_col == splits[i])
    ncells <- ncol(tmp)
    t_mat <- as.data.frame(table(tmp$Cell))
    t_mat$overall_freq <- as.data.frame(table(data$Cell))$Freq
    t_mat$pct <- t_mat$Freq / t_mat$overall_freq
    t_mat$sample <- splits[i]
    colnames(t_mat) <- c('Cell', 'Freq', 'Overall Freq', 'Pct', 'Sample')
    
    mat <- rbind(mat, t_mat)
  }
  
  bp_pal <- randomcoloR::distinctColorPalette(length(splits))
  ggplot(mat,
         aes(x = Cell,
             y = Pct,
             fill = Sample,)) +
    geom_bar(position = 'fill', stat = 'identity', ) +
    rotate_x_text(angle = 45) +
    scale_fill_manual(values = bp_pal) + 
    ggtitle('Percentage of cells per sample') +
    geom_hline(yintercept = 0.5, colour = "grey", linetype = "dashed", size=0.5, alpha=0.5) 
  ggsave(filename = paste0(opt$output_path, "/barplot_cell_pct_", split_col, ".pdf"))
  
  
}

plot_cell_pct(seurat, "sample")

cat('\nDone with visualizations.\n')


######################
### MARKER HEATMAP ###
######################
markers<- wilcoxauc(seurat, 'Cell', assay = 'data')
write.csv(markers, file = paste0(opt$output_path, "/markers_all_celltype.csv"))
t_markers<- top_markers(markers, n = 20, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
write.csv(t_markers, file = paste0(opt$output_path, "/top_markers_celltype.csv"))

all_markers<- t_markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]
cols <- viridis::viridis(100)[c(1, 50, 100)]

mat <- seurat[["RNA"]]@data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))

cluster_anno<- seurat@meta.data$Cell

col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#440154FF", "#2D708EFF", "#FDE725FF"))

pdf(file = paste0(opt$output_path, '/complex_heat_map.pdf'), height = 12, width = 12)
ComplexHeatmap::Heatmap(mat, name = "Expression",  
                        column_split = cluster_anno,
                        cluster_columns = FALSE,
                        show_column_dend = FALSE,
                        cluster_column_slices = TRUE,
                        column_title_gp = grid::gpar(fontsize = 8),
                        column_gap = unit(0.5, "mm"),
                        cluster_rows = FALSE,
                        show_row_dend = FALSE,
                        col = col_fun,
                        row_names_gp = grid::gpar(fontsize = 4),
                        column_title_rot = 90,
                        top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_block(gp = grid::gpar(fill = col))),
                        show_column_names = FALSE,
                        use_raster = TRUE,
                        show_row_names = FALSE,
                        raster_quality = 4)
dev.off()
