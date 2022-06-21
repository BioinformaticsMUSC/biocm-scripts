#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis performs analysis and visualization of differential expression - only to be run after 07_differential_expression.\n")

opt = list(seurat_file_path='/Users/bryanwgranger/biocm/projects/dubois/seurat_objects/07_dubois_de.rds',
           project_name='dubois',
           output_path='/Users/bryanwgranger/biocm/projects/dubois/main_analysis/markers',
           seurat_save_name='07_dubois_de.rds',
           de_table_file='none',
           replicate_col='tissue',
           plot_title="sample vs sample")

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
           
#new integrated DGE
##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(ggrepel)
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(cowplot)
  library(EnhancedVolcano)
  library(patchwork)
  library(dplyr)
  library(plyr)
}))

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

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat('\nLoading seurat object and preparing for Libra...')

seurat <- readRDS(opt$seurat_file_path)
DefaultAssay(seurat) <- "SCT"

viz_path <- paste0(opt$output_path, '/viz')
if(!dir.exists(viz_path)){dir.create(viz_path,recursive = T)}


#############
##ANALYSIS###
#############

de_table <- read.table(file = opt$de_table_file)

col <- randomcoloR::distinctColorPalette(21)

dge <- de_table %>%
  mutate(Abs = abs(avg_logFC)) %>%
  filter(p_val_adj < 0.05 & Abs > 0.3) %>%
  arrange(desc(Abs))

df_celltype <- dge %>%
  mutate(LOG = -log10(p_val_adj)) %>%
  mutate(Threshold = if_else(p_val_adj < 0.05 & Abs > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(avg_logFC > 0.3 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.3 & p_val_adj < 0.05 ~ "DownReg"))

ggscatter(df_celltype, 
          x = "avg_logFC", 
          y = "LOG",
          color = "cell_type",
          palette=col,
          size = 1,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+ 
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) + 
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = df_celltype, 
                  mapping = aes(label = gene), 
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,75) + xlim(-5,+5)
ggsave(filename = paste0(viz_path, "/volcano_plot.pdf"), height = 12, width = 12)

EnhancedVolcano(df_celltype,
                lab = df_celltype$gene,
                x = "avg_logFC", 
                y = "p_val_adj",
                title = opt$plot_title,
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                col=c('grey', 'black', 'blue', 'red3'),
                colAlpha = 1,
                legendLabels=c('NotSig','Log2FC','p-val',
                               'p-val&Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 5.0) +
  ylim(0,300) + xlim(-5,+5)
ggsave(filename = paste0(viz_path, "/volcano_plot_enhanced.pdf"), height = 12, width = 12)

## Volcano plots for each cell type
for (cell in unique(df_celltype$cell_type)) {
  new_df <- subset(df_celltype, cell_type == cell)
  max_y <- max(-log10(filter(new_df, p_val_adj != 0)$p_val_adj))
  EnhancedVolcano(new_df,
                  lab = new_df$gene,
                  x = "avg_logFC", 
                  y = "p_val_adj",
                  title = paste0(opt$output_path, cell),
                  pCutoff = 10e-32,
                  FCcutoff = 0.5,
                  pointSize = 3.0,
                  labSize = 6.0,
                  col=c('grey', 'black', 'blue', 'red3'),
                  colAlpha = 1,
                  legendLabels=c('NotSig','Log2FC','p-val',
                                 'p-val&Log2FC'),
                  legendPosition = 'bottom',
                  legendLabSize = 16,
                  legendIconSize = 5.0) +
    ylim(0,max_y+20) + xlim(-5,+5)
  ggsave(filename = paste0(viz_path, "/volcano_plot_enhanced_",cell,".pdf"), height = 12, width = 12)
}

##Violin plots for each cell type
##Note: this subsets the seurat for each cell type and then plots from that
for (cell in unique(dge$cell_type)){
  top10 <- dge %>%
    filter(cell_type == cell) %>%
    arrange(-Abs) %>%
    head(10)

  if (nrow(top10) >= 3){
  Idents(seurat) <- "sample"
  tmp <- subset(seurat, Cell == cell)
  p <- VlnPlot(tmp, features=top10$gene[1:2], cols = biocm_colors(length(unique(seurat$sample))), combine = FALSE)
  (p[[1]] + theme(legend.position = "none") +
      p[[2]] + theme(legend.position = "none")) +
      #p[[3]] + theme(legend.position = "none")) + 
    plot_annotation(title = paste0(cell, " - top 2 differentially expressed genes"),
    caption = paste0("N = ", nrow(filter(dge, cell_type == cell))))
  ggsave(filename = paste0(viz_path, "/sub_Vln_", cell, ".pdf"), width = 12)

  #Idents(seurat) <- "HTO_maxID"
  #VlnPlot(seurat, features=top10$gene[1:3])
  #ggsave(filename = paste0(save_dir, "de/viz/hto/Ridge_", cell, "_hto.pdf"), width = 12)

  }
}
#DotPlot- top 20 (or fewer) features,
Idents(seurat) <- "Cell"
overall_top <- dge %>%
  arrange(-Abs) %>%
  head(20)

for (rep in unique(seurat[[opt$replicate_col]])){
  tmp <- subset(seurat, opt$replicate_col == rep)
  DotPlot(tmp, features = unique(overall_top$gene), assay = 'SCT', split.by = 'condition') + rotate_x_text(angle = 45) +
    theme(legend.position = "right") + ggtitle(rep)
  ggsave(paste0(viz_path, "/dotplot_topfeats_", rep, ".pdf"), width = 8)
}
DotPlot(seurat, features = unique(overall_top$gene), assay = 'SCT', split.by = 'condition') + rotate_x_text(angle = 45) +
  theme(legend.position = "right") + ggtitle('Overall')
ggsave(paste0(viz_path, "/dotplot_topfeats.pdf"), width = 8)

cat('Done with visualizations.')


