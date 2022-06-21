#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script begins VDJ analysis using scRepertoire.\n")
option_list = list(
  make_option(c("-c", "--combined_contig_file"), type = "character",   metavar="character",   default='none',  help="Path for combined contig file."),
  make_option(c("-i", "--seurat_file_path"), type = "character",   metavar="character",   default='none',  help="Path for vdj input data."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default='seurat_vdj_highlights.rds',  help="nCount_RNA threshold.")
)

opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

viz_dir <- paste0(opt$output_path, "/viz")
if(!dir.exists(viz_dir)){dir.create(viz_dir, recursive = T)}
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


#################
### LOAD DATA ###
#################
cat('\nLoading data ...')

seurat <- readRDS(opt$seurat_file_path)
combined <- readRDS(opt$combined_contig_file)
cat('done!\n')

################
### VIZ TIME ###
################

cat('\nCreating visualizations...\n')
#contig analysis
#quantContig - number of unique clonotypes per sample
png(filename = paste0(viz_dir, '/1_quantContig_percent.png'), width = 800, height = 600)
quantContig(combined, cloneCall="gene+nt", scale = T)
invisible(dev.off())

png(filename = paste0(viz_dir, '/2_quantContig_total.png'), width = 800, height = 600)
quantContig(combined, cloneCall="gene+nt", scale = F)
invisible(dev.off())

#abundanceContig
png(filename = paste0(viz_dir, '/3_abundanceContig.png'), width = 800, height = 600)
abundanceContig(combined)
invisible(dev.off())

#length of the clonotypes
png(filename = paste0(viz_dir, '/4_lengthContig_aa.png'), width = 800, height = 600)
lengthContig(combined, cloneCall = 'aa')
invisible(dev.off())

png(filename = paste0(viz_dir, '/4_lengthContig_nt.png'), width = 800, height = 600)
lengthContig(combined, cloneCall = 'nt')
invisible(dev.off())

#clonalHomeostasis
png(filename = paste0(viz_dir, '/5_clonalHomeostasis_gene-nt.png'), width = 800, height = 600)
clonalHomeostasis(combined, cloneCall = 'gene+nt')
invisible(dev.off())

png(filename = paste0(viz_dir, '/6_clonalHomeostasis_aa.png'), width = 800, height = 600)
clonalHomeostasis(combined, cloneCall = 'aa')
invisible(dev.off())

#compareClonotypes
png(filename = paste0(viz_dir, '/7_clonalOverlap.png'), width = 800, height = 600)
compareClonotypes(combined, cloneCall = 'gene+nt', numbers = 10, samples = c('pbmc_s', 'tumor_s'))
invisible(dev.off())

pbmc <- filter(seurat@meta.data, sample == 'pbmc')
tumor <- filter(seurat@meta.data, sample == 'tumor')
overlaps <- intersect(tumor[,'CTstrict'], pbmc[,'CTstrict'])
overlaps <- overlaps[!is.na(overlaps)]

png(filename = paste0(viz_dir, '/7a_clonalOverlap_alloverlaps.png'), width = 800, height = 600)
compareClonotypes(combined, cloneCall = 'gene+nt', clonotypes = overlaps, samples = c('pbmc_s', 'tumor_s'))
invisible(dev.off())


#clonal proportion
png(filename = paste0(viz_dir, '/8_clonalProportion_gene-nt.png'), width = 800, height = 600)
clonalProportion(combined, cloneCall = "gene+nt", split = c(10, 100, 1000, 10000, 30000, 1e+05))
invisible(dev.off())

png(filename = paste0(viz_dir, '/9_clonalProportion_aa.png'), width = 800, height = 600)
clonalProportion(combined, cloneCall = "aa", split = c(10, 100, 1000, 10000, 30000, 1e+05))
invisible(dev.off())

#single cell analysis
png(filename = paste0(viz_dir, '/10_dimplot_clonetypes.png'), width = 800, height = 600)
DimPlot(seurat, group.by = 'cloneType')
invisible(dev.off())

png(filename = paste0(viz_dir, '/11_dimplot_clonetypes_split.png'), width = 800, height = 600)
DimPlot(seurat, group.by = 'cloneType', split.by = 'sample')
invisible(dev.off())

#highlight top 10 clonotypes
n <- 10
top_highlights <- head(sort(table(seurat@meta.data$CTstrict), decreasing = TRUE), n)

seurat <- highlightClonotypes(seurat, cloneCall = 'gene+nt', sequence = rownames(top_highlights))

png(filename = paste0(viz_dir, '/12_clonetype_highlights.png'))
DimPlot(seurat, group.by = 'highlight')
invisible(dev.off())

png(filename = paste0(viz_dir, '/13_clonetype_highlights_split.png'))
DimPlot(seurat, group.by = 'highlight', split.by = 'sample')
invisible(dev.off())

#occupiedscrep
png(filename = paste0(viz_dir, '/14_sc_clonotype_by_cluster.png'))
occupiedscRepertoire(seurat, x.axis = "cluster")
invisible(dev.off())

png(filename = paste0(viz_dir, '/15_sc_clonotype_by_celltype.png'))
if ("Cell" %in% colnames(seurat@meta.data)) {
  
  occupiedscRepertoire(seurat, x.axis ='Cell')
  
}
invisible(dev.off())

cat('\nsaving new seurat object...')
#save seurat obj
saveRDS(seurat, file = paste0(opt$output_path, "/", opt$seurat_save_name))