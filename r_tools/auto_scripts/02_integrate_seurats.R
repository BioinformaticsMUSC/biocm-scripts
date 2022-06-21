#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nThis script integrates multiple Seurat objects and applies SCTransform to them\n")
option_list = list(
  make_option(c("-i", "--seurat_data_dir"), type = "character",   metavar="character",   default='none',  help="Path to directory with seurat objects."),
  make_option(c("-p", "--project_name"),     type = "character",   metavar="character",   default='seurat_project',   help="Overall project name."),
  make_option(c("-o", "--output_path"),      type = "character",   metavar="character",   default='none',  help="Output DIRECTORY."),
  make_option(c("-m", "--mito_thresh"),      type = "integer",   metavar="integer",   default=20,   help="Mitochondrial removal threshold."),
  make_option(c("-n", "--ncount_RNA_thresh"),      type = "integer",   metavar="integer",   default=10000,  help="nCount_RNA threshold."),
  make_option(c("-s", "--seurat_save_name"),      type = "character",   metavar="character",   default="integrated_seurat_object.rds",  help="Seurat save name - must have .rds at end.")
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
  library(dplyr)
  library(plyr)
  library(ggplot2)
  library(cowplot)
}))

####################################
### LOAD DATA INTO SEURAT OBJECT ###
####################################
cat("\nReading directory for RDS files ... ")
files <- list.files(opt$seurat_data_dir, pattern = 'rds')
#seurat_data_dir <- '/Users/bryanwgranger/biocm/projects/dubois/seurat_objects/raw_for_int
#files <- list.files(seurat_data_dir, pattern = 'rds')
if (length(files) != 0) {
  cat(length(files), " files found:\n\n")
  for (f in files) {
    cat(f,"\n")
  } 
  if (length(files) == 1) {
    stop("Only 1 file found. This script includes integration and expects multiple files.\nPlease use the script for single files.")
  }
} else {
    stop('\n\nNo rds files were found in the given directory.')
}

cat("\nLoading Seurat objects ...")
seurat_list <- list()
for (i in 1:length(files)) {
  seurat_obj <- readRDS(paste0(opt$seurat_data_dir, '/', files[i]))
  seurat_list[[i]] <-seurat_obj
}
cat('done!\n')

##merge
cat("\nMerging Seurat Objects ... ")
seurat <- seurat_list[[1]]
for (j in 1:(length(seurat_list)-1)){
  seurat <- merge(seurat, y=seurat_list[[j+1]])
}
cat("done!\n")

#percent.mt
seurat[['pMito']] <- PercentageFeatureSet(seurat, pattern = "^mt-")

if(sum(seurat[['pMito']]) == 0) {
  seurat[['pMito']] <- PercentageFeatureSet(seurat, pattern = "^MT-")
}

#save as unfiltered object
cat("\nSaving combined unfiltered seurat object ...")
saveRDS(seurat, file = paste0(opt$output_path, '/combined_unfiltered_seurat.rds'))
cat("done!\n")


###########################
### QUALITY CONTROL VIZ ###
###########################
cat("\nCreating Quality Control vizualizations ... ")
qcviz_path = paste0(opt$output_path, '/viz_qc')

if(!dir.exists(qcviz_path)){dir.create(qcviz_path,recursive = T)}

#QC plot
pdf(filename = paste0(qcviz_path, '/01_qc_plot_rna.png'))
VlnPlot(seurat, features = c('nCount_RNA', 'nFeature_RNA', 'pMito'),
        group.by = 'sample')
invisible(dev.off())

#QC plot
pdf(filename = paste0(qcviz_path, '/02_qc_plot_umi-gene.png'))
FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA", group.by = "sample", pt.size = 0)
invisible(dev.off())

pdf(filename = paste0(qcviz_path, '/03_qc_plot_numi-pmito.png'))
FeatureScatter(seurat, "nCount_RNA", "pMito", group.by = "sample", pt.size = 0)
invisible(dev.off())

#QC plot -figure out error with this
pdf(filename = paste0(qcviz_path, '/04_qc_plot_by_geno.png'))
ggplot(seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=sample)) +
  geom_point(size=0.5) 
#theme_classic()
invisible(dev.off())

# Visualize the correlation between genes detected and number of UMIs. D
# Determine whether strong presence of cells with low numbers of genes/UMIs
#png(filename = paste0(qcviz_path, '/05_qc_plot_correlation.png'))
scatter <- seurat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=pMito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
#invisible(dev.off())
ggsave(file = paste0(qcviz_path, "/05_qc_plot_correlation.pdf"))
cat

###############################################
## CREATE TABLE OF COUNTS TO BE FILTERED OUT ##
###############################################

#percent mito > 20
seurat@meta.data$mito_filter_out <- seurat@meta.data$pMito > opt$mito_thresh
filt <- as.data.frame(table(seurat@meta.data$mito_filter_out))
colnames(filt) <- c('Var1', paste0('percent_mito > ', opt$mito_thresh, '%'))

#nUMI >10000
seurat@meta.data$numi_out <- seurat@meta.data$nCount_RNA > opt$ncount_RNA_thresh
u_filt <- as.data.frame(table(seurat@meta.data$numi_out))
colnames(u_filt) <- c('Var1', paste0('nUMI > ', opt$ncount_RNA_thresh))

filtering <- merge(filt, u_filt, by='Var1')

#all cells to be filtered out
seurat@meta.data$all_out <- (seurat@meta.data$nCount_RNA > opt$ncount_RNA_thresh) | (seurat@meta.data$pMito > opt$mito_thresh)
all_filt <- as.data.frame(table(seurat@meta.data$all_out))
colnames(all_filt) <- c('Var1', 'Total')

filtering <-merge(filtering, all_filt, by ='Var1')

filtering$Var1 <- mapvalues(filtering$Var1, 
                            from = c('FALSE', 'TRUE'),
                            to = c('Remaining', 'To be filtered out'))
rownames(filtering) <- filtering$Var1
filtering <- filtering[,-1]

write.csv(filtering, file = paste0(opt$output_path, "/pm_numi_all_counts.csv"))

#####################
## Remove MT genes ##
#####################

cat('\nRemoving cells with high mitochondrial expression ...')
mito_filtered = seurat@assays$RNA@counts[-grep("^mt-", rownames(seurat@assays$RNA@counts)),]

seurat_final <- CreateSeuratObject(counts = mito_filtered, project = opt$project_name)

seurat_final <- AddMetaData(seurat_final, metadata = as.data.frame(seurat@meta.data))

seurat_filtered <- subset(seurat_final, 
                          subset = nCount_RNA < opt$ncount_RNA_thresh & pMito < opt$mito_thresh)

rm(seurat, seurat_final, mito_filtered)

df <- seurat_filtered@meta.data %>% as.data.frame()
tmp <- table(df$sample) %>% as.data.frame()

ggplot(tmp,aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle("Gene Expression by Sample") +
  xlab("Sample") +
  ylab("Gene Expression") +
  theme(legend.position="none")
ggsave(file = paste0(qcviz_path, "/06_qc_gex_bar_plot.pdf"))

cat('done!\n')

###################
### INTEGRATION ###
###################
cat('\nBeginning the integration process.\n')
#split object into list of seurat objects
seurat_split <- SplitObject(seurat_filtered, split.by = 'sample')

cat('\nPerforming SCTransform on each of the seurat objects ...\n')
#perform sctransform on each of the split seurat objects
for (i in 1:length(seurat_split)) {
  seurat_split[[i]] <- SCTransform(seurat_split[[i]], 
                                   vars.to.regress = c("nCount_RNA","pMito"), 
                                   verbose = T)
}
cat('\nSCTransform done!\n\nPreparing for integration.\n')
#find integration features
integ_features <- SelectIntegrationFeatures(object.list = seurat_split, nFeatures= 4000)

#run prep sc integration
seurat_split <- PrepSCTIntegration(object.list = seurat_split, 
                                   anchor.features = integ_features)

#find integration anchors
integ_anchors <- FindIntegrationAnchors(object.list = seurat_split, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

#integrate
cat('\nStarting to integrate data ...')
seurat_integrated <- IntegrateData(
  anchorset = integ_anchors,
  new.assay.name = "integrated",
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 100,
  sd.weight = 1,
  eps = 0.5,
  verbose = TRUE
)

DefaultAssay(seurat_integrated) <- 'integrated'

cat('done!\n')

saveRDS(seurat_integrated, file = paste0(opt$output_path, '/', opt$seurat_save_name))