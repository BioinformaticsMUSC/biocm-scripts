#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script integrates multiple Seurat objects and applies SCTransform to them\n")

opt = list(seurat_data_dir = "{{seurat_data_dir}}",
           project_name = "{{project_name}}",
           seurat_save_name = "{{seurat_save_name}}",
           seurat_output_path = "{{seurat_output_path}}",
           output_path = "{{output_path}}",
           mito_thresh = {{mito_thresh}},
           nCount_RNA_thresh= {{nCount_RNA_thresh}})

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
  library(stringr)
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
  cat('working on merge', j, '\n')
  seurat <- merge(seurat, y=seurat_list[[j+1]])
}
cat("done!\n")

#percent.mt
seurat[['pMito']] <- PercentageFeatureSet(seurat, pattern = "^[Mm][Tt]-")

if(sum(seurat[['pMito']]) == 0) {
  seurat[['pMito']] <- PercentageFeatureSet(seurat, pattern = "^MT-")
}

##UPDATE METADATA #### UPDATE THIS AS NEEDED
seurat@meta.data <- seurat@meta.data %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) 

#save as unfiltered object
cat("\nSaving combined unfiltered seurat object ...")
saveRDS(seurat, file = paste0(opt$seurat_output_path, '/combined_unfiltered_seurat.rds'))
cat("done!\n")


###########################
### QUALITY CONTROL VIZ ###
###########################
cat("\nCreating Quality Control vizualizations ... ")
qcviz_path = paste0(opt$output_path, '/viz_qc')

if(!dir.exists(qcviz_path)){dir.create(qcviz_path,recursive = T)}

#QC plot
pdf(file = paste0(qcviz_path, '/01_qc_plot_rna.pdf'), width = 15)
VlnPlot(seurat, features = c('nUMI', 'nGene', 'pMito'),
        group.by = 'sample')
invisible(dev.off())

pdf(file = paste0(qcviz_path, '/01a_qc_plot_rna.pdf'), width = 15)
VlnPlot(seurat, features = c('nUMI', 'nGene', 'pMito'), pt.size = 0,
        group.by = 'sample') + theme(legend.position = 'none')
invisible(dev.off())

#density plot
pdf(file = paste0(qcviz_path, '/01b_density_plot_rna.pdf'))
seurat@meta.data %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density")
invisible(dev.off())

#QC plot
pdf(file = paste0(qcviz_path, '/02_qc_plot_umi-gene.pdf'))
FeatureScatter(seurat, "nUMI", "nGene", group.by = "sample", pt.size = 0)
invisible(dev.off())

pdf(file = paste0(qcviz_path, '/03_qc_plot_numi-pmito.pdf'))
FeatureScatter(seurat, "nUMI", "pMito", group.by = "sample", pt.size = 0)
invisible(dev.off())

#QC plot -figure out error with this
pdf(file = paste0(qcviz_path, '/04_qc_plot_by_geno.pdf'))
ggplot(seurat@meta.data, aes(x=nUMI, y=nGene, color=sample)) +
  geom_point(size=0.5) 
#theme_classic()
invisible(dev.off())

# Visualize the correlation between genes detected and number of UMIs. D
# Determine whether strong presence of cells with low numbers of genes/UMIs
#png(filename = paste0(qcviz_path, '/05_qc_plot_correlation.png'))
scatter <- seurat@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=pMito)) + 
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

## Create table of counts for each sample
for (samp in unique(seurat$sample)){
  tmp <- as.data.frame(seurat@meta.data) %>%
    filter(sample == samp)
  
  #percent mito > 20
  tmp <- tmp %>% mutate(pmito_out = pMito > opt$mito_thresh)
  filt <- as.data.frame(table(tmp$pmito_out))
  colnames(filt) <- c('Var1', paste0('percent_mito > ', opt$mito_thresh, '%'))
  
  #nUMI >10000
  tmp <- tmp %>% mutate(numi_out = tmp$nUMI > opt$nCount_RNA_thresh)
  u_filt <- as.data.frame(table(tmp$numi_out))
  colnames(u_filt) <- c('Var1', paste0('nUMI > ', opt$ncount_RNA_thresh))
  
  filtering <- merge(filt, u_filt, by='Var1')
  
  #all cells to be filtered out
  tmp <- tmp %>% mutate(all_out = nUMI > opt$nCount_RNA_thresh | pMito > opt$mito_thresh)
  all_filt <- as.data.frame(table(tmp$all_out))
  colnames(all_filt) <- c('Var1', 'Total')
  
  filtering <-merge(filtering, all_filt, by ='Var1')
  
  filtering$Var1 <- plyr::mapvalues(filtering$Var1, 
                              from = c('FALSE', 'TRUE'),
                              to = c('Remaining', 'To be filtered out'))
  rownames(filtering) <- filtering$Var1
  filtering <- filtering[,-1]
  
  write.csv(filtering, file = paste0(qcviz_path, "/pm_numi_all_counts_", samp, ".csv"))
}

## Create table of counts for all combined

#percent mito > 20
seurat@meta.data <- seurat@meta.data %>% mutate(mito_filter_out = pMito > opt$mito_thresh)
filt <- as.data.frame(table(seurat@meta.data$mito_filter_out))
colnames(filt) <- c('Var1', paste0('percent_mito > ', opt$mito_thresh, '%'))

#nUMI >10000
seurat@meta.data <- seurat@meta.data %>% mutate(numi_out = nUMI > opt$nCount_RNA_thresh)
u_filt <- as.data.frame(table(seurat@meta.data$numi_out))
colnames(u_filt) <- c('Var1', paste0('nUMI > ', opt$ncount_RNA_thresh))

filtering <- merge(filt, u_filt, by='Var1')

#all cells to be filtered out
seurat@meta.data <- seurat@meta.data %>% 
  mutate(all_out = (nUMI > opt$nCount_RNA_thresh) | (pMito > opt$mito_thresh))
all_filt <- as.data.frame(table(seurat@meta.data$all_out))
colnames(all_filt) <- c('Var1', 'Total')

filtering <-merge(filtering, all_filt, by ='Var1')

filtering$Var1 <- mapvalues(filtering$Var1, 
                            from = c('FALSE', 'TRUE'),
                            to = c('Remaining', 'To be filtered out'))
rownames(filtering) <- filtering$Var1
filtering <- filtering[,-1]

write.csv(filtering, file = paste0(qcviz_path, "/pm_numi_all_counts.csv"))

#####################
## Remove MT genes ##
#####################

nUMI_thresh <- opt$nCount_RNA_thresh
mito_thresh <- opt$mito_thresh
cat('\nRemoving cells with high mitochondrial expression ...')

mito_filtered = seurat@assays$RNA@counts[str_subset(rownames(seurat@assays$RNA@counts),
                                                    pattern = "^[Mm][Tt]",
                                                    negate = TRUE),]


seurat_final <- CreateSeuratObject(counts = mito_filtered, project = opt$project_name)

seurat_final <- AddMetaData(seurat_final, metadata = as.data.frame(seurat@meta.data))

seurat_filtered <- subset(seurat_final, 
                          subset = nUMI < nUMI_thresh & pMito < mito_thresh)

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

saveRDS(seurat_filtered, file = paste0(opt$seurat_output_path, "/filtered_seurat_pre_int.rds"))

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
                                   vars.to.regress = c("nUMI","pMito"), 
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

saveRDS(seurat_integrated, file = paste0(opt$seurat_output_path, '/', opt$seurat_save_name))
