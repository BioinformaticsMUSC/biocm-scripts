#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################

cat("\nThis script runs ArchR basic analysis.\n")

#set the cr_outs_path to the directory that contains all the samples (each of which has the 'outs' directory within)

opt = list(cr_outs_path = "/mnt/bg/asd/cellranger_out", 
           project_name = "asd_multiome",
           output_path = "/mnt/bg/asd/archr",
           genome = "hg38")

print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
suppressMessages(suppressWarnings({
  library(ArchR)
  library(stringr)
  library(patchwork)
  library(ggplot2)
  library(stringr)
}))

############################
### CREATE ARCHR PROJECT ###
############################

inputFiles = list.files(opt$cr_run_path, recursive=T, pattern='atac_fragments.tsv.gz$')
sampleNames = str_split_i(inputFiles, pattern="/", 1)

addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = opt$output_path,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

#### PLOT UMAPS ####
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
wrap_plots(p1, p2)
ggsave('archr/joint_umaps.pdf', width = 10)

#### PLOT QC ####
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "01_QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)

plotPDF(p1,p2, name = "02_QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- saveArchRProject(ArchRProj = proj, load = FALSE)

#### CREATE MARKER GENE SCORES ####
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  transpose = TRUE
)

plotPDF(heatmapGS, name = "03_GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

cat('\nSaving marker gene scores...')
saveRDS(markersGS, paste0(opt$output_path, "/markersGS.rds"))

#### call peaks ####

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

proj <- addReproduciblePeakSet(proj)

proj <- addPeakMatrix(proj)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "04_Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

saveRDS(markersPeaks, paste0(opt$output_path, "/markersPeaks.rds"))

proj <- saveArchRProject(ArchRProj = proj, load = FALSE)


