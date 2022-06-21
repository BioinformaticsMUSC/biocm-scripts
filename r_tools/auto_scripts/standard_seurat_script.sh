#!/bin/bash

# Template for standard Seurat workflow

##########################
### DEFINE DIRECTORIES ###
##########################
# specify main project directory
main='/Users/bryanwgranger/biocm/projects/dubois'
r_tools='/Users/bryanwgranger/Documents/bioCM/r_tools'
project_name='dubois'

###########################
### PREPARE ENVIRONMENT ###
###########################

##############################
### CREATE LIST OF SAMPLES ###
##############################
#searches input_data directory for all samples
sample_list=(`ls -d $main/input_data/count/*/ | xargs -n 1 basename`)


#############################
### LOAD DATA INTO SEURAT ###
#############################

for sample in ${sample_list[@]}
do
	file=(`ls $main/input_data/count/$sample/*.h5`)
	Rscript $r_tools/01_load_seurat_data.R \
	-g $file \
	-p $project_name \
	-s $sample \
	-o /Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/out
done		

################################
### INTEGRATE SEURAT OBJECTS ###
################################

Rscript 02_integrate_seurats.R \
-i $main/out \
-p $project_name \
-o $main \
-m 20 \
-n 10000 \
-s seurat_integrated.rds

#########################
### STANDARD WORKFLOW ###
#########################

Rscript 03_standard_processing.R \
-i $main/seurat_integrated.rds \
-p $project_name \
-o $main \
-s seurat_processed.rds \
-t true \
-r none

#######################
### REMOVE DOUBLETS ###
#######################

RScript 03a_remove_doublets.R \
-i $main/seurat_processed.rds \
-p $project_name \
-o $main \
-s seurat_processed_filt.rds \
-t true \
-r 0.6

#################################
### FIND MARKERS FOR CLUSTERS ###
#################################

RScript 03b_find_cluster_markers.R \
-i $main/seurat_processed.rds \
-p $project_name \
-o $main \
-s seurat_processed_filt.rds \
-t true \
-r 0.6

###################################
### CLASSIFY CELLS WITH GARNETT ###
###################################

RScript 04_cell_classify_GARNETT.R \
-i /Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/seurat_processed_filt.rds \
-p scBCR \
-o /Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/markers \
-s seurat_processed_celltypes.rds \
-m /Users/bryanwgranger/Documents/bioCM/scBCRseq/hsPBMC_20191017.RDS \
-f /Users/bryanwgranger/Documents/bioCM/scBCRseq/hsPBMC_markers.txt \
-d human

