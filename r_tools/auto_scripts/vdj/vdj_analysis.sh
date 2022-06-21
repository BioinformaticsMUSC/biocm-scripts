#!/bin/bash

# Template for standard Seurat/VDJ workflow

##########################
### DEFINE DIRECTORIES ###
##########################
# specify main project directory
main='/Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/'
vdj='/Users/bryanwgranger/Documents/bioCM/scBCRseq/new_analysis/vdj'
r_tools='/Users/bryanwgranger/Documents/bioCM/r_tools'
project_name='scBCR'

###########################
### PREPARE ENVIRONMENT ###
###########################

##############################
### CREATE LIST OF SAMPLES ###
##############################
#searches input_data directory for all samples
sample_list=(`ls -d $main/input_data/vdj/*/ | xargs -n 1 basename`)

#########################################
### MERGE VDJ DATA WITH SEURAT OBJECT ###
#########################################

Rscript vdj01_prep_scRep.R \
-i $main/seurat_processed_celltypes.rds \
-v $main/input_data/vdj \
-p $project_name \
-o $vdj \
-r false \
-s seurat_vdj.rds