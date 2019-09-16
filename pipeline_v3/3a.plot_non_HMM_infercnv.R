#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
numcores=as.numeric(args[1])
print(paste0("number cores = ", numcores))
sample_name <- args[2]
print(sample_name)
subset_data <- as.logical(args[3])
print(subset_data)
combined_sample <- as.logical(args[6])
print(combined_sample)
subproject_name <- args[7]
print(subproject_name)

#numcores=30
#sample_name <- "CID4463"
#subset_data <- FALSE
#combined_sample <- FALSE
#subproject_name <- "brca_mini_atlas_060819"


library(Seurat)
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)
library(purrr)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)

project_name <- "identify_epithelial"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/pipeline_v3/functions/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")
sample_dir <- paste0(seurat_path, "/seurat_", sample_name, "/Output/")
print(sample_dir)
setwd(sample_dir)
new_seurat_dir <- paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/supervised_clustering/non_HMM/")
  system(paste0("mkdir -p ", new_seurat_dir))

out_dir <- paste0("InferCNV/supervised_clustering/", run_mode, "/")

input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))

print(paste0("Running InferCNV pipeline on ", sample_name, ", filtering out ", 
  "non-epithelial cells..."))