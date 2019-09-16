#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
numcores=as.numeric(args[1])
print(paste0("number cores = ", numcores))
sample_name <- args[2]
print(sample_name)
subproject_name <- args[3]
print(subproject_name)
subset_data <- as.logical(args[4])
print(subset_data)
temp_resolution <- args[5]
print(subproject_name)
supervised_clustering <- args[6]
print(supervised_clustering)
remove_cells_from_plot <- args[7]

numcores=10
sample_name <- "CID4463"
subproject_name <- "brca_mini_atlas_060819"
subset_data <- FALSE
temp_resolution <- "PC_A_res.1"
supervised_clustering <- TRUE
remove_cells_from_plot <- "CAF|nknown|nassigned"

library(Seurat)
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)
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
  "/Output/Rdata/non_HMM/supervised_clustering/")
  system(paste0("mkdir -p ", new_seurat_dir))

print(paste0("Running InferCNV pipeline on ", sample_name, ", filtering out ", 
  "non-epithelial cells..."))


################################################################################
### 0. Define functions ###
################################################################################

temp_png_function <- dget(paste0(func_dir, "temp_png_function.R"))
prepare_infercnv_metadata_non_HMM <- dget(paste0(func_dir, 
  "prepare_infercnv_metadata_non_HMM.R"))
create_infercnv_object <- dget(paste0(func_dir, "create_infercnv_object.R"))
run__non_HMM_infercnv <- dget(paste0(func_dir, "run_non_HMM_infercnv.R"))
create_group_annotation <- dget(paste0(func_dir, "create_group_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, 
  "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
annotate_PAM50_CNV <- dget(paste0(func_dir, "annotate_PAM50_CNV.R"))
create_CNV_genes_annotation <- dget(paste0(func_dir, 
  "create_CNV_genes_annotation.R"))
create_infercnv_level_annotation <- dget(paste0(func_dir, 
  "create_infercnv_level_annotation.R"))
create_normal_call_non_HMM_annotation <- dget(paste0(func_dir, 
  "create_normal_call_non_HMM_annotation.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))


################################################################################
### 1. Load seurat objects ###
################################################################################

# load seurat object:
seurat_10X <- readRDS(paste0(seurat_path, "seurat_", sample_name, 
"/Output/Rdata/03_seurat_object_processed.RData"))
Idents(seurat_10X) <- paste0(sample_name, "_", 
  eval(parse(text = paste0("seurat_10X@meta.data$", temp_resolution))))


################################################################################
### 2. Generate input matrix and metadata files ###
################################################################################

if (supervised_clustering) {
  if (subset_data) {
    out_dir <- paste0("InferCNV/non_HMM/supervised_clustering/subset/")
  } else {
  	out_dir <- paste0("InferCNV/non_HMM/supervised_clustering/")
  }
} else {
  if (subset_data) {
    out_dir <- paste0("InferCNV/non_HMM/unsupervised_clustering/subset/")
  } else {
  	out_dir <- paste0("InferCNV/non_HMM/unsupervised_clustering/")
  }
}
input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))


# create raw matrix input file
count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
if (subset_data) {
  count_df <- count_df[1:500, 1:500]
}

# prepare infercnv metadata and write to files
print("Creating inferCNV metadata file...")
infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, combined_sample,
  temp_resolution, subset_data=subset_data, count_df, supervised_clustering)
names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
seurat_10X <- infercnv_metadata$seurat

# only keep cells in metadata df:
count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]

# create raw matrix input file
if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
  print("Creating inferCNV raw counts file...")
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
  sep="\t", col.names=T, row.names=T)
}

# generate cluster metric plots for epithelial clusters:
epithelial_clusters <- grep("pithelial", levels(Idents(seurat_10X)), value=T)
if (!file.exists("InferCNV/metrics_by_epithelial_cluster.png")) {
  temp_png_function("InferCNV/metrics_by_epithelial_cluster.png")
    temp_violinplot <- VlnPlot(
      object = seurat_10X,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
      pt.size = 1.5,
      idents = epithelial_clusters
    )
    print(temp_violinplot)
  dev.off()
}

print(paste0("Idents are: ", levels(Idents(seurat_10X))))
write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
  quote=F, sep="\t", col.names=F, row.names=F)
write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
  "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
saveRDS(seurat_10X, paste0(new_seurat_dir, 
  "04_seurat_object_annotated.RData"))

normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
  unique(Idents(seurat_10X)[names(Idents(seurat_10X)) %in% colnames(count_df)]), value=T, 
  invert=T)

print(paste0("Normals are: ", normals))


################################################################################
### 3. Run InferCNV ###
################################################################################
  
if ( length(list.files(paste0(seurat_path, "seurat_", 
  sample_name, "/Output/InferCNV/supervised_clustering/sample/"), 
  pattern = "infercnv.12_denoised.observations.txt", full.names = T)) == 0 ) {
  print("Creating inferCNV object...")
  raw_path <- paste0(input_dir, "input_matrix.txt")
  annotation_path <- paste0(input_dir, "metadata.txt")
  gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- create_infercnv_object(raw_path, annotation_path,
    gene_path, normals)
  
  print("InferCNV object created, running inferCNV...")
  infercnv_output <- run_infercnv(initial_infercnv_object, numcores-1, out_dir, 0.1, 
    101, 3, 1.3)
}


################################################################################
### 4. Create heatmap and annotations ###
################################################################################

infercnv_output_filename <- list.files(out_dir, 
   pattern = "infercnv.12_denoised.observations.txt", full.names = T)

infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))

# remove CAFs from infercnv_output:
cells_to_remove <- 
infercnv_metadata$metadata$cell_ids[
  grep(remove_cells_from_plot, infercnv_metadata$metadata$cell_type)
]
#infercnv_metadata$metadata <- infercnv_metadata$metadata[
#  !(infercnv_metadata$metadata$cell_ids %in% cells_to_remove),
#]
print(paste0("Cell no. in data pre-unwanted cell removal: ", nrow(infercnv_output)))
heatmap_df <- infercnv_output[!(rownames(infercnv_output) %in% cells_to_remove),]
print(paste0("Cell no. in data post-unwanted cell removal: ", nrow(heatmap_df)))

# create group annotation df
group_annotation <- create_group_annotation(heatmap_df, infercnv_metadata$metadata)
# ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
heatmap_df <- heatmap_df[m,]
chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)

# add infercnv level annotation:
infercnv_level_annotation <- create_infercnv_level_annotation(heatmap_df)
# add normal annotation:
normal_call_annotation <- create_normal_call_non_HMM_annotation(heatmap_df, 
  infercnv_metadata$metadata, infercnv_level_annotation$infercnv_level,
  0.3, 0.2, 0.05)

# save infercnv levels:
if (!file.exists(paste0(out_dir, "infercnv_level.txt"))) {
  write.table(infercnv_level_annotation$infercnv_level, 
  	paste0(out_dir, "infercnv_level.txt"), 
    quote=F, row.names=F)
}
# save correlation scores:
if (!file.exists(paste0(out_dir, "top_5%_cancer_correlation_scores.txt"))) {
  write.table(normal_call_annotation$correlation_df, 
  	paste0(out_dir, "top_5%_cancer_correlation_scores.txt"), quote=F, row.names=F)
}

QC_annotation <- create_QC_annotation(seurat_10X, heatmap_df)

# create heatmap annotation for genome-wide PAM50 subtype CNV frequency
PAM50_subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")
METABRIC_CNV_frequencies <- read.table(paste0(ref_dir, "infercnv_metabric_cnv.txt"), header=T, as.is=T, fill=T)
for ( i in 1:length(PAM50_subtypes) ) {
    print(paste0("Generating ", PAM50_subtypes[i], " CNV plot..."))
    if (i==1) {
      metabric_plots <- 
        list(annotate_PAM50_CNV(heatmap_df, METABRIC_CNV_frequencies, 
        PAM50_subtypes[i], chr_data$ends, chr_data$lengths))
    } else {
      metabric_plots[[i]] <- 
        annotate_PAM50_CNV(heatmap_df, METABRIC_CNV_frequencies, 
        PAM50_subtypes[i], chr_data$ends, chr_data$lengths)
    }
}
names(metabric_plots) <- PAM50_subtypes
# create heatmap annotation for gain and loss-associated genes, 
#collated by Niantao
# read in CNV_genes
CNV_genes <- read.table(paste0(ref_dir, 
  "./infercnv_brca_genes_associated_with_CNVs.txt"), header = T, as.is = T)
# create CNV_genes annotation:
print("Annotating CNV-associated genes...")
CNV_genes_annotation <- create_CNV_genes_annotation(heatmap_df, CNV_genes)

plot_object <- heatmap_df



