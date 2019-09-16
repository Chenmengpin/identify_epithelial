numcores=3
sample_name <- "CID4515"
subset_data <- FALSE
HMM <- TRUE
subpop_mode <- TRUE
combined_sample <- FALSE

library(Seurat)
library(DropletUtils)
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)
library(HiddenMarkov, lib.loc=lib_loc)
library(purrr)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)
library(scales)

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_030719"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/pipeline_v3/functions/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")
sample_dir <- paste0(seurat_path, "/seurat_", sample_name, "/Output/")
setwd(sample_dir)
new_seurat_dir <- paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/supervised_clustering/")
  system(paste0("mkdir -p ", new_seurat_dir))

out_dir <- paste0("InferCNV/supervised_clustering/")
if (subpop_mode) {
  out_dir <- paste0(out_dir, "/subpop_mode/")
}

input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))

print(paste0("Running InferCNV pipeline on ", sample_name, "..."))


################################################################################
### 0. Define functions ###
################################################################################

temp_png_function <- dget(paste0(func_dir, "temp_png_function.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
create_infercnv_object <- dget(paste0(func_dir, "create_infercnv_object.R"))
run_infercnv <- dget(paste0(func_dir, "run_infercnv.R"))
run_infercnv <- possibly(run_infercnv, otherwise = "error_occurred")
create_group_annotation <- dget(paste0(func_dir, "create_group_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
annotate_PAM50_CNV <- dget(paste0(func_dir, "annotate_PAM50_CNV.R"))
create_CNV_genes_annotation <- dget(paste0(func_dir, "create_CNV_genes_annotation.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))


################################################################################
### 1. Load seurat objects ###
################################################################################

# load seurat object:
if (combined_sample) {
  seurat_10X <- readRDS(paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/04_seurat_object_combined.RData"))
} else {
  seurat_10X <- readRDS(paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/03_seurat_object_processed.RData"))
  cell_ids <- names(Idents(seurat_10X))
  Idents(seurat_10X) <- paste0(sample_name, "_", 
    seurat_10X@meta.data$PC_A_res.1)
}


################################################################################
### 2. Generate input matrix and metadata files ###
################################################################################

if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
  print("Creating inferCNV raw counts file...")
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  if (subset_data) {
    count_df <- count_df[1:500, 1:500]
  }
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
  sep="\t", col.names=T, row.names=T)
} else {
  count_df <- read.table(paste0(input_dir, "input_matrix.txt"))
}

print("Creating inferCNV metadata file...")
infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, combined_sample,
  "PC_A_res.1", subset_data=subset_data, count_df)
names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
seurat_10X <- infercnv_metadata$seurat

normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
    unique(Idents(seurat_10X)[names(Idents(seurat_10X)) %in% colnames(count_df)]), value=T, 
    invert=T)

print(paste0("Normals are: ", normals))

cells_to_remove <- 
infercnv_metadata$metadata$cell_ids[grep("CAF", infercnv_metadata$metadata$cell_type)]
infercnv_metadata$metadata <- infercnv_metadata$metadata[!(infercnv_metadata$metadata$cell_ids %in% 
  cells_to_remove),]
print(dim(infercnv_metadata$metadata))

if (HMM) {
  if (subpop_mode) {
    if (subset_data) {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/subpop_mode/subset/"), 
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt", 
      full.names = T)
    } else {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/subpop_mode/"), 
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt", 
      full.names = T)
    }
  } else {
    if (subset_data) {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/subset/"), 
      pattern = "infercnv.12_denoised.observations.txt", full.names = T)
    } else {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/"), 
      pattern = "infercnv.12_denoised.observations.txt", full.names = T)
    }
  }
} else {
  if (subpop_mode) {
    if (subset_data) {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/subset/"), 
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt", 
      full.names = T)
    } else {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/"), 
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt", 
      full.names = T)
    }
  } else {
    if (subset_data) {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/subset/"), 
      pattern = "infercnv.12_denoised.observations.txt", full.names = T)
    } else {
      infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      sample_name, "/Output/InferCNV/supervised_clustering/"), 
      pattern = "infercnv.12_denoised.observations.txt", full.names = T)
    }
  }
}

infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))
heatmap_df <- infercnv_output[!(rownames(infercnv_output) %in% cells_to_remove),]
print(dim(heatmap_df))
# remove 'luminal' from epithelial cells as Garnett labels everything luminal:
infercnv_metadata$metadata$cell_type <- gsub("Luminal_", "", infercnv_metadata$metadata$cell_type)
# create group annotation df
group_annotation <- create_group_annotation(heatmap_df, infercnv_metadata$metadata)
# ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
heatmap_df <- heatmap_df[m,]
chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)

create_GIN_annotation <- function(df, normal_cells) {

  # adjust HMM values
  # 0 = complete loss, change to 2
  # 0.5 = loss of one copy, change to 1
  # 1 = neutral, change to 0
  # 1.5 = addition of one copy, change to 1
  # 2 = addition of two copies, change to 2
  # 3 = addition of three or more copies, change to 3
  old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
  new_scores <- c(2, 1, 0, 1, 2, 3)
  for (s in 1:length(new_scores)) {
    df[df == old_scores[s]] <- new_scores[s]
  }

  # sum total CNV levels for each cell:
  GIN_levels <- apply(df, 1, function(x) {
    x[is.na(x)] <- 0
    return(sum(x))
  })

  # normalise to number of genes in plot and scale to 1:100:
  GIN_levels <- GIN_levels/ncol(df)
  max_GIN_level <- 3*ncol(df)/ncol(df)

  GIN_levels <- round((GIN_levels/max_GIN_level)*100, 1)

  
  # create GIN annotation:
  GIN_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      GIN_levels,
      gp = gpar(
        col = "#AF548E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )

  # record maximum levels of GIN in normals:
  max_normal_GIN <- max(GIN_levels[normal_cells])

  # create GIN levels data frame:
  GIN_levels <- data.frame(cell_id=rownames(df), GIN=GIN_levels)

  result_list <- list(GIN_annotation, GIN_levels, max_normal_GIN)
  names(result_list) <- c("GIN_annotation", "GIN_levels", "max_normal_GIN")

  return(result_list)
}


if (HMM) {
  # determine normal cell ids:
  normal_cells <- infercnv_metadata$metadata$cell_ids[
    infercnv_metadata$metadata$cell_type %in% normals
  ]
  # add GIN annotation:
  GIN_annotation <- create_GIN_annotation(heatmap_df, normal_cells)
}

create_correlation_annotation <- function(df, metadata, GIN_levels) {
  # determine top 5% of tumour cells:
  epithelial_cells <- metadata$cell_ids[grep("pithelial", metadata$cell_type)]
  epithelial_GIN <- GIN_levels[GIN_levels$cell_id %in% epithelial_cells,]
  epithelial_GIN <- epithelial_GIN[order(epithelial_GIN$GIN, decreasing=T),]
  top_GIN <- head(epithelial_GIN, nrow(epithelial_GIN)*0.05)
  # find average genome-wide CNV predictions across genome:
  top_GIN_CNV_average <- apply(df[top_GIN$cell_id,], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  GIN_correlations <- apply(df, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_GIN_CNV_average)
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  cor_df <- do.call("rbind", GIN_correlations)
  cor_df <- cbind(cor_df, GIN_levels$GIN)
  return(cor_df)
}

all_correlations <- create_correlation_annotation(heatmap_df, 
  infercnv_metadata$metadata, GIN_annotation$GIN_levels)
colnames(all_correlations) <- c("cor.estimate", "cor.p.value",
  "GIN_levels")

normal_epithelial_cluster <- "Epithelial_17"
normal_epithelial_ids <- infercnv_metadata$metadata$cell_ids[
  infercnv_metadata$metadata$cell_type == normal_epithelial_cluster
]
normal_epi_correlation <- all_correlations[normal_epithelial_ids,]
normal_epi_correlation <- normal_epi_correlation[
  normal_epi_correlation$cor.p.value < 0.05, 
]
normal_correlation_range <- round(range(normal_epi_correlation$cor.estimate), 2)
normal_GIN_range <- range(normal_epi_correlation$GIN_levels)

other_epi_correlation <- all_correlations[
  !(rownames(all_correlations) %in% normal_epithelial_ids),
]
other_epi_correlation <- other_epi_correlation[
  other_epi_correlation$cor.p.value < 0.05, 
]
other_correlation_range <- round(range(other_epi_correlation$cor.estimate), 2)
other_GIN_range <- range(other_epi_correlation$GIN_levels)

result_df <- data.frame(rbind(normal_GIN_range, normal_correlation_range,
  other_GIN_range, other_correlation_range))

write.table(result_df, paste0(out_dir, "/GIN_and_correlation_ranges.txt"), 
  row.names=T, col.names=F, quote=F)


