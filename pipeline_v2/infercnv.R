#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#args = commandArgs(trailingOnly=TRUE)
#print(args)
#numcores=as.numeric(args[1])
#sample_name <- args[2]
#subset_data <- as.logical(args[3])
#HMM <- as.logical(args[4])

numcores=12
sample_name <- "integrated"
subset_data <- FALSE
HMM <- FALSE

if ( length(grep("PDX", sample_name)) != 0 ) {
  sample_type <- "PDX"
} else if ( length(grep("IND", sample_name)) != 0 ) {
  sample_type <- "IND"
} else {
  sample_type <- "primary"
}

date()

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

project_name <- "identify_epithelial"
subproject_name <- "v2/CAF_paper/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/pipeline_v2/functions/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")
sample_dir <- paste0(seurat_path, "/seurat_", sample_name, "/Output/")
setwd(sample_dir)
new_seurat_dir <- paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/HMM/")
  system(paste0("mkdir -p ", new_seurat_dir))
if (HMM) {
  out_dir <- paste0("InferCNV/HMM/")
} else {
  out_dir <- paste0("InferCNV/")
}

input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))

print(paste0("Running InferCNV pipeline on ", sample_name, ", filtering out ", 
  "non-epithelial cells..."))


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
create_GIN_annotation <- dget(paste0(func_dir, "create_GIN_annotation.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))


################################################################################
### 1. Load seurat objects ###
################################################################################

# load seurat object:
#seurat_10X <- readRDS(paste0(seurat_path, sample_name, 
#"/Output/Rdata/03_seurat_object_processed.RData"))
#cell_ids <- names(seurat_10X@ident)
#cell_types <- gsub("_.*[0-9]$", "", seurat_10X@meta.data$celltype)
#seurat_10X@ident <- factor(paste0(cell_types, 
#  "_", seurat_10X@ident))
seurat_10X <- readRDS("/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/v2/CAF_paper/original/Rdata/seurat_CCA_aligned_processed_filtered_annotated_doublet_filtered_newstromalIDs_v2.Rdata")
cell_ids <- names(seurat_10X@ident)
seurat_10X@ident <- seurat_10X@meta.data$celltype_3
names(seurat_10X@ident) <- cell_ids


################################################################################
### 2. Generate input matrix and metadata files ###
################################################################################

if (subset_data) {
  if (HMM) {
    out_dir <- paste0("InferCNV/HMM/subset")
    } else {
      out_dir <- paste0("InferCNV/subset")
    }
  
  input_dir <- paste0(out_dir, "/input_files/")
  system(paste0("mkdir -p ", input_dir))
}

# create raw matrix input file
if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
  print("Creating inferCNV raw counts file...")
  count_df <- as.matrix(seurat_10X@raw.data)
  if (subset_data) {
    count_df <- count_df[1:500, 1:500]
  }
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
  sep="\t", col.names=T, row.names=T)
} else {
  count_df <- read.table(paste0(input_dir, "input_matrix.txt"))
}

# prepare infercnv metadata and write to files
#print("Creating inferCNV metadata file...")
#
#infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, "celltype_3", 
#  subset_data=F, count_df)
#
#names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
#seurat_10X <- infercnv_metadata$seurat

infercnv_metadata <- read.table("/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/v2/CAF_paper/integrated/Output/InferCNV/input_files/metadata.txt")
colnames(infercnv_metadata) <- c("cell_ids", "cell_type")

# collapse stromals into CAFs, immune and endothelial groups:
metadata_cell_type <- as.character(infercnv_metadata$cell_type)
metadata_cell_type[grep("CAF", metadata_cell_type)] <- paste0(
  gsub(
    "_.*$", "", metadata_cell_type[grep("CAF", metadata_cell_type)]
  ), "_CAF"
)

metadata_cell_type[
  grep("Plasma|Myeloid|[t,T]_[c,C]ell|[b,B]_[c,C]ell|[t,T]_Reg", metadata_cell_type)
] <- paste0(
  gsub(
    "_.*$", "", metadata_cell_type[
      grep("Plasma|Myeloid|[t,T]_[c,C]ell|[b,B]_[c,C]ell|[t,T]_Reg", metadata_cell_type)
    ]
  ), "_Immune"
)

# generate cluster metric plots for epithelial clusters:
#epithelial_clusters <- grep("pithelial", levels(seurat_10X@ident), value=T)
#temp_violinplot <- VlnPlot(
#    object = seurat_10X,
#    features.plot = c("nGene",
#                      "nUMI",
#                      "percent.mito"),
#    nCol = 3, 
#    group.by = "PCs_1_0.4",
#    ident.include = epithelial_clusters
#  )
#png("InferCNV/metrics_by_epithelial_cluster.png", width = 14, 
#      height = 8, res = 300, units = 'in')
#  print(temp_violinplot)
#dev.off()

# if no epithelial or CAF clusters present, run without normals:
#if ( length(grep("pithelial",  
#  infercnv_metadata$number_per_group$Var1)) < 1 ) {
#  print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
#
#}
#print(paste0("Idents are: ", levels(seurat_10X@ident)))
#write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
#  quote=F, sep="\t", col.names=F, row.names=F)
#write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
#  "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
#saveRDS(seurat_10X, paste0(new_seurat_dir, 
#  "05_seurat_object_annotated.RData"))
## define normals as stromal groups, excluding epithelial and unassigned 
## cells:
#
#normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
#    unique(seurat_10X@ident[names(seurat_10X@ident) %in% 
#    colnames(count_df)]), value=T, invert=T)

normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
    unique(infercnv_metadata$cell_type[infercnv_metadata$cell_ids %in% 
    colnames(count_df)]), value=T, invert=T)


################################################################################
### 3. Run InferCNV ###
################################################################################
  
if ( length(list.files(paste0(seurat_path, 
    sample_name, "/Output/InferCNV/"), 
    pattern = "infercnv.12_denoised.observations.txt", full.names = T)) == 0 ) {

  print("Creating inferCNV object...")
  raw_path <- paste0(input_dir, "input_matrix.txt")
  annotation_path <- paste0(input_dir, "metadata.txt")
  gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- create_infercnv_object(raw_path, annotation_path,
    gene_path, normals)
  
  print("InferCNV object created, running inferCNV...")
  infercnv_output <- run_infercnv(initial_infercnv_object, numcores-1, out_dir, 0.1, 
    101, 3, 1.3, HMM, "i6")

}
  

################################################################################
### 4. Create heatmap and annotations ###
################################################################################

# fetch vector of all genes:
all_genes <- as.character(read.table(paste0(ref_dir, "/infercnv_gene_order.txt"))$V1)

# remove CAFs from heatmap df and metadata:
#cells_to_remove <- 
#infercnv_metadata$cell_ids[grep("pithelial", infercnv_metadata$metadata$cell_type, 
#  invert=T)]
#infercnv_metadata <- infercnv_metadata[!(infercnv_metadata$cell_ids %in% 
#  cells_to_remove),]
#print(dim(infercnv_metadata))

if (HMM) {
  if (subset_data) {
  infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
    sample_name, "/Output/InferCNV/HMM/subset/"), 
    pattern = "infercnv.14_HMM_predHMMi6.hmm_mode-samples.repr_intensities.observations.txt", full.names = T)
  } else {
    infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
    sample_name, "/Output/InferCNV//HMM/"), 
    pattern = "infercnv.14_HMM_predHMMi6.hmm_mode-samples.repr_intensities.observations.txt", full.names = T)
  }
} else {
  if (subset_data) {
    infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
    sample_name, "/Output/InferCNV//HMM/subset/"), 
    pattern = "denoised.observations.txt", full.names = T)
    infercnv_stromal_filename <- list.files(paste0(seurat_path, "seurat_", 
    sample_name, "/Output/InferCNV//HMM/subset/"), 
    pattern = "infercnv.12_denoised.references", full.names = T)
  } else {
    infercnv_output_filename <- list.files(paste0(seurat_path, 
    sample_name, "/Output/InferCNV/"), 
    pattern = "denoised.observations.txt", full.names = T)
    infercnv_stromal_filename <- list.files(paste0(seurat_path, 
    sample_name, "/Output/InferCNV/"), 
    pattern = "infercnv.12_denoised.references", full.names = T)
  }
}

infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))
infercnv_refs <- as.data.frame(t(read.table(infercnv_stromal_filename)))
heatmap_df <- rbind(infercnv_output, infercnv_refs)
#heatmap_df <- infercnv_output[!(rownames(infercnv_output) %in% cells_to_remove),]
print(dim(heatmap_df))
# remove 'luminal' from epithelial cells as Garnett labels everything luminal:
#infercnv_metadata$metadata$cell_type <- gsub("Luminal_", "", infercnv_metadata$metadata$cell_type)
# create group annotation df
group_annotation <- create_group_annotation(heatmap_df, infercnv_metadata)
# ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
heatmap_df <- heatmap_df[m,]
chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)

## add GIN annotation:
#GIN_annotation <- create_GIN_annotation(heatmap_df)
  
QC_annotation <- create_QC_annotation(seurat_10X, heatmap_df)

## create heatmap annotation for genome-wide PAM50 subtype CNV frequency
#PAM50_subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")
#METABRIC_CNV_frequencies <- read.table(paste0(ref_dir, "infercnv_metabric_cnv.txt"), header=T, as.is=T, fill=T)
#for ( i in 1:length(PAM50_subtypes) ) {
#    print(paste0("Generating ", PAM50_subtypes[i], " CNV plot..."))
#    if (i==1) {
#      metabric_plots <- 
#        list(annotate_PAM50_CNV(heatmap_df, METABRIC_CNV_frequencies, 
#        PAM50_subtypes[i], chr_data$ends, chr_data$lengths))
#    } else {
#      metabric_plots[[i]] <- 
#        annotate_PAM50_CNV(heatmap_df, METABRIC_CNV_frequencies, 
#        PAM50_subtypes[i], chr_data$ends, chr_data$lengths)
#    }
#}
#names(metabric_plots) <- PAM50_subtypes
## create heatmap annotation for gain and loss-associated genes, 
##collated by Niantao
## read in CNV_genes
#CNV_genes <- read.table(paste0(ref_dir, 
#  "./infercnv_brca_genes_associated_with_CNVs.txt"), header = T, as.is = T)
## create CNV_genes annotation:
#print("Annotating CNV-associated genes...")
#CNV_genes_annotation <- create_CNV_genes_annotation(heatmap_df, CNV_genes)

plot_object <- heatmap_df
#old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
#  new_scores <- c(-2, -1, 0, 1, 2, 3)
#  for (s in 1:length(new_scores)) {
#    df[df == old_scores[s]] <- new_scores[s]
#  }
colnames(plot_object) <- rep("la", ncol(plot_object))

print("Generating final heatmap...")

na_less_vector <- unlist(plot_object)
na_less_vector <- na_less_vector[!is.na(na_less_vector)]
# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = colorRamp2(c(min(na_less_vector), 1, max(na_less_vector)), 
    c("#00106B", "white", "#680700"), space = "sRGB"),
  cluster_columns = F, cluster_rows = F,
  split = group_annotation$group_annotation_df$group,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  #bottom_annotation = CNV_genes_annotation, bottom_annotation_height = unit(2, "cm"),
  #gap = unit(1, "cm"),
  heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
  grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
  use_raster = T, raster_device = c("png")
)
# determine co-ordinates of horizontal lines at group borders:
spl_groups <- split(group_annotation$group_annotation_df$group, 
group_annotation$group_annotation_df$group)
spl_groups <- spl_groups[unique(group_annotation$group_annotation_df$group)]
if (length(spl_groups) > 1) {
  for ( n in 1:(length(spl_groups)-1) ) {
    if (n==1) {
      hlines <- c(length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group))
    } else {
      hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group)
    }
  }
  hlines <- 1-hlines
} else {
  hlines <- c(length(spl_groups[[1]])/length(group_annotation$group_annotation_df$group))
}


################################################################################
### 5. Generate heatmap and annotations ###
################################################################################

#if ( length(unique(GIN_annotation$GIN_levels)) > 1 ) {
#  ht_list <- group_annotation$group_annotation + final_heatmap + 
#  GIN_annotation$GIN_annotation + QC_annotation$nUMI_annotation + 
#  QC_annotation$nGene_annotation
#} else {
  ht_list <- group_annotation$group_annotation + final_heatmap + 
  QC_annotation$nUMI_annotation + 
  QC_annotation$nGene_annotation
#}

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
x_coord <- longest_cluster_name*0.0037

pdf(paste0(out_dir, "final_infercnv_heatmap.pdf"), 
  height = 14.5, width = 18)
  grid.newpage()
  # plot Normal subtype:
  pushViewport(viewport(x = x_coord+0.003, y = 0.090,
                        width = 0.835+0.019, height = 0.08, just = c("left", "top")))
  grid.draw(metabric_plots[[5]])
  popViewport()
  
  # plot Basal subtype:
  pushViewport(viewport(x = x_coord+0.008, y = 0.169,
                        width = 0.835+0.014, height = 0.08, just = c("left", "top")))
  grid.draw(metabric_plots[[4]])
  popViewport()
  
  # plot Her2 subtype:
  pushViewport(viewport(x = x_coord+0.01, y = 0.251, 
                        width = 0.835+0.012, height = 0.08, just = c("left", "top")))
  grid.draw(metabric_plots[[3]])
  popViewport()
  
  # plot LumB subtype:
  pushViewport(viewport(x = x_coord+0.007, y = 0.333, 
                        width = 0.835+0.015, height = 0.08, just = c("left", "top")))
  grid.draw(metabric_plots[[2]])
  popViewport()
  
  # plot LumA subtype:
  pushViewport(viewport(x = x_coord+0.007, y = 0.411, 
                        width = 0.835+0.015, height = 0.08, just = c("left", "top")))
  grid.draw(metabric_plots[[1]])
  popViewport()
  pushViewport(viewport(x = 0, y = 0.4, 
                        width = 0.99, height = 0.6, just = c("left", "bottom")))
  grid.draw(annotated_heatmap)
  decorate_heatmap_body("hm", {

    for ( e in 1:length(chr_data$end_pos) ) {
      grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 1, 
        col = "#383838"))
      grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
        unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=8))
    }
    for ( m in 1:length(hlines) ) {
      grid.lines(c(0, 1), c(hlines[m], hlines[m]), gp = gpar(lwd = 1, col = "#383838"))
    }
  })

  popViewport()
  pushViewport(viewport(x=x_coord + 0.865, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
  #grid.draw(lollipop)
  grid.text("GIN", rot=55)
  popViewport()
  pushViewport(viewport(x=x_coord + 0.88, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
  #grid.draw(lollipop)
  grid.text("nUMI", rot=55)
  popViewport()
  pushViewport(viewport(x=x_coord + 0.895, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
  #grid.draw(lollipop)
  grid.text("nGene", rot=55)
  popViewport()
  
dev.off()

# convert pdf to png:
system(paste0("for p in ", out_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
  "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
  "convert -density 150 ", out_dir, "$f -quality 90 ", out_dir, "$new; done"))

print(paste0("Group heatmap created, output in ", out_dir))


#pdf(paste0(out_dir, "filtered_epithelial_infercnv_heatmap.pdf"), 
#  height = 14.5, width = 17)
#   # plot heatmap:
#  pushViewport(viewport(x = 0.5, y = 0.1, 
#                        width = 1, height = 0.8, just = "bottom"))
#  grid.draw(annotated_heatmap)
#  decorate_heatmap_body("hm", {
#
#    for ( e in 1:length(chr_data$end_pos) ) {
#      grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 1, 
#        col = "#383838"))
#      grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
#        unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=8))
#    }
#    for ( m in 1:length(hlines) ) {
#      grid.lines(c(0, 1), c(hlines[m], hlines[m]), gp = gpar(lwd = 1, col = "#383838"))
#    }
#  })
#
#  popViewport()
#
#dev.off()


  
date()
  
  