#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

#args = commandArgs(trailingOnly=TRUE)
#subproject_name <- args[1]
#sample_names <- args[2]
#subset_data <- as.logical(args[3])

subproject_name <- "identify_normals"
sample_names <- c( "CID4463", "CID4513", "CID4515", "CID45171", "CID4523")
subset_data <- FALSE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample names = ", sample_names))
print(paste0("Subset data? ", as.character(subset_data)))

Rstudio = F

if (Rstudio) {
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(fpc)
} else {
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  library(cluster, lib.loc = lib_loc)
  library(ComplexHeatmap, lib.loc=lib_loc)
  library(circlize, lib.loc = lib_loc)
  library(scales, lib.loc = lib_loc)
  library(fpc, lib.loc = lib_loc)
}

library(Seurat)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)

project_name <- "identify_epithelial"
if (Rstudio) {
  home_dir <- "/Users/jamestorpy/clusterHome/"
} else {
  home_dir <- "/share/ScratchGeneral/jamtor/"
}
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/", subproject_name, "/functions/")
results_dir <- paste0(project_dir, "results/", subproject_name, "/")
out_dir <- paste0(results_dir, "group_infercnv/")
Robject_dir <- paste0(out_dir, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))

print(paste0("Sample directory = ", sample_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))

print("Running InferCNV identify normals pipeline on: ")
print(sample_names)


################################################################################
### 0. Define functions ###
################################################################################

create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))


################################################################################
### 1. Load heatmap and annotation data ###
################################################################################

# determine order of samples by subtype:
luminal <- c("CID3941", "CID3948", "CID4067", "CID4290A", "CID4461", "CID4463", 
  "CID4465", "CID4471", "CID4530", "CID4535")
HER2 <- c("CID3586", "CID3921", "CID3963", "CID4066", "CID4398", "CID45171")
TNBC <- c("CID44041", "CID4495", "CID44971", "CID44991", "CID4515")
metaplastic <- c("CID4513", "CID4523")
sample_names <- c(luminal, HER2, TNBC, metaplastic)
sample_names <- sample_names[
  sample_names %in% c( "CID4463", "CID4513", "CID4515", "CID45171", "CID4523")
]


for (s in 1:length(sample_names)) {
  print(s)

  heatmap_data <- readRDS(
  	paste0(results_dir, "seurat_", sample_names[[s]], 
  	  "/Output/InferCNV/Rdata/final_combined_heatmap_object.RData")
  )

  if (s==1) {
  	heatmap_dfs <- list(heatmap_data$heatmap_df)
  	names(heatmap_dfs) <- sample_names[s]
  	gene_list <- colnames(heatmap_data$heatmap_df)
  	group_annotations <- list(heatmap_data$group_annotation)
  	names(group_annotations) <- sample_names[s]
  	QC_annotations <- list(heatmap_data$QC_annotation)
  	names(QC_annotations) <- sample_names[s]
  } else {
  	heatmap_dfs[[s]] <- heatmap_data$heatmap_df
  	names(heatmap_dfs)[s] <- sample_names[s]
  	gene_list <- c(gene_list, colnames(heatmap_data$heatmap_df))
  	group_annotations[[s]] <- heatmap_data$group_annotation
  	names(group_annotations)[s] <- sample_names[s]
  	QC_annotations[[s]] <- heatmap_data$QC_annotation
  	names(QC_annotations)[s] <- sample_names[s]
  }

}

# get complete gene list from samples and order:
gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"))
gene_list <- unique(gene_list)
gene_list <- as.character(gene_order$V1[gene_order$V1 %in% gene_list])

# add all missing genes as columns to all heatmap dfs:
complete_heatmap_dfs <- lapply(heatmap_dfs, function(x) {
  
  missing_genes <- gene_list[!(gene_list %in% colnames(x))]
  missing_genes_df <- data.frame(matrix(NA, nrow = nrow(x), ncol = length(missing_genes)))
  colnames(missing_genes_df) <- missing_genes

  complete_df <- cbind(x, missing_genes_df)
  m <- match(gene_list, colnames(complete_df))
  complete_df <- complete_df[,m]

  return(complete_df)

})

# collate all dfs and check rows line up:
group_heatmap_df <- do.call(rbind, complete_heatmap_dfs)
rownames(group_heatmap_df) <- gsub("^.*\\.C", "C", rownames(group_heatmap_df))
print("Are all rows present in group heatmap df?")
identical(nrow(group_heatmap_df), sum(unlist(lapply(complete_heatmap_dfs, nrow))))

# collate all group annotations and add sample and subtype columns:
group_annotation_df <- do.call(rbind, group_annotations)
print("Are rows of group_annotation_df ordered the same as group_heatmap_df?")
identical(nrow(group_annotation_df), nrow(group_heatmap_df))
rownames(group_annotation_df) <- gsub("^.*\\.C", "C", rownames(group_annotation_df))
group_annotation_df$sample <- gsub("_.*$", "", rownames(group_annotation_df))
group_annotation_df$sample <- gsub("\\.log2", "", group_annotation_df$sample)
group_annotation_df$subtype <- NA
group_annotation_df$subtype[group_annotation_df$sample %in% luminal] <- "luminal"
group_annotation_df$subtype[group_annotation_df$sample %in% HER2] <- "HER2"
group_annotation_df$subtype[group_annotation_df$sample %in% TNBC] <- "TNBC"
group_annotation_df$subtype[group_annotation_df$sample %in% metaplastic] <- "metaplastic"

# reorder group_annotations_df by subtype:
group_annotation_df <- rbind(
  group_annotation_df[group_annotation_df$subtype=="luminal",],
  group_annotation_df[group_annotation_df$subtype=="HER2",],
  group_annotation_df[group_annotation_df$subtype=="TNBC",],
  group_annotation_df[group_annotation_df$subtype=="metaplastic",]
)

# collate all QC annotations:
group_QC_df <- do.call(rbind, QC_annotations)
print("Are rows of group_QC_df ordered the same as group_heatmap_df?")
identical(nrow(group_QC_df), nrow(group_heatmap_df))
rownames(group_QC_df) <- gsub("^.*\\.C", "C", rownames(group_QC_df))

# collate all annotations:
all_annotation_df <- merge(group_annotation_df, group_QC_df, by="row.names")
m <- match(rownames(group_annotation_df), all_annotation_df$Row.names)
all_annotation_df <- all_annotation_df[m,]
all_annotation_df$Row.names <- as.character(all_annotation_df$Row.names)

# create group annotation object:
extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  col_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
    "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
    "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
    "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
    "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
    "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
subtype_annotation_df <- data.frame(row.names=all_annotation_df$Row.names, 
  subtype=all_annotation_df$subtype)
cluster_number <- length(unique(subtype_annotation_df$subtype))
cluster_cols <- col_palette[1:cluster_number]
subtype_annotation <- Heatmap(subtype_annotation_df, col = cluster_cols, 
  name = "subtype_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(title = "Subtype", 
  title_gp = gpar(fontsize = 8, fontface = "bold"), 
  labels_gp = gpar(fontsize = 8)))

sample_annotation_df <- data.frame(row.names=all_annotation_df$Row.names, 
  sample=all_annotation_df$sample)
cluster_number <- length(unique(sample_annotation_df$sample))
cluster_cols <- col_palette[1:cluster_number]
sample_annotation <- Heatmap(sample_annotation_df, col = cluster_cols, 
  name = "sample_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(title = "Sample", 
  title_gp = gpar(fontsize = 8, fontface = "bold"), 
  labels_gp = gpar(fontsize = 8)))

group_annotation_df <- data.frame(row.names=all_annotation_df$Row.names, 
  group=all_annotation_df$group)
cluster_number <- length(unique(group_annotation_df$group))
cluster_cols <- col_palette[1:cluster_number]
group_annotation <- Heatmap(group_annotation_df, col = cluster_cols, 
  name = "group_annotation", 
  width = unit(4, "mm"), 
  show_row_names = F, show_column_names = F,
  heatmap_legend_param = list(title = "Cell type", 
  title_gp = gpar(fontsize = 8, fontface = "bold"), 
  labels_gp = gpar(fontsize = 8)))

# order group_heatmap_df:
group_heatmap_df <- group_heatmap_df[all_annotation_df$, ]

# create QC annotation objects:
group_QC_annotation <- create_QC_annotation(group_QC_df, group_heatmap_df)

# prepare df for plotting:
plot_object <- group_heatmap_df
colnames(plot_object) <- rep("la", ncol(plot_object))
# define heatmap colours:
heatmap_cols <- colorRamp2(c(-1, 0, 1), c("#00106B", "white", "#680700"), space = "sRGB")

print("Generating final heatmap...")
# create main CNV heatmap:
final_heatmap <- Heatmap(
  plot_object, name = paste0("hm"), 
  col = heatmap_cols,
  cluster_columns = F, cluster_rows = F,
  split = all_annotation_df$group,
  show_row_names = F, show_column_names = T,
  column_names_gp = gpar(col = "white"),
  show_row_dend = FALSE,
  heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
  grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
  title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
  use_raster = T, raster_device = c("png")
)
# determine co-ordinates of vertical lines at chromosome borders:
chr_data <- fetch_chromosome_boundaries(group_heatmap_df, ref_dir)
# determine co-ordinates of horizontal lines at group borders:
all_annotation_df$group_and_sample <- paste0(all_annotation_df$group, "_", 
	all_annotation_df$sample)
spl_groups <- split(all_annotation_df$group_and_sample, 
all_annotation_df$group_and_sample)
spl_groups <- spl_groups[unique(all_annotation_df$group_and_sample)]
if (length(spl_groups) > 1) {
  for ( n in 1:(length(spl_groups)-1) ) {
    if (n==1) {
      hlines <- c(length(spl_groups[[n]])/length(all_annotation_df$group_and_sample))
    } else {
      hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(all_annotation_df$group_and_sample)
    }
  }
  hlines <- 1-hlines
} else {
  hlines <- c(length(spl_groups[[1]])/length(all_annotation_df$group_and_sample))
}
  
  
################################################################################
### 5. Plot heatmap and annotations ###
################################################################################

ht_list <- subtype_annotation + sample_annotation + group_annotation + final_heatmap + 
  group_QC_annotation$nUMI_annotation + group_QC_annotation$nGene_annotation

annotated_heatmap <- grid.grabExpr(
  draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
)

# determine where starting co-ordinates for heatmap are based upon longest cluster name
# (0.00604 units per character):
longest_cluster_name <- max(nchar(unique(as.character(all_annotation_df$group))))
x_coord <- longest_cluster_name*0.0037
# save final heatmap objects with measurements needed to create it:
annotated_heatmap_and_measurements <- list(
	annotated_heatmap = annotated_heatmap,
	group_heatmap_df = group_heatmap_df,
  group_annotation = group_annotation_df,
  QC_annotation = group_QC_df,
	chr_data = chr_data, 
	hlines = hlines, 
	x_coord = x_coord
)
saveRDS(annotated_heatmap_and_measurements, paste0(Robject_dir, "final_combined_heatmap_object.RData"))

# plot final annotated heatmap:
pdf(paste0(plot_dir, "final_combined_infercnv_heatmap_subset.pdf"), height = 14.5, width = 18)   
grid.newpage()
  pushViewport(viewport(x = 0, y = 0.15, 
                       width = 0.99, height = 0.8, just = c("left", "bottom")))
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
   pushViewport(viewport(x=x_coord + 0.917, y=0.1, width = 0.1, height = 0.1, just = "bottom"))
     grid.text("nUMI", rot=65)
   popViewport()
   pushViewport(viewport(x=x_coord + 0.935, y=0.1, width = 0.1, height = 0.1, just = "bottom"))
     grid.text("nGene", rot=65)
   popViewport()
    
dev.off()
