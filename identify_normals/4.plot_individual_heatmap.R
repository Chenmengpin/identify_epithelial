#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript
subproject_name <- "identify_normals"
sample_name <- "CID4463"
subset_data <- FALSE
draw_bulk_annotations <- FALSE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Draw bulk annotations? ", as.character(draw_bulk_annotations)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(dplyr)

project_name <- "identify_epithelial"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")
func_dir <- paste0(project_dir, "/scripts/", subproject_name, "/functions/")
results_dir <- seurat_path <- paste0(project_dir, "results/", subproject_name, "/")
sample_dir <- paste0(results_dir, "/seurat_", sample_name, "/Output/")
setwd(sample_dir)

out_dir <- paste0("InferCNV/")
input_dir <- paste0(out_dir, "/input_files/")
system(paste0("mkdir -p ", input_dir))
Robject_dir <- paste0(out_dir, "/Rdata/")
system(paste0("mkdir -p ", Robject_dir))
plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "/tables/")
system(paste0("mkdir -p ", table_dir))
integrated_dir <- paste0("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/",
	"Jun2019/04_reclustering_analysis/run06_v1.2.1/output/Epithelial/02_Rdata/")

print(paste0("Sample directory = ", sample_dir))
print(paste0("Reference directory = ", ref_dir))
print(paste0("R function directory = ", func_dir))
print(paste0("Output directory = ", out_dir))
print(paste0("R object directory = ", Robject_dir))
print(paste0("Plot directory = ", plot_dir))
print(paste0("Plot directory = ", table_dir))

print(paste0("Running InferCNV identify normals pipeline on ", sample_name))


################################################################################
### 0. Define functions ###
################################################################################

temp_png_function <- dget(paste0(func_dir, "temp_png_function.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
create_group_annotation <- dget(paste0(func_dir, "create_group_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))
create_array_CNV_annotation <- dget(paste0(func_dir, "create_array_CNV_annotation.R"))


################################################################################
### 1. Load InferCNV output and create heatmap and metadata dfs ###
################################################################################
  
if (file.exists(paste0(out_dir, "infercnv.12_denoised.observations.txt")) | 
  file.exists(paste0(out_dir, "infercnv.15_denoised.observations.txt"))) {

  if (!file.exists(paste0(Robject_dir, "/heatmap_df.RData")) & 
  	!file.exists(paste0(Robject_dir, "/metadata_df.RData"))) {
    # load InferCNV output:
    print("Loading InferCNV output files...")
    if (file.exists(paste0(out_dir, "infercnv.12_denoised.observations.txt"))) {
      infercnv_output <- as.data.frame(t(read.table(paste0(out_dir, 
    	  "infercnv.12_denoised.observations.txt"))))
    } else if (file.exists(paste0(out_dir, "infercnv.15_denoised.observations.txt"))) {
    	infercnv_output <- as.data.frame(t(read.table(paste0(out_dir, 
    	  "infercnv.15_denoised.observations.txt"))))
    }

    # determine the epithelial cells and only include these in heatmap:
    print(paste0("Number of heatmap rows before non-epithelial thrown: ", nrow(infercnv_output)))
    epithelial_ids <- metadata_df$cell_ids[grep("pithelial", metadata_df$cell_type)]
    heatmap_df <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
    print(paste0("Number of heatmap rows after non-epithelial thrown: ", nrow(heatmap_df)))
  
    # check all epithelial cells are present:
    if (!file.exists(paste0(input_dir, "integrated_object_epithelial_cell_ids.txt"))) {
      integrated_object <- readRDS(paste0(integrated_dir, 
        "/01_seurat_subset_Epithelial.Rdata"))
      integrated_epi_ids <- names(Idents(integrated_object))[grep(sample_name, 
        names(Idents(integrated_object)))]
      write.table(integrated_epi_ids, paste0(input_dir, 
        "integrated_object_epithelial_cell_ids.txt"), quote=F, sep="\n", col.names=F,
        row.names=F)
    } else {
      integrated_epi_ids <- read.table(paste0(input_dir, "integrated_object_epithelial_cell_ids.txt"),
        header=F, sep="\n", as.is=T)[,1]
    }
    missing_epis <- integrated_epi_ids[!(integrated_epi_ids %in% epithelial_ids)]
    print(paste0("Total no. epithelial cells in integrated data but not in InferCNV output = ",
      length(missing_epis)))
  
    # remove normal epithelial cells by throwing anything with CNA value < 0.4 or 
    #correlation value < 0.2:
    infercnv_cna_values <- read.table(paste0(ref_dir, 
    	"/infercnv_measures_mean_of_scaled_squares_df.txt"), sep = "\t", header = T)
    sample_cna_values <- infercnv_cna_values[
      grep(sample_name, rownames(infercnv_cna_values)),
    ]
    sample_cna_values <- subset(sample_cna_values, select = c(infercnv_levels, 
    	significant_infercnv_correlation_0.05))
    to_keep <- rownames(sample_cna_values)[sample_cna_values$infercnv_levels > 0.2 | 
      sample_cna_values$significant_infercnv_correlation_0.05 > 0.4]
    heatmap_df <- heatmap_df[to_keep,]

    # create cluster metadata df:
    seurat_10X <- readRDS(paste0(sample_dir, "Rdata/03_seurat_object_processed.RData"))
    Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  
  	metadata <- prepare_infercnv_metadata(seurat_10X, subset_data = F, 
  		as.data.frame(t(heatmap_df)), for_infercnv=F)
  	metadata_df <- metadata$metadata_df
  	metadata_df <- metadata_df[rownames(heatmap_df),]

    # scale heatmap df to values between -1 and 1:
    heatmap_df <- as.data.frame(rescale(as.matrix(heatmap_df), c(-1, 1)))

    saveRDS(heatmap_df, paste0(Robject_dir, "/heatmap_df.RData"))
  	saveRDS(metadata_df, paste0(Robject_dir, "/metadata_df.RData"))

  } else {
  	print("Loading heatmap and metadata dfs...")
  	heatmap_df <- readRDS(paste0(Robject_dir, "/heatmap_df.RData"))
  	metadata_df <- readRDS(paste0(Robject_dir, "/metadata_df.RData"))
  }


  ################################################################################
  ### 2. Create heatmap annotations ###
  ################################################################################

   # create group annotation df:
  if (!file.exists(paste0(Robject_dir, "group_annotation.RData"))) {
  	group_annotation <- create_group_annotation(heatmap_df, metadata_df)
  	saveRDS(group_annotation, paste0(Robject_dir, "group_annotation.RData"))
  } else {
  	group_annotation <- readRDS(paste0(Robject_dir, "group_annotation.RData"))
  }
  
  # fetch chromosome boundary co-ordinates:
  if (!file.exists(paste0(Robject_dir, "chromosome_data.RData"))) {
  	chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)
  	saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.RData"))
  } else {
  	chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.RData"))
  }
  
  # create nUMI and nGene annotations:
  if (!file.exists(paste0(Robject_dir, "QC_annotation.RData"))) {
  	if(!exists("seurat_10X")) {
  	  seurat_10X <- readRDS(paste0(sample_dir, "Rdata/03_seurat_object_processed.RData"))
  	}
  	qc_df <- data.frame(
	  seurat_10X@meta.data$nCount_RNA,
	  seurat_10X@meta.data$nFeature_RNA,
	  row.names = as.character(names(Idents(seurat_10X)))
    )
    colnames(qc_df) <-  c("nUMI", "nGene")
    QC_annotation <- create_QC_annotation(qc_df, heatmap_df)
    saveRDS(QC_annotation, paste0(Robject_dir, "QC_annotation.RData"))
  } else {
  	QC_annotation <- readRDS(paste0(Robject_dir, "QC_annotation.RData"))
  }

  # create array CNV annotation:
  all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
  if (which(colnames(all_array_CNVs) %in% sample_name) > 0) {
    array_CNV_annotation <- create_array_CNV_annotation(heatmap_df, all_array_CNVs)
    saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.RData"))
    grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
      heatmap_legend_side = "left"))
  }

  # prepare df for plotting:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))

  # define heatmap colours:
  heatmap_cols <- colorRamp2(c(-1, 0, 1), c("#00106B", "white", "#680700"), space = "sRGB")
  
  print("Generating final heatmap...")

  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = heatmap_cols,
    cluster_columns = F, cluster_rows = F,
    split = group_annotation$group_annotation_df$group,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = F,
    show_heatmap_legend = F,
    heatmap_legend_param = list(labels_gp = gpar(col = "red", fontsize = 12)),
    #heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    #grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    #title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
    use_raster = T, raster_device = c("png")
  )
  
  
  ################################################################################
  ### 3. Generate annotated heatmap###
  ################################################################################
  
  ht_list <- group_annotation$group_annotation + final_heatmap + 
  QC_annotation$nUMI_annotation + QC_annotation$nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
  )
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0037

  # save final heatmap object with measurements needed to create it:
  if (!file.exists(paste0(Robject_dir, "final_heatmap_object.RData"))) {
  	annotated_heatmap_and_measurements <- list(
    	annotated_heatmap = annotated_heatmap, 
    	chr_data = chr_data,
    	x_coord = x_coord
    )
  
    saveRDS(annotated_heatmap_and_measurements, paste0(Robject_dir, "final_heatmap_object.RData"))
  }
  
  # plot final annotated heatmap:
  pdf(paste0(plot_dir, "test156.pdf"), height = 12, width = 18)   
    grid.newpage()
 	  pushViewport(viewport(x = 0.027, y = 0.05, width = 0.95, height = 0.76, 
 	  	just = c("left", "bottom")))
        grid.draw(annotated_heatmap)
        decorate_heatmap_body("hm", {
          for ( e in 1:length(chr_data$end_pos) ) {
          grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), 
          	gp = gpar(lwd = 1, col = "#383838"))
          grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
            unit(0, "npc") + unit(-3.5, "mm"), gp=gpar(fontsize=12))
          }
        })
      popViewport()
      pushViewport(viewport(x=x_coord + 0.907, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nUMI", rot=65)
      popViewport()
      pushViewport(viewport(x=x_coord + 0.925, y=0.0001, width = 0.1, height = 0.1, 
        just = "bottom"))
        grid.text("nGene", rot=65)
      popViewport()
      pushViewport(viewport(x = 0.941, y = 0.835, 
       	width = 0.857, height = 0.16, just = c("right", "bottom")))
        grid.draw(grid_array_heatmap)
      popViewport()
      
  dev.off()
}
print(paste0("Heatmap created, output in ", plot_dir))
