#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
subproject_name <- args[1]
sample_name <- args[2]
subset_data <- as.logical(args[3])
draw_bulk_annotations <- as.logical(args[4])

#subproject_name <- "identify_normals"
#sample_name <- "CID4535"
#subset_data <- FALSE
#draw_bulk_annotations <- FALSE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Draw bulk annotations? ", as.character(draw_bulk_annotations)))

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

  if (!file.exists(paste0(Robject_dir, "/combined_heatmap_df.RData")) & 
  	!file.exists(paste0(Robject_dir, "/combined_metadata_df.RData"))) {
    # load InferCNV output:
    print("Loading InferCNV output files...")
    if (file.exists(paste0(out_dir, "infercnv.12_denoised.observations.txt"))) {
      infercnv_output <- as.data.frame(t(read.table(paste0(out_dir, 
    	  "infercnv.12_denoised.observations.txt"))))
    } else if (file.exists(paste0(out_dir, "infercnv.15_denoised.observations.txt"))) {
    	infercnv_output <- as.data.frame(t(read.table(paste0(out_dir, 
    	  "infercnv.15_denoised.observations.txt"))))
    }
  
    # load metadata:
    metadata_df <- read.table(paste0(input_dir, "metadata.txt"), header=F,
    	  sep="\t", as.is=T)
    	colnames(metadata_df) <- c("cell_ids", "cell_type")
  
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
  
    # remove normal epithelial cells by throwing anything with < 0.4 CNA value:
    infercnv_cna_values <- read.table(paste0(ref_dir, 
    	"/infercnv_measures_mean_of_scaled_squares_df.txt"), sep = "\t", header = T)
    sample_cna_values <- infercnv_cna_values[
      grep(sample_name, rownames(infercnv_cna_values)),
    ]
    sample_cna_values <- subset(sample_cna_values, select = c(infercnv_levels, 
    	significant_infercnv_correlation_0.05))
    to_keep <- rownames(sample_cna_values)[sample_cna_values$infercnv_levels > 0.2 & 
      sample_cna_values$significant_infercnv_correlation_0.05 > 0.4]
    heatmap_df <- heatmap_df[to_keep,]
    metadata_df <- metadata_df[metadata_df$cell_ids %in% to_keep,]

    # scale heatmap df to values between -1 and 1:
    heatmap_df <- as.data.frame(rescale(as.matrix(heatmap_df), c(-1, 1)))
  
    
    ################################################################################
    ### 2. Load SNP array CNV and add to heatmap and metadata dfs ###
    ################################################################################
   
    # load SNP array CNV values if they exist:
    all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
    if (any(colnames(all_array_CNVs) %in% sample_name)) {
      print("Adding SNP array CNVs to heatmap df...")
      array_CNVs <- all_array_CNVs[,colnames(all_array_CNVs) %in% sample_name]
      names(array_CNVs) <- rownames(all_array_CNVs)
  
      ## convert allele number to total copy number values and take the exponential
      # to reverse the log transformations:
      array_CNVs <- array_CNVs/2
      # array_df <- data.frame(CNV=array_CNVs, sample=rep("CID4463", length(array_CNVs)))
      # p <- ggplot(array_df, aes(x=sample, y=CNV))
      # p <- p + geom_jitter()
      # 
      # scaled_array_CNVs <- rescale(array_CNVs, c(-1,1))
      # scaled_array_df <- data.frame(CNV=scaled_array_CNVs, sample=rep("CID4463", length(scaled_array_CNVs)))
      # p <- ggplot(scaled_array_df, aes(x=sample, y=CNV))
      # p <- p + geom_jitter()
      # 
      # scaled_log2_array_CNVs <- rescale(log2_array_CNVs, c(-1,1))
      
      log2_array_CNVs <- log2(array_CNVs)
      # log2_array_df <- data.frame(CNV=log2_array_CNVs, sample=rep("CID4463", length(log2_array_CNVs)))
      # p <- ggplot(log2_array_df, aes(x=sample, y=CNV))
      # p <- p + geom_jitter()
      
      # select only genes present in heatmap_df:
      log2_array_CNVs <- log2_array_CNVs[colnames(heatmap_df)]
      print(paste0(
      	"No. genes without CNV information = ", length(which(is.na(log2_array_CNVs)))
      ))
      names(log2_array_CNVs) <- colnames(heatmap_df)

      # scale array CNVs to values between -1 and 1:
      # array_CNVs <- rescale(array_CNVs, c(-1, 1))
    
      # add array CNVs as multiple rows at start of heatmap and update metadata_df:
      no_array_rows <- round(nrow(heatmap_df)/10, 0)
      array_CNV_df <- as.data.frame(t(data.frame(log2_array_CNVs)))
      array_CNV_df <- array_CNV_df[rep(1, no_array_rows),]
      heatmap_df <- rbind(array_CNV_df, heatmap_df)
    
      array_metadata <- data.frame(cell_ids = rownames(array_CNV_df), 
      	cell_type = "SNP_array")
      metadata_df <- rbind(array_metadata, metadata_df)
    
      saveRDS(heatmap_df, paste0(Robject_dir, "/combined_heatmap_df.RData"))
      saveRDS(metadata_df, paste0(Robject_dir, "/combined_metadata_df.RData"))
    } else {
      saveRDS(heatmap_df, paste0(Robject_dir, "/combined_heatmap_df.RData"))
      saveRDS(metadata_df, paste0(Robject_dir, "/combined_metadata_df.RData"))
    }
  } else {
  	print("Loading combined InferCNV and SNP array CNV dataframes...")
  	heatmap_df <- readRDS(paste0(Robject_dir, "/combined_heatmap_df.RData"))
    metadata_df <- readRDS(paste0(Robject_dir, "/combined_metadata_df.RData"))
  }


  ################################################################################
  ### 3. Create heatmap and annotations ###
  ################################################################################

  # create group annotation df:
  if (!file.exists(paste0(Robject_dir, "combined_group_annotation.RData"))) {
  	group_annotation <- create_group_annotation(heatmap_df, metadata_df)
  	saveRDS(group_annotation, paste0(Robject_dir, "combined_group_annotation.RData"))
  } else {
  	group_annotation <- readRDS(paste0(Robject_dir, "combined_group_annotation.RData"))
  }
  
  # fetch chromosome boundary co-ordinates:
  if (!file.exists(paste0(Robject_dir, "combined_chromosome_data.RData"))) {
  	chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)
  	saveRDS(chr_data, paste0(Robject_dir, "combined_chromosome_data.RData"))
  } else {
  	chr_data <- readRDS(paste0(Robject_dir, "combined_chromosome_data.RData"))
  }
  
  # create nUMI and nGene annotations:
  if (!file.exists(paste0(Robject_dir, "combined_QC_annotation.RData"))) {
  	seurat_10X <- readRDS(paste0(sample_dir, "Rdata/03_seurat_object_processed.RData"))
  	qc_df <- data.frame(
	  seurat_10X@meta.data$nCount_RNA,
	  seurat_10X@meta.data$nFeature_RNA,
	  row.names = as.character(names(Idents(seurat_10X)))
    )
    colnames(qc_df) <-  c("nUMI", "nGene")
    QC_annotation <- create_QC_annotation(qc_df, heatmap_df)
    saveRDS(QC_annotation, paste0(Robject_dir, "combined_QC_annotation.RData"))
  } else {
  	QC_annotation <- readRDS(paste0(Robject_dir, "combined_QC_annotation.RData"))
  }

  # generate histograms of final array and inferCNV vale distributions:
  final_infercnv <- unlist(heatmap_df[grep("array", rownames(heatmap_df), invert=T),])
  pdf(paste0(plot_dir, "final_invercnv_histogram"))
    hist(final_infercnv)
  dev.off()

  if (length(grep("array", rownames(heatmap_df))) > 0) {
    final_array <- unlist(heatmap_df[grep("array", rownames(heatmap_df)),])
    pdf(paste0(plot_dir, "final_array_histogram"))
      print(hist(final_array))
    dev.off()
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
    show_row_dend = FALSE,
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
  ### 5. Plot heatmap and annotations ###
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

  # save final heatmap objects with measurements needed to create it:
  annotated_heatmap_and_measurements <- list(
  	annotated_heatmap = annotated_heatmap,
  	heatmap_df = heatmap_df,
    group_annotation = group_annotation$group_annotation_df,
    QC_annotation = QC_annotation$qc_df,
  	chr_data = chr_data, 
  	hlines = hlines, 
  	x_coord = x_coord
  )

  saveRDS(annotated_heatmap_and_measurements, paste0(Robject_dir, "final_combined_heatmap_object.RData"))

  # plot final annotated heatmap:
  pdf(paste0(plot_dir, "final_combined_infercnv_heatmap.pdf"), height = 14.5, width = 18)   
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
}
print(paste0("Heatmap created, output in ", plot_dir))
