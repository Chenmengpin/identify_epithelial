#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
subproject_name <- args[1]
sample_name <- args[2]
subset_data <- as.logical(args[3])
numcores=as.numeric(args[4])
draw_bulk_annotations <- as.logical(args[5])

#subproject_name <- "identify_normals"
#sample_name <- "CID4513"
#subset_data <- FALSE
#numcores=6
#draw_bulk_annotations <- FALSE

print(paste0("Subproject name = ", subproject_name))
print(paste0("Sample name = ", sample_name))
print(paste0("Subset data? ", as.character(subset_data)))
print(paste0("Number cores = ", numcores))
print(paste0("Draw bulk annotations? ", as.character(draw_bulk_annotations)))

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(cluster, lib.loc = lib_loc)
library(Seurat)
library(infercnv, lib.loc=lib_loc)
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
create_infercnv_object <- dget(paste0(func_dir, "create_infercnv_object.R"))
run_infercnv <- dget(paste0(func_dir, "run_infercnv.R"))
create_group_annotation <- dget(paste0(func_dir, "create_group_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))
create_array_CNV_annotation <- dget(paste0(func_dir, "create_array_CNV_annotation.R"))


################################################################################
### 1. Generate input matrix and metadata files ###
################################################################################

# set up subset output directory if necessary:
if (subset_data) {
  out_dir <- paste0(out_dir, "subset")
  input_dir <- paste0(out_dir, "/input_files/")
  system(paste0("mkdir -p ", input_dir))
}

if (!file.exists(paste0(out_dir, "infercnv.12_denoised.observations.txt")) & 
	!file.exists(paste0(out_dir, "infercnv.15_denoised.observations.txt"))) {

  # load seurat object:
  seurat_10X <- readRDS(paste0(sample_dir, "Rdata/03_seurat_object_processed.RData"))
  Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1
  
  # create raw matrix input file and subset if necessary:
  count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
  if (subset_data) {
    count_df <- count_df[1:500, 1:500]
  }
  
  ###### prepare infercnv metadata with annotated cell types and update seurat object
  # adjust to include metaplastic epithelial cells ###
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, subset_data=subset_data, 
    count_df)
  seurat_10X <- infercnv_metadata$seurat
  print(paste0("Cell types are: ", unique(infercnv_metadata$metadata$cell_type)))
  
  # only keep cells in metadata df:
  print(paste0("No cells in count df before filtering for those in metadata df = ", 
      ncol(count_df)))
  count_df <- count_df[,colnames(count_df) %in% infercnv_metadata$metadata$cell_ids]
  print(paste0("No cells in count df after filtering for those in metadata df = ", 
      ncol(count_df)))
  
  # create raw matrix input file
  if (!file.exists(paste0(input_dir, "input_matrix.txt"))) {
    print("Creating inferCNV raw counts file...")
    write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
    sep="\t", col.names=T, row.names=T)
  }
  
  # generate cluster metric plots for epithelial cluster:
  epithelial_clusters <- grep("pithelial", unique(infercnv_metadata$metadata$cell_type), value=T)
  print(paste0("Epithelial cluster = ", epithelial_clusters))
  if (!file.exists(paste0(plot_dir, "metrics_by_epithelial_cluster.png"))) {
    temp_png_function(paste0(plot_dir, "metrics_by_epithelial_cluster.png"))
      temp_violinplot <- VlnPlot(
        object = seurat_10X,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        pt.size = 1.5,
        idents = epithelial_clusters
      )
      print(temp_violinplot)
    dev.off()
  }
  
  # remove cluster information for epithelial cells:
  infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)] <- 
  gsub("_[0-9].*$", "", 
    infercnv_metadata$metadata$cell_type[grep("pithelial", infercnv_metadata$metadata$cell_type)])
  
  # check all epithelial cells from integrated data are present:
  individual_epi_ids <- names(Idents(seurat_10X))[grep("pithelial", Idents(seurat_10X))]
  
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
  
  missing_epis <- integrated_epi_ids[!(integrated_epi_ids %in% individual_epi_ids)]
  print(paste0("Total no. epithelial cells in integrated data but not individual = ",
    length(missing_epis)))
  
  # if no epithelial or CAF clusters present, abort:
  if (length(epithelial_clusters) < 1) {
    print(paste0("No epithelial/myoepithelial clusters detected, aborting..."))
  } else {
    # write metadata files and new seurat object:
    if (!file.exists(paste0(input_dir, "metadata.txt"))) {
      write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
        quote=F, sep="\t", col.names=F, row.names=F)
      write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
        "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
      saveRDS(seurat_10X, paste0(sample_dir, "Rdata/04_seurat_object_annotated.RData"))
    }
    
    # define normals which will act as InferCNV reference cells:
    normals <- grep(
    	"[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned|[u,U]nknown|[t,T]umour|[t,T]umor", 
      unique(infercnv_metadata$metadata$cell_type[
        infercnv_metadata$metadata$cell_ids %in% colnames(count_df)
      ]), value=T, 
      invert=T
    )
    
    print(paste0("Normal is: ", normals))
  
    print("Creating inferCNV object...")
    raw_path <- paste0(input_dir, "input_matrix.txt")
    annotation_path <- paste0(input_dir, "metadata.txt")
    gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
    initial_infercnv_object <- create_infercnv_object(raw_path, annotation_path,
      gene_path, normals)
    
    print("InferCNV object created, running inferCNV...")
    system.time(
    infercnv_output <- run_infercnv(initial_infercnv_object, numcores-1, out_dir, 
      0.1, 101, 3, 1.3)
    )
  }
} else {
  print("InferCNV output already exists!")

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
  ######
}

  
################################################################################
### 4. Create heatmap and annotations ###
################################################################################
  
if (file.exists(paste0(out_dir, "infercnv.12_denoised.observations.txt")) | 
  file.exists(paste0(out_dir, "infercnv.15_denoised.observations.txt"))) {
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
  if (!exists("infercnv_metadata")) {
  	metadata_df <- read.table(paste0(input_dir, "metadata.txt"), header=F,
  	  sep="\t", as.is=T)
  	colnames(metadata_df) <- c("cell_ids", "cell_type")
  } else {
  	metadata_df <- infercnv_metadata$metadata
  }

  #### load seurat object - isolate nUMI and nGene information from this so
  # it doesn't need loading each time ###
  seurat_10X <- readRDS(paste0(sample_dir, "Rdata/03_seurat_object_processed.RData"))

  # determine the epithelial cells and only include these in heatmap:
  epithelial_ids <- metadata_df$cell_ids[grep("pithelial", metadata_df$cell_type)]
  heatmap_df <- infercnv_output[rownames(infercnv_output) %in% epithelial_ids,]
  print(dim(heatmap_df))

  # check all epithelial cells are present:
  if (!exists("integrated_epi_ids")) {
  	integrated_epi_ids <- read.table(paste0(input_dir, "integrated_object_epithelial_cell_ids.txt"),
      header=F, sep="\n", as.is=T)[,1]
  }
  missing_epis <- integrated_epi_ids[!(integrated_epi_ids %in% epithelial_ids)]
  print(paste0("Total no. epithelial cells in integrated data but not in InferCNV output = ",
    length(missing_epis)))

  # create group annotation df
  group_annotation <- create_group_annotation(heatmap_df, metadata_df)
  saveRDS(group_annotation, paste0(Robject_dir, "group_annotation.RData"))
  
  # fetch chromosome boundary co-ordinates:
  chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)
  saveRDS(chr_data, paste0(Robject_dir, "chromosome_data.RData"))
  
  # create nUMI and nGene annotations:
  QC_annotation <- create_QC_annotation(seurat_10X, heatmap_df)
  saveRDS(QC_annotation, paste0(Robject_dir, "QC_annotation.RData"))

  # create array CNV annotation:
  all_array_CNVs <- read.table(paste0(ref_dir, "all_array_CNVs.txt"))
  array_CNV_annotation <- create_array_CNV_annotation(heatmap_df, all_array_CNVs)
  saveRDS(array_CNV_annotation, paste0(Robject_dir, "array_CNV_annotation.RData"))
  grid_array_heatmap <- grid.grabExpr(draw(array_CNV_annotation$array_CNV_heatmap, 
    heatmap_legend_side = "left"))

  ###### create array PDX CNV annotation - inbuild this into function ######
  pdx_name <- gsub("CID", "PDX", sample_name)
  PDX_CNVs <- all_array_CNVs[,colnames(all_array_CNVs) %in% pdx_name]
  names(PDX_CNVs) <- rownames(all_array_CNVs)

  PDX_CNVs[PDX_CNVs > 6] <- 6
  PDX_CNVs <- PDX_CNVs/2
  PDX_CNVs <- PDX_CNVs[colnames(heatmap_df)]
  PDX_CNVs[is.na(PDX_CNVs)] <- 2
  names(PDX_CNVs) <- colnames(heatmap_df)
  PDX_CNVs <- as.data.frame(PDX_CNVs)
  PDX_CNVs <- as.data.frame(t(PDX_CNVs))

  saveRDS(PDX_CNVs, paste0(Robject_dir, "PDX_CNVs.RData"))

  PDX_heatmap <- Heatmap(
  	PDX_CNVs, name=paste0("array_heatmap"),
  	col = colorRamp2(c(min(PDX_CNVs[1,]), 1, max(PDX_CNVs[1,])), 
      c("#00106B", "white", "#680700"), space = "sRGB"),
  	  cluster_columns = F, cluster_rows = F, show_row_dend = FALSE,
  	  show_row_names = F, show_column_names = F,
  	  heatmap_legend_param = list(title = "SNP array\nPDX CNV", legend_direction = "horizontal"),
  	  use_raster = T, raster_device = c("png")
  )
  grid_PDX_heatmap <- grid.grabExpr(draw(PDX_heatmap, heatmap_legend_side = "left"))

  # prepare df for plotting:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))

  saveRDS(plot_object, paste0(Robject_dir, "final_heatmap_object.RData"))
  
  print("Generating final heatmap...")

  ######
#  group_annotation <- readRDS(paste0(Robject_dir, "group_annotation.RData"))
#  chr_data <- readRDS(paste0(Robject_dir, "chromosome_data.RData"))
#  QC_annotation <- readRDS(paste0(Robject_dir, "QC_annotation.RData"))
#  array_CNVs <- readRDS(paste0(Robject_dir, "array_CNVs.RData"))
#  plot_object <- readRDS(paste0(Robject_dir, "final_heatmap_object.RData"))
  ######
  
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
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0037

  ht_list <- group_annotation$group_annotation + final_heatmap + 
  QC_annotation$nUMI_annotation + QC_annotation$nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
  )
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0037

  ###### add if PDX exists for PDX plots! ######
  ###### find out how to vertically plot complex heatmaps - may need to update package? ######
  pdf(paste0(plot_dir, "infercnv_heatmap_with_array_and_PDX_CNV_annotation.pdf"), height = 14.5, width = 18)   
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

      pushViewport(viewport(x = 0.0245, y = 0.097, 
      	width = 0.93, height = 0.05, just = c("left", "bottom")))
        grid.draw(grid_array_heatmap)
      popViewport()

      pushViewport(viewport(x = 0.0245, y = 0.047, 
      	width = 0.93, height = 0.05, just = c("left", "bottom")))
        grid.draw(grid_PDX_heatmap)
      popViewport()
      
  dev.off()
}
print(paste0("Group heatmap created, output in ", plot_dir))








  ### why is the range different than previously for the following vectors? ###
  
  # calculate CNA values as the mean of the scaled squares of all CNV values across the genome for each cell:
  CNA_values <- apply(heatmap_df, 1, function(y) {
    y[is.na(y)] <- 0
    #return(round(sum(abs(y))/ncol(x), 6))
    #return(mean(y^2))
    scaled_y <- rescale(y, c(-1, 1))
    return(mean(scaled_y^2))
  })
    



    infercnv_measures <- readRDS(paste0(seurat_path, "infercnv_measures_mean_of_scaled_squares.RData"))
    infercnv_measures <- round(eval(parse(text=paste0("infercnv_measures$", sample_name))), 6)
  
    print(paste0("Do all epithelial cells have infercnv values? ", 
    length(grep(sample_name, integrated_epi, value=T)) == nrow(infercnv_measures)))
  
    # create density plot of infercnv values:
    density_plot <- density(infercnv_measures$infercnv_levels, bw="SJ")
    
    #  prepare data for quad plots:
    infercnv_level_vs_correlation <- subset(infercnv_measures, select = c(infercnv_levels, significant_infercnv_correlation_0.05))
    
    # scale data and visualise:
    scaled_infercnv_level_vs_correlation <- scale(infercnv_level_vs_correlation) %>% as.data.frame()
  
    # run silhouette cluster analysis to determine clusters and thresholds:
    pamk_result <- pamk(scaled_infercnv_level_vs_correlation, krange=1:4)
    pamk_result$nc
    silhouette_result <- pam(scaled_infercnv_level_vs_correlation, 
                             pamk_result$nc)
    
    # use cell values passing silhouette cutoff of 0.3, 0.4, 0.5 for intercept designation:
    scuts <- c(0.4)
    for (l in 1:length(scuts)) {
      print(l)
      if (pamk_result$nc > 1) {
        sil_values <- as.data.frame(silhouette_result$silinfo$widths)
        #barplot(sil_values$sil_width)
    
        infercnv_level_vs_correlation <- merge(
          infercnv_measures,
          data.frame(row.names=names(silhouette_result$clustering),
                     cluster=silhouette_result$clustering),
          by="row.names"
        )
        
        # determine normal cluster:
        cluster_split <- split(infercnv_level_vs_correlation, 
                               infercnv_level_vs_correlation$cluster)
        cluster_max_levels <- lapply(cluster_split, function(x) max(x$significant_infercnv_correlation_0.05))
        max_vals <- do.call("c", cluster_max_levels)
        normal_cluster <- which(max_vals == min(max_vals))
      
        normal_ids <- rownames(sil_values)[sil_values$cluster == normal_cluster]
        normal_sil <- sil_values$sil_width[rownames(sil_values) %in% normal_ids]
        names(normal_sil) <- normal_ids
        normal_non_outliers <- normal_sil[normal_sil > scuts[l]]
        max_normal_infercnv_level <- max(infercnv_level_vs_correlation$infercnv_levels[
          infercnv_level_vs_correlation$Row.names %in% names(normal_non_outliers)
          ])
        max_normal_cor <- max(infercnv_level_vs_correlation$significant_infercnv_correlation_0.05[
          infercnv_level_vs_correlation$Row.names %in% names(normal_non_outliers)
          ])
          
        # define min values of cancer cells with > cutoff silhouette score:
        cancer_sil <- sil_values$sil_width[!(rownames(sil_values) %in% normal_ids)]
        names(cancer_sil) <- rownames(sil_values)[sil_values$cluster != normal_cluster]
        cancer_non_outliers <- cancer_sil[cancer_sil > scuts[l]]
        min_cancer_infercnv_level <- min(infercnv_level_vs_correlation$infercnv_levels[
          infercnv_level_vs_correlation$Row.names %in% names(cancer_non_outliers)
          ])
        min_cancer_cor <- min(infercnv_level_vs_correlation$significant_infercnv_correlation_0.05[
          infercnv_level_vs_correlation$Row.names %in% names(cancer_non_outliers)
          ])
        
        # define quad separators as values halfway between max normal and min cancer values:
        x_int <- max(c(max_normal_infercnv_level, min_cancer_infercnv_level)) - abs(max_normal_infercnv_level -                                                                               min_cancer_infercnv_level)
        y_int <- max(c(max_normal_cor, min_cancer_cor)) - abs(max_normal_cor - min_cancer_cor)
      
        # define normal and cancer cells:
        normal_call_df <- infercnv_level_vs_correlation
        normal_call_df$normal_call <- "cancer"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels < x_int & normal_call_df$significant_infercnv_correlation_0.05 < y_int
        ] <- "normal"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels < x_int & normal_call_df$significant_infercnv_correlation_0.05 > y_int
        ] <- "unassigned"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels > x_int & normal_call_df$significant_infercnv_correlation_0.05 < y_int
        ] <- "unassigned"
        rownames(normal_call_df) <- normal_call_df$Row.names
      
        if (pamk_result$nc == 4) {
  
          # create quad plot:
            p <- ggplot(normal_call_df, 
                        aes(x=infercnv_levels, y=significant_infercnv_correlation_0.05, color=as.factor(normal_call)))
            p <- p + geom_point()
            p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                      labels=c("Cancer", "Normal", "Unassigned"))
            p <- p + xlab("Infercnv level")
            p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
            p <- p + theme(legend.title = element_blank())
            p <- p + geom_vline(xintercept = x_int)
            p <- p + geom_hline(yintercept = y_int)
            p
            quad_plot <- p
          
        } else if (pamk_result$nc == 3) {
            
            # create quad plot:
            p <- ggplot(normal_call_df, 
                        aes(x=infercnv_levels, y=significant_infercnv_correlation_0.05, color=as.factor(normal_call)))
            p <- p + geom_point()
            p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                      labels=c("Cancer", "Normal", "Unassigned"))
            p <- p + xlab("Infercnv level")
            p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
            p <- p + theme(legend.title = element_blank())
            p <- p + geom_vline(xintercept = x_int)
            p <- p + geom_hline(yintercept = y_int)
            p
            quad_plot <- p
            
          } else if (pamk_result$nc == 2) {
            
            # create quad plot:
            p <- ggplot(normal_call_df, 
                        aes(x=infercnv_levels, y=significant_infercnv_correlation_0.05, color=as.factor(normal_call)))
            p <- p + geom_point()
            p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                      labels=c("Cancer", "Normal", "Unassigned"))
            p <- p + xlab("Infercnv level")
            p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
            p <- p + theme(legend.title = element_blank())
            p <- p + geom_vline(xintercept = x_int)
            p <- p + geom_hline(yintercept = y_int)
            p
            quad_plot <- p
            
          }
          system(paste0("echo ", sample_name, "_x_y_int:,", round(x_int, 3), ",", round(y_int, 3), " >> ", seurat_path, "collated_x_and_y_ints.txt"))
          
          png(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_", scuts[l], 
            "_silhouette_cutoff.png"), width = 430, height = 200)
            print(quad_plot)
          dev.off()
  
      } else if (pamk_result$nc == 1) {
        # establish cutoffs as mean cutoffs of all samples which did cluster:
        x_int <- 0.2
        y_int <- 0.35
    
        # define normal and cancer cells:
        normal_call_df <- infercnv_level_vs_correlation
        normal_call_df$normal_call <- "cancer"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels < x_int & normal_call_df$significant_infercnv_correlation_0.05 < y_int
        ] <- "normal"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels < x_int & normal_call_df$significant_infercnv_correlation_0.05 > y_int
        ] <- "unassigned"
        normal_call_df$normal_call[
          normal_call_df$infercnv_levels > x_int & normal_call_df$significant_infercnv_correlation_0.05 < y_int
        ] <- "unassigned"
        
        # create quad plot:
        p <- ggplot(normal_call_df, 
                    aes(x=infercnv_levels, y=significant_infercnv_correlation_0.05, color=as.factor(normal_call)))
        p <- p + geom_point()
        p <- p + scale_color_manual(values=c("black", "#74add1", "#b2182b"), 
                                      labels=c("Cancer", "Normal", "Unassigned"))
        p <- p + xlab("Infercnv level")
        p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
        p <- p + theme(legend.title = element_blank())
        p <- p + geom_vline(xintercept = x_int)
        p <- p + geom_hline(yintercept = y_int)
        p
        quad_plot <- p
  
        if (!file.exists(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_fixed_thresholds.png"))) {
          png(paste0(plot_dir, "normal_call_quad_plot_mean_of_scaled_squares_fixed_thresholds.png"), width = 430, height = 200)
            print(quad_plot)
          dev.off()
        }
      }
      
      if (!file.exists(paste0(plot_dir, "infercnv_level_mean_of_scaled_squares_distributions.png"))) {
        png(paste0(plot_dir, "infercnv_level_mean_of_scaled_squares_distributions.png"), width = 860, height = 400)
          par(mfrow=c(1,2))
          hist(infercnv_measures$infercnv_levels, main = NULL, xlab = "InferCNV levels")
          plot(density_plot, main=NA, xlab = "Infercnv value")
        dev.off()
      }
  
      if (pamk_result$nc > 1) {
        if (!file.exists(paste0(plot_dir, "cluster_silhouette_scores.png"))) {
          png(paste0(plot_dir, "cluster_silhouette_scores.png"), width = 430, height = 200)
            barplot(sil_values$sil_width)
          dev.off()
        }
      }
      
      # create infercnv level annotation:
      heatmap_df <- heatmap_df[rownames(heatmap_df) %in% rownames(normal_call_df),]
      m <- match(rownames(heatmap_df), rownames(normal_call_df))
      normal_call_df <- normal_call_df[m,]
      infercnv_level_annot_vector <- normal_call_df$infercnv_levels
      names(infercnv_level_annot_vector) <- rownames(normal_call_df)
  
      if (l==1) {
        final_normal_call_df <- normal_call_df
        colnames(final_normal_call_df) <- gsub(
          "normal_call", paste0("normal_call_", scuts[l], "_silhouette_cutoff"),
          colnames(final_normal_call_df)
        )
      } else {
        final_normal_call_df <- merge(
          final_normal_call_df, 
          data.frame(rownames=rownames(normal_call_df), normal_call_df$normal_call))
        colnames(final_normal_call_df[ncol(final_normal_call_df)]) <- paste0(
          "normal_call_", scuts[l], "_silhouette_cutoff"
        )
      }
    
      infercnv_level_annotation <- rowAnnotation(
        infercnv_level_annotation = anno_barplot(
          normal_call_df$infercnv_levels,
          gp = gpar(
            col = "#AF548E", 
            width = unit(4, "cm")
          ), 
          border = FALSE, 
          which = "row", 
          axis = F
        )
      )
    
      correlation_annotation <- rowAnnotation(
        correlation_annotation = anno_barplot(
          normal_call_df$significant_infercnv_correlation_0.05,
          gp = gpar(
            col = "#5628CE", 
            width = unit(4, "cm")
          ), 
          border = FALSE, 
          which = "row", 
          axis = F
        )
      )
    
      # create normal call annotation:
      normal_call_annot_df <- data.frame(row.names = rownames(normal_call_df), normal_call = normal_call_df$normal_call)
      normal_call_annot_df$normal_call <- as.character(normal_call_annot_df$normal_call)
      normal_call_annotation <- Heatmap(normal_call_annot_df, 
        col = c("unassigned" = "#E7E4D3", "normal" = "#1B7837", "cancer" = "#E7298A"), 
        name = "normal_call_annotation", width = unit(6, "mm"), 
        show_row_names = F, show_column_names = F, 
        heatmap_legend_param = list(title = "Normal epithelial calls", title_gp = gpar(fontsize = 8, 
        fontface = "bold"), labels_gp = gpar(fontsize = 6)))
  		
  	  # create nUMI and nGene annotations:
      QC_annotation <- create_QC_annotation(seurat_10X, heatmap_df)

      # prepare df for heatmap:
      plot_object <- heatmap_df
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
      
      ht_list <- group_annotation$group_annotation + final_heatmap + 
      normal_call_annotation + infercnv_level_annotation + QC_annotation$nUMI_annotation + 
      QC_annotation$nGene_annotation
      
      annotated_heatmap <- grid.grabExpr(
        draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
      )
      
      # determine where starting co-ordinates for heatmap are based upon longest cluster name
      # (0.00604 units per character):
      longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
      x_coord <- longest_cluster_name*0.0037



    
      if (draw_bulk_annotations) {
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
      }

      plot_object <- heatmap_df
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
      
      ht_list <- group_annotation$group_annotation + final_heatmap + 
      normal_call_annotation + infercnv_level_annotation + QC_annotation$nUMI_annotation + 
      QC_annotation$nGene_annotation
      
      annotated_heatmap <- grid.grabExpr(
        draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
      )
      
      # determine where starting co-ordinates for heatmap are based upon longest cluster name
      # (0.00604 units per character):
      longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
      x_coord <- longest_cluster_name*0.0037
      
      
      if (draw_bulk_annotations) {
        pdf(paste0(plot_dir, "final_infercnv_heatmap_all_annotated_", scuts[l], "_silhouette_cutoff.pdf"), 
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
          grid.text("nUMI", rot=55)
          popViewport()
          pushViewport(viewport(x=x_coord + 0.895, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
          #grid.draw(lollipop)
          grid.text("nGene", rot=55)
          popViewport()
          
        dev.off()
    
      } else {
    
        pdf(paste0(plot_dir, "final_infercnv_heatmap_normals_annotated_", scuts[l], 
          "_silhouette_cutoff.pdf"), 
          height = 14.5, width = 18)
        
          grid.newpage()
          pushViewport(viewport(x = 0, y = 0.1, 
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
          pushViewport(viewport(x=x_coord + 0.865, y=0.1, width = 0.1, height = 0.1, just = "bottom"))
          #grid.draw(lollipop)
          grid.text("nUMI", rot=55)
          popViewport()
          pushViewport(viewport(x=x_coord + 0.895, y=0.1, width = 0.1, height = 0.1, just = "bottom"))
          #grid.draw(lollipop)
          grid.text("nGene", rot=55)
          popViewport()
          
        dev.off()
      }
      print(paste0("Group heatmap created, output in ", plot_dir))
    }
  }
  print(paste0("Original dimensions: ", dim(infercnv_measures)))
  print(paste0("Final dimensions: ", dim(final_normal_call_df)))
  dims <- data.frame(dim(infercnv_measures), dim(final_normal_call_df))
  write.table(dims, paste0(table_dir, "dims_df.txt"),
    quote = F, sep = "\t")
  
  write.table(final_normal_call_df, paste0(table_dir, "infercnv_scaled_squared_mean_normal_call_values.txt"),
    quote = F, sep = "\t")
  
  # convert pdf to png:
  system(paste0("for p in ", plot_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
  "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
  "convert -density 150 ", plot_dir, "$f -quality 90 ", plot_dir, "$new; done"))

}

}

