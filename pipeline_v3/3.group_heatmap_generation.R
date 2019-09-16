# bash call of script:
subtype="primary_vs_metastasis"

#mkdir -p /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/logs/$subtype/
#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/3.group_heatmap_generation.R
#qsub -wd /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/logs/$subtype/ -pe smp 6 -N $subtype -b y -j y -V -P TumourProgression "${R} CMD BATCH  --no-save /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/3.group_heatmap_generation.R"

#mkdir -p /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/logs/$subtype/
#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/3.group_heatmap_generation.R
#qsub -wd /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/logs/$subtype/ -pe smp 6 -N $subtype -b y -j y -V -P TumourProgression "${R} CMD BATCH  --no-save '--args $a' /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/pipeline_v3/3.group_heatmap_generation.R"

library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)
library(HiddenMarkov, lib.loc=lib_loc)
library(ComplexHeatmap, lib.loc = lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(grid)
library(RColorBrewer)

#temp_args <-
#  commandArgs(trailingOnly = T)
#include_annotations <- as.logical(temp_args[1])
#print(paste0("include annotations: ", include_annotations))

## ER+
#sample_ids <- c("CID3948", "CID3941", "CID4067", "CID4290A", "CID4461", 
#  "CID4463", "CID4465", "CID4471", "CID4530", "CID4535")
#metastases <- "none"

## HER2+
#sample_ids <- c("CID3586", "CID4066", "CID4398", "CID3921", "CID3963", 
#  "CID45171", "CID45172")
#metastases <- c("CID45172")

## TNBC
#sample_ids <- c("CID44041", "CID4495", "CID44971", "CID44991", "CID44992", 
#  "CID4515", "CID43862", "CID43863", "CID44972")
#metastases <- c("CID43862", "CID43863", "CID44972")

## metastases
#sample_ids <- c("CID43862", "CID43863")
#metastases <- c("CID43862", "CID43863")

# primary_vs_metastasis
sample_ids <- c("CID45171", "CID45172")
metastases <- c("CID45172")

#sample_ids <- c("CID3921", "CID3941", "CID3948", "CID3963", 
#  "CID4066", "CID4067", "CID4290A", "CID4386", "CID43862", 
#  "CID43863", "CID4398", "CID44041", "CID44042", "CID4408", 
#  "CID4409", "CID4461", "CID4463", "CID4465", "CID4471", 
#  "CID4495", "CID44971", "CID44972", "CID44991", "CID44992", 
#  "CID4513", "CID4515", "CID45171", "CID45172", "CID4520N", 
#  "CID4523", "CID4523N", "CID4530", "CID4530N", "CID4535" )

include_annotations <- TRUE
subset <- FALSE

subproject <- "brca_mini_atlas_030719"
missing_genes_colour <- "white"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
func_dir <- paste0(project_dir, "/scripts/pipeline_v3/functions/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject, "/")

if (include_annotations) {
  out_path <- paste0(seurat_path, paste(sample_ids, collapse = "_"), "/")
} else {
  out_path <- paste0(seurat_path, paste(sample_ids, collapse = "_"), 
    "_no_annotations/")
}
ref_dir <- paste0(project_dir, "/refs/")
system(paste0("mkdir -p ", out_path))


#########################################################################################
### 0. Define functions ###
#########################################################################################

temp_png_function <- dget(paste0(func_dir, "temp_png_function.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
create_infercnv_object <- dget(paste0(func_dir, "create_infercnv_object.R"))
run_infercnv <- dget(paste0(func_dir, "run_infercnv.R"))
create_group_annotation_groupplot <- dget(paste0(func_dir, "create_group_annotation_groupplot.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
annotate_PAM50_CNV <- dget(paste0(func_dir, "annotate_PAM50_CNV.R"))
create_CNV_genes_annotation <- dget(paste0(func_dir, "create_CNV_genes_annotation.R"))
create_GIN_df <- dget(paste0(func_dir, "create_GIN_df.R"))
create_QC_annotation_group <- dget(paste0(func_dir, "create_QC_annotation_group.R"))


#########################################################################################
### 1. Load data ###
#########################################################################################

# fetch vector of all genes:
all_genes <- as.character(read.table(paste0(ref_dir, "/infercnv_gene_order.txt"))$V1)

# load samples, metadata and fetch vector of genes contained in samples:
for (i in 1:length(sample_ids)) {
  print(i)
  metadata_dir <- paste0(seurat_path, "/seurat_", sample_ids[i], 
    "/Output/InferCNV/supervised_clustering/subpop_mode/input_files/")
  infercnv_metadata <- read.table(paste0(metadata_dir, "metadata.txt"),
    header=F)
  colnames(infercnv_metadata) <- c("cell_ids", "cell_type")
  # remove CAFs from heatmap df and metadata:
  cells_to_remove <- 
  infercnv_metadata$cell_ids[grep("pithelial", infercnv_metadata$cell_type, invert=T)]
  infercnv_metadata <- infercnv_metadata[!(infercnv_metadata$cell_ids %in% cells_to_remove),]
  print(dim(infercnv_metadata))

  if ( length(grep(sample_ids[i], infercnv_metadata$cell_type)) < 1 ) {
    infercnv_metadata$cell_type <- paste0(sample_ids[i], "_", infercnv_metadata$cell_type)
  }
  
  infercnv_output_filename <- list.files(
      paste0(seurat_path, "seurat_", 
      sample_ids[i], 
      "/Output/InferCNV/supervised_clustering/subpop_mode/"), 
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt",
      full.names = T
    )

  if (length(infercnv_output_filename) != 0) {
    infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))
    infercnv_output <- infercnv_output[!(rownames(infercnv_output) %in% cells_to_remove),]
    print(dim(infercnv_output))
  
    if (subset) {
      infercnv_output <- infercnv_output[,1:300]
    }
  
    output_genes <- colnames(infercnv_output)
  
    if (i==1) {
  
      metadata_df <- infercnv_metadata
      infercnv_output_list <- list(infercnv_output)
      infercnv_output_genes <- output_genes
  
    } else {

      if (exists("metadata_df")) {
        metadata_df <- rbind(metadata_df, infercnv_metadata)
        infercnv_output_list[[i]] <- infercnv_output
        infercnv_output_genes <- unique(c(infercnv_output_genes, output_genes))
      } else {
        metadata_df <- infercnv_metadata
        infercnv_output_list <- list(infercnv_output)
        infercnv_output_genes <- output_genes
      }
    }

  } else {

    if (i==1) {
      sample_ids_to_remove <- sample_ids[i]
    } else {
      if (exists("sample_ids_to_remove")) {
        sample_ids_to_remove <- c(sample_ids_to_remove, sample_ids[i])
      } else {
        sample_ids_to_remove <- sample_ids[i]
      }
    }
  }
}

if (exists("sample_ids_to_remove")) {
  sample_ids <- sample_ids[!(sample_ids_to_remove == sample_ids)]
}

GIN_dfs <- lapply(infercnv_output_list, create_GIN_df)
all_GIN <- do.call("rbind", GIN_dfs)

o=1
create_normal_call_df <- function(df, metadata, GIN_dfs) {
  # determine top 5% of tumour cells:
  epithelial_cells <- rownames(df)
  GIN_df <- GIN_dfs[[o]]
  epithelial_GIN <- GIN_df[as.character(GIN_df$cell_id) %in% epithelial_cells,]
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
  cor_df <- cbind(cor_df, GIN_df$GIN)
  sig_cors <- cor_df[cor_df$cor.p.value < 0.05,]
  sig_cors <- data.frame(cell_id=rownames(sig_cors), cor=sig_cors$cor.estimate)
  normal_cors <- sig_cors$cell_id[sig_cors$cor < 0.3]
  normal_GIN <- epithelial_GIN$cell_id[epithelial_GIN$GIN < 0.2]
  normal_calls <- normal_cors[normal_cors %in% normal_GIN]
  normal_call_df <- data.frame(row.names = epithelial_cells, 
    call = rep("other", length(epithelial_cells)), stringsAsFactors=F)
  normal_call_df$call[rownames(normal_call_df) %in% normal_calls] <- "normal"
  
  o <<- o+1
  return(normal_call_df)
}

o=1
normal_call_dfs <- lapply(infercnv_output_list, create_normal_call_df, 
  metadata_df, GIN_dfs)
all_normal_calls <- do.call("rbind", normal_call_dfs)

# add nUMI and nGene barplot annotations:
for (i in 1:length(sample_ids)) {
  print(paste0("Loading ", sample_ids[i], " seurat object..."))
  if (i==1) {
    if ( length(grep("PDX|IND", sample_ids[i])) == 0 ){
      seurat_list <- list(readRDS(paste0(seurat_path, "seurat_", 
      sample_ids[i], "/Output/Rdata/03_seurat_object_processed.RData")))
    } else {
      seurat_list <- list(readRDS(paste0(seurat_path, "seurat_", 
      sample_ids[i], "/Output/Rdata/04_seurat_object_combined.RData")))
    }
  } else {
    if ( length(grep("PDX|IND", sample_ids[i])) == 0 ){
      seurat_list[[i]] <- readRDS(paste0(seurat_path, "seurat_", 
      sample_ids[i], "/Output/Rdata/03_seurat_object_processed.RData"))
    } else {
      seurat_list[[i]] <- readRDS(paste0(seurat_path, "seurat_", 
      sample_ids[i], "/Output/Rdata/04_seurat_object_combined.RData"))
    }
  }
}

p=1
QC_dfs <- lapply(infercnv_output_list, create_QC_annotation_group, seurat_list)
QC_df <- do.call("rbind", QC_dfs)

infercnv_output_dfs <- lapply(infercnv_output_list, function(x) {

  # create 'missing_genes' dataframe of NA values:
  missing_genes <- infercnv_output_genes[!(infercnv_output_genes %in% colnames(x))]
  
  if ( length(missing_genes) > 0 ) {

    missing_genes_list <- c(
      rep(
        list( rep(NA, nrow(x)) ), length(missing_genes)
      )
    )
    missing_genes_df <- do.call("cbind", missing_genes_list)
    
    rownames(missing_genes_df) <- rownames(x)
    colnames(missing_genes_df) <- missing_genes

    # cbind to infercnv output heatmap:
    non_ordered_result_df <- cbind(x, missing_genes_df)
  
    # order columns:
    common_genes <- all_genes[all_genes %in% colnames(non_ordered_result_df)]
    m <- match(common_genes, colnames(non_ordered_result_df))
    ordered_result_df <- non_ordered_result_df[,m]
    print(paste0("Checking if columns of extended data frame are in correct order... ", 
      identical(colnames(ordered_result_df), common_genes)))
  
    return(ordered_result_df)

  } else {
    return(x)
  }
})
heatmap_df <- do.call("rbind", infercnv_output_dfs)

# create group annotation df
group_annotation <- create_group_annotation_groupplot(heatmap_df, metadata_df, metastases)
# ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
heatmap_df <- heatmap_df[m,]

m <- match(rownames(heatmap_df), rownames(all_GIN))
all_GIN <- all_GIN[m,]
GIN_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    all_GIN$GIN,
    gp = gpar(
      col = "#AF548E", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)

m <- match(rownames(heatmap_df), rownames(all_normal_calls))
all_normal_calls <- data.frame(calls = all_normal_calls[m,], 
  row.names = rownames(all_normal_calls)[m])
normal_call_annotation <- Heatmap(all_normal_calls, 
      col = c("other" = "#E7E4D3", "normal" = "#1B7837"), 
      name = "normal_call_annotation", width = unit(6, "mm"), 
      show_row_names = F, show_column_names = F, 
      heatmap_legend_param = list(title = "Normal epithelial calls", title_gp = gpar(fontsize = 8, 
      fontface = "bold"), labels_gp = gpar(fontsize = 6)))

m <- match(rownames(heatmap_df), rownames(QC_df))
QC_df <- QC_df[m,]
nUMI_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      QC_df$nUMI,
      gp = gpar(
        col = "#D8B72E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )

nGene_annotation <- rowAnnotation(
  correlation_annotation = anno_barplot(
    QC_df$nGene,
    gp = gpar(
      col = "#9ECAE1", 
      width = unit(4, "cm")
    ), 
    border = FALSE, 
    which = "row", 
    axis = F
  )
)

# define chromosome lengths in terms of number of filtered genes
chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)

# determine proportion of genes used for plot:
gene_proportion <- paste0(
  as.character(
    round(
      nrow(heatmap_df)/length(all_genes)*100, 2
    )
  ), "%"
)

if (!include_annotations) {

  # replace column labels with two letter label for positioning:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  print("Generating final heatmap...")
  
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = colorRamp2(c(0, 0.5, 1, 1.5, 2, 3), 
      c("#00106B", "#9191CC", "white", "#DDB6B6", "#AB4848", "#930707"), 
      space = "sRGB"),  
    na_col = missing_genes_colour, 
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
  for ( n in 1:(length(spl_groups)-1) ) {
    if (n==1) {
      hlines <- c(length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group))
    } else {
      hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group)
    }
  }
  hlines <- 1-hlines

  ht_list <- group_annotation$type_annotation + group_annotation$sample_annotation + 
    group_annotation$group_annotation + 
    final_heatmap + normal_call_annotation + GIN_annotation + 
    nUMI_annotation + nGene_annotation
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0045
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, heatmap_legend_side = "left")
  )
  
  if (subset) {
    pdf(paste0(out_path, "subset_infercnv_prediction_heatmap_non_annotated.pdf"), 
      height = 14.5, width = 17)
  } else {
    pdf(paste0(out_path, "final_infercnv_prediction_heatmap_non_annotated.pdf"), 
      height = 14.5, width = 17)
  }
    # plot heatmap:
    pushViewport(viewport(x = 0.5, y = 0.1, 
                          width = 1, height = 0.8, just = "bottom"))
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

    pushViewport(viewport(x=x_coord + 0.824, y=0.038, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("normal calls", rot=55)
    popViewport()

    pushViewport(viewport(x=x_coord + 0.854, y=0.05, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("GIN", rot=55)
    popViewport()

    pushViewport(viewport(x=x_coord + 0.871, y=0.05, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nUMI", rot=55)
    popViewport()

    pushViewport(viewport(x=x_coord + 0.887, y=0.05, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nGene", rot=55)
    popViewport()

    pushViewport(viewport(x=0.036, y=0.25, width = 0.2, height = 0.1, just = "bottom"))
    #grid.draw(lollipop)
    grid.text("% genes used:", gp = gpar(fontface = "bold"))
    popViewport()

    pushViewport(viewport(x=0.025, y=0.23, width = 0.2, height = 0.1, just = "bottom"))
    #grid.draw(lollipop)
    grid.text(gene_proportion)
    popViewport()

  dev.off()

} else {
  

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
  
  # replace column labels with two letter label for positioning:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  print("Generating final heatmap...")
  
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]

  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = colorRamp2(c(0, 0.5, 1, 1.5, 2, 3), 
      c("#00106B", "#9191CC", "white", "#DDB6B6", "#AB4848", "#930707"), 
      space = "sRGB"), 
    na_col = missing_genes_colour, 
    cluster_columns = F, cluster_rows = F,
    split = group_annotation$group_annotation_df$group,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = FALSE,
    bottom_annotation = CNV_genes_annotation, bottom_annotation_height = unit(2, "cm"),
    gap = unit(1, "cm"),
    heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), grid_width = unit(1.5, "cm"), legend_direction = "horizontal",
    title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
    use_raster = T, raster_device = c("png")
  )

  
  # determine co-ordinates of horizontal lines at group borders:
  spl_groups <- split(group_annotation$group_annotation_df$group, 
  group_annotation$group_annotation_df$group)
  spl_groups <- spl_groups[unique(group_annotation$group_annotation_df$group)]
  for ( n in 1:(length(spl_groups)-1) ) {
    if (n==1) {
      hlines <- c(length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group))
    } else {
      hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group)
    }
  }
  hlines <- 1-hlines

  ht_list <- group_annotation$type_annotation + group_annotation$sample_annotation + 
    group_annotation$group_annotation + 
    final_heatmap + normal_call_annotation + GIN_annotation + 
    nUMI_annotation + nGene_annotation
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0045
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, heatmap_legend_side = "left")
  )
  
  if (subset) {
    pdf(paste0(out_path, "subset_infercnv_prediction_heatmap.pdf"), 
        height = 14.5, width = 17)
  } else {
      pdf(paste0(out_path, "final_infercnv_prediction_heatmap.pdf"), 
        height = 14.5, width = 17)
  }
    grid.newpage()
    # plot Normal subtype:
    pushViewport(viewport(x = x_coord+0.003, y = 0.090,
                          width = 0.845+0.019, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[5]])
    popViewport()
    
    # plot Basal subtype:
    pushViewport(viewport(x = x_coord+0.008, y = 0.169,
                          width = 0.845+0.014, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[4]])
    popViewport()
    
    # plot Her2 subtype:
    pushViewport(viewport(x = x_coord+0.01, y = 0.251, 
                          width = 0.845+0.012, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[3]])
    popViewport()
    
    # plot LumB subtype:
    pushViewport(viewport(x = x_coord+0.007, y = 0.333, 
                          width = 0.845+0.015, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[2]])
    popViewport()
    
    # plot LumA subtype:
    pushViewport(viewport(x = x_coord+0.007, y = 0.411, 
                          width = 0.845+0.015, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[1]])
    popViewport()
    pushViewport(viewport(x = 0, y = 0.4, 
                          width = 0.98, height = 0.6, just = c("left", "bottom")))
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

    pushViewport(viewport(x=x_coord + 0.885, y=0.045, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("GIN", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.905, y=0.045, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nUMI", rot=65)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.92, y=0.045, width = 0.1, height = 0.1, just = "bottom"))
    grid.text("nGene", rot=65)
    popViewport()
    
  dev.off()

}


# convert pdf to png:
system(paste0("for p in ", out_path, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
  "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
  "convert -density 150 ", out_path, "$f -quality 90 ", out_path, "$new; done"))

print(paste0("Group heatmap created, output in ", out_path))
