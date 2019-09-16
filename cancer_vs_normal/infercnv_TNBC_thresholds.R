#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap, lib.loc=lib_loc)
library(circlize, lib.loc = lib_loc)
library(reshape2)
library(ggplot2)

sample_ids <- c("CID44041", "CID44971", "CID44972", "CID44991", 
  "CID44992", "CID4515", "CID43862", "CID43863" )
normal_clusters <- list(c("Epithelial_7"), c("Epithelial_18"))

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
in_path <- paste0(project_dir, "results/seurat/brca_mini_atlas_030719/")
func_dir <- paste0(project_dir, "/scripts/normal_threshold/functions/")
ref_dir <- paste0(project_dir, "/refs/")


################################################################################
### 0. Define functions ###
################################################################################

temp_png_function <- dget(paste0(func_dir, "temp_png_function.R"))
prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
create_group_annotation <- dget(paste0(func_dir, "create_group_annotation.R"))
fetch_chromosome_boundaries <- dget(paste0(func_dir, "fetch_chromosome_boundaries.R"))
gg_color_hue <- dget(paste0(func_dir, "gg_color_hue.R"))
annotate_PAM50_CNV <- dget(paste0(func_dir, "annotate_PAM50_CNV.R"))
create_CNV_genes_annotation <- dget(paste0(func_dir, "create_CNV_genes_annotation.R"))
create_GIN_annotation <- dget(paste0(func_dir, "create_GIN_annotation.R"))
create_QC_annotation <- dget(paste0(func_dir, "create_QC_annotation.R"))


#########################################################################################
### 1. Load data ###
#########################################################################################

# load samples:
for (j in 1:length(sample_ids)) {
  seurat_dir <- paste0(in_path, "seurat_", sample_ids[j], 
    "/Output/Rdata")
  in_dir <- paste0(in_path, "seurat_", sample_ids[j], 
    "/Output/InferCNV/supervised_clustering/subpop_mode/")
  plot_dir <- paste0(in_path, "seurat_", sample_ids[j], 
    "/Output/Plots/")

  # load seurat object:
  seurat_10X <- readRDS(paste0(seurat_dir, "/03_seurat_object_processed.RData"))

  # prepare infercnv metadata and write to files
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, "PC_A_res.1")
  names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
  seurat_10X <- infercnv_metadata$seurat
  
  normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
    unique(Idents(seurat_10X)), value=T, invert=T)
  print(paste0("Normals are: ", normals))
  
  
  ################################################################################
  ### 2. Create heatmap and annotations ###
  ################################################################################
  
  # remove CAFs from heatmap df and metadata:
  cells_to_remove <- 
  infercnv_metadata$metadata$cell_ids[grep("CAF", 
    infercnv_metadata$metadata$cell_type)]
  infercnv_metadata$metadata <- 
    infercnv_metadata$metadata[!(infercnv_metadata$metadata$cell_ids %in% 
    cells_to_remove),]
  print(dim(infercnv_metadata$metadata))
  
  infercnv_output_filename <- list.files(paste0(in_dir), 
    pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt", 
    full.names = T)
  
  infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))
  heatmap_df <- infercnv_output[!(rownames(infercnv_output) %in% cells_to_remove),]
  print(dim(heatmap_df))
  # create group annotation df
  group_annotation <- create_group_annotation(heatmap_df, infercnv_metadata$metadata)
  # ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
  m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
  heatmap_df <- heatmap_df[m,]
  chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)
  
  # determine normal cell ids:
  normal_cells <- infercnv_metadata$metadata$cell_ids[
    infercnv_metadata$metadata$cell_type %in% normal_clusters[[j]]
  ]
  # add GIN annotation:
  GIN_annotation <- create_GIN_annotation(heatmap_df, normal_cells)
  # save GIN scores:
  if (!file.exists(paste0(in_dir, "GIN_scores.txt"))) {
    write.table(GIN_annotation$GIN_levels, paste0(out_dir, 
      "GIN_scores.txt"), quote=F, row.names=F)
  }
  # plot GIN ranges:
  pdf(paste0(plot_dir, "GIN_range_plot.pdf"))
    print(GIN_annotation$GIN_range_plot)
  dev.off()

  ######
  # if length of GIN annotation > 4, cancer cells exist and the following 
  # correlation can be calculated:
  if (length(GIN_annotation) > 4 ) {
    # add correlation annotation:
    create_correlation_data <- function(df, metadata, GIN_levels, normal_cells) {
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
      sig_cors <- cor_df[cor_df$cor.p.value < 0.05,]
      sig_cors <- data.frame(cell_id=rownames(df), cor=sig_cors$cor.estimate)

      # record maximum levels of GIN in normals:
      max_normal_cor <- round(max(sig_cors$cor[sig_cors$cell_id %in% normal_cells]), 2)
      min_normal_cor <- round(min(sig_cors$cor[sig_cors$cell_id %in% normal_cells]), 2)
    
      if ( any(!(rownames(df) %in% normal_cells)) ) {
    
        max_other_cor <- 
          round(max(sig_cors$cor[!(rownames(df) %in% normal_cells)]), 2)
        min_other_cor <- 
          round(min(sig_cors$cor[!(rownames(df) %in% normal_cells)]), 2)
    
        sig_cors$type <- "other"
        sig_cors$type[sig_cors$cell_id %in% normal_cells] <- "normal"
    
        p <- ggplot(sig_cors, aes(x=type, y=cor, fill=type))
        p <- p + geom_violin()
        p <- p + scale_fill_manual(values = c("#FED976", "#9B59B6"))
        p <- p + annotate("text", x=1, y=max_normal_cor, 
          label=paste0("Max: ", max_normal_cor),
          vjust = -1)
        p <- p + annotate("text", x=1, y=min_normal_cor, 
          label=paste0("Min: ", min_normal_cor),
          vjust = 1.5)
        p <- p + annotate("text", x=2, y=max_other_cor, 
          label=paste0("Max: ", max_other_cor),
          vjust = -1)
        p <- p + annotate("text", x=2, y=min_other_cor, 
          label=paste0("Min: ", min_other_cor),
          vjust = 1.5)
        p <- p + theme_bw() + 
          theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text=element_text(size=12)
          )
        p <- p + ylab("correlation with top 5% cancer cells")
    
        result_list <- list(cor_df, p)
        names(result_list) <- c("correlation_df",  
          "correlation_range_plot")

      } else {

        sig_cors$type <- "normal"
        p <- ggplot(sig_cors, aes(x=type, y=cor, fill=type))
        p <- p + geom_violin()
        p <- p + scale_fill_manual(values = "#FED976")
        p <- p + annotate("text", x=1, y=max_normal_cor, 
          label=paste0("Max: ", max_normal_cor),
          vjust = -1)
        p <- p + annotate("text", x=1, y=min_normal_cor, 
          label=paste0("Min: ", min_normal_cor),
          vjust = 1.5)
        p <- p + theme_bw() + 
          theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text=element_text(size=12)
          )
        p <- p + ylab("correlation with top 5% cancer cells")
    
        result_list <- list(cor_df, p)
        names(result_list) <- c("correlation_df",  
          "correlation_range_plot")

      }
#      # generate quadrant plot for all cell values:
#      quad_df <- cor_df[cor_df$cor.estimate != "no_CNVs_recorded",]
#      quad_df$cor.estimate <- as.numeric(quad_df$cor.estimate)
#      quad_df$cor.p.value <- as.numeric(quad_df$cor.p.value)
#      quad_df$significance <- FALSE
#      quad_df$significance[quad_df$cor.p.value < 0.05] <- TRUE
#      colnames(quad_df) <- c("correlation", "p.value", "GIN", "significance")
#      quad_df <- subset(quad_df, select=c(correlation, GIN, significance))
#      quad_df$call <- "undetermined"
#      quad_df$call[quad_df$correlation > 0.4 & quad_df$GIN > 0.3] <- "malignant"
#      quad_df$call[quad_df$correlation < 0.4 & quad_df$GIN < 0.3] <- "non-malignant"
#      quad_df$call[!(rownames(quad_df) %in% epithelial_cells)] <- "stromal"
#      
#      p <- ggplot(quad_df, aes(x=correlation, y=GIN, col=call))
#      p <- p + geom_point()
#      p <- p + theme_minimal()
#      p <- p + coord_fixed()
#      p <- p + geom_vline(xintercept = 0.4) + geom_hline(yintercept = 0.3) 
      return(result_list)
    }
    correlation_annotation <- create_correlation_data(heatmap_df, 
      infercnv_metadata$metadata, GIN_annotation$GIN_levels, normal_cells)
  
    pdf(paste0(plot_dir, "cancer_correlation_ranges.pdf"))
      print(correlation_annotation$correlation_range_plot)
    dev.off()
  }
  ######
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
  
  old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
   new_scores <- c(-2, -1, 0, 1, 2, 3)
   for (s in 1:length(new_scores)) {
     plot_object[plot_object == old_scores[s]] <- new_scores[s]
   }
  
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  print("Generating final heatmap...")
  
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
   plot_object, name = paste0("hm"), 
   col = colorRamp2(c(-2, -1, 0, 1, 2, 3), 
     c("#00106B", "#9191CC", "white", "#DDB6B6", "#AB4848", "#930707"), 
     space = "sRGB"),
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
  GIN_annotation$GIN_annotation + QC_annotation$nUMI_annotation + 
  QC_annotation$nGene_annotation
  
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
  )
  
  # determine where starting co-ordinates for heatmap are based upon longest cluster name
  # (0.00604 units per character):
  longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
  x_coord <- longest_cluster_name*0.0037


  # if length of GIN annotation > 4, cancer cells exist and can be plotted:
  pdf(paste0(plot_dir, "final_infercnv_predictions_heatmap.pdf"), 
    height = 14.5, width = 18)
    grid.newpage()
    # plot Normal subtype:
    pushViewport(viewport(x = x_coord+0.003, y = 0.090,
                          width = 0.855+0.019, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[5]])
    popViewport()
    
    # plot Basal subtype:
    pushViewport(viewport(x = x_coord+0.008, y = 0.169,
                          width = 0.855+0.014, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[4]])
    popViewport()
    
    # plot Her2 subtype:
    pushViewport(viewport(x = x_coord+0.01, y = 0.251, 
                          width = 0.855+0.012, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[3]])
    popViewport()
    
    # plot LumB subtype:
    pushViewport(viewport(x = x_coord+0.007, y = 0.333, 
                          width = 0.855+0.015, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[2]])
    popViewport()
    
    # plot LumA subtype:
    pushViewport(viewport(x = x_coord+0.007, y = 0.411, 
                          width = 0.855+0.015, height = 0.08, just = c("left", "top")))
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
    pushViewport(viewport(x=x_coord + 0.895, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
    #grid.draw(lollipop)
    grid.text("GIN", rot=55)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.91, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
    #grid.draw(lollipop)
    grid.text("nUMI", rot=55)
    popViewport()
    pushViewport(viewport(x=x_coord + 0.925, y=0.395, width = 0.1, height = 0.1, just = "bottom"))
    #grid.draw(lollipop)
    grid.text("nGene", rot=55)
    popViewport()
    
  dev.off()

}

print(paste0("Group heatmap created, output in ", out_dir))

# convert pdf to png:
system(paste0("for p in ", out_dir, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
"new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
"convert -density 150 ", out_dir, "$f -quality 90 ", out_dir, "$new; done"))

