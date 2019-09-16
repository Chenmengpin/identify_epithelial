lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"

library(cluster, lib.loc = lib_loc)
library(ggplot2)
library(scales, lib.loc = lib_loc)
library(fpc, lib.loc = lib_loc)
library(dplyr)

project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_130819"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")

plot_dir <- paste0(out_dir, "/plots/")
system(paste0("mkdir -p ", plot_dir))
table_dir <- paste0(out_dir, "/tables/")
system(paste0("mkdir -p ", table_dir))

# load infercnv dataframe:
infercnv_data <- readRDS(paste0(seurat_path, 
  "/infercnv_measures_mean_of_scaled_squares.RData"))

# check total length of infercnv dataframe:
total_length <- sum(unlist(lapply(infercnv_data, nrow)))

for (i in 1:length(names(infercnv_data))) {

  print(i)

  infercnv_dir <- paste0(seurat_path, "seurat_", names(infercnv_data),
  	"/Output/InferCNV/")

  infercnv_measures <- infercnv_data[[i]]

  infercnv_output_filename <- list.files(paste0(seurat_path, "seurat_", 
      names(infercnv_data)[i], "/Output/InferCNV/supervised_clustering/subpop/"), 
      pattern = "infercnv.1[2,5]_denoised.observations.txt", 
      full.names = T)

  infercnv_output <- as.data.frame(t(read.table(infercnv_output_filename)))

  print(paste0("No rows in infercnv_output: ", nrow(infercnv_output)))
  print(paste0("No rows in corresponding infercnv_data: ", nrow(infercnv_measures)))

  #  prepare data for quad plots:
  infercnv_level_vs_correlation <- subset(infercnv_measures, select = c(infercnv_levels, significant_infercnv_correlation_0.05))
  
  # scale data and visualise:
  scaled_infercnv_level_vs_correlation <- scale(infercnv_level_vs_correlation) %>% as.data.frame()

  # run silhouette cluster analysis to determine clusters and thresholds:
  pamk_result <- pamk(scaled_infercnv_level_vs_correlation, krange=1:4)
  pamk_result$nc
  silhouette_result <- pam(scaled_infercnv_level_vs_correlation, 
                           pamk_result$nc)

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