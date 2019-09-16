#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R

library(ggplot2)
library(Seurat)
library(scales)
library(fpc)
library(dplyr)
library(cluster)

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
seurat_path <- paste0(project_dir, "results/seurat/brca_mini_atlas_130819/")
Rdata_dir <- paste0(seurat_path, "/Rdata/")
system(paste0("mkdir -p ", Rdata_dir))
seurat_filename <- "04_seurat_object_annotated.RData"
include_calls <- FALSE
include_HMM <- FALSE

sample_ids <- gsub("seurat_", "", list.files(seurat_path, pattern="CID"))

######
#sample_ids <- grep("CID4290A", sample_ids, invert=T, value=T)
#sample_ids <- c("CID3948", "CID3941", "CID4067", "CID4461", "CID4463", "CID4465", 
#  "CID4398", "CID3921", "CID3963", "CID45171", "CID44041", "CID4523", "CID3586", 
#  "CID44971")
######

#########################################################################################
### 0. Load data ###
#########################################################################################

# load integrated object to add scores to:
integrated_object <- readRDS(
  "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Jun2019/02_integration/output/CCA_CCA23Jun2019/Output/Rdata/03_seurat_CCA_aligned_processed.Rdata"
)

# 37382 epithelial cells in integrated object

# fetch cell ids of epithelial cells in integrated object:
integrated_epi <- names(Idents(integrated_object))[
  grep("pithelial", integrated_object@meta.data$garnett_call_ext_major)
]
# fetch corresponding CNV measures:
epi_CNV_measures <- HMM_CNV_measures_df[
  rownames(HMM_CNV_measures_df) %in% integrated_epi,
]

# only select samples specified in sample_ids:
custom_integrated_epi <- integrated_epi[gsub("\\_.*", "", integrated_epi) %in% sample_ids]
# check if any epithelial cell values are in integrated_epi but not epi_CNV_measures:
custom_integrated_epi[which(!(custom_integrated_epi %in% rownames(epi_CNV_measures)))]

# load samples:
for (i in 1:length(sample_ids)) {
  print(i)
  if (i==1) {
    infercnv_filename <- grep(
      "Target", list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
        "/Output/InferCNV/supervised_clustering/subpop"),
        pattern = "infercnv.15_denoised.observations.txt",
        recursive=T, full.names=T), invert=T, value=T
    )
    if (length(infercnv_filename) < 1) {
      infercnv_filename <- list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
      "/Output/InferCNV/supervised_clustering/sample"),
      pattern = "infercnv.12_denoised.observations.txt",
      recursive=T, full.names=T)
    }
    print(infercnv_filename)
    infercnv_filenames <- list(infercnv_filename)
  	HMM_filenames <- list(
  	  list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
      "/Output/InferCNV/supervised_clustering/subpop"),
      pattern = "infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.repr_intensities.observations.txt",
      recursive=T, full.names=T)
    )
  } else {
    infercnv_filename <- grep(
      "Target", list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
        "/Output/InferCNV/supervised_clustering/subpop"),
        pattern = "infercnv.15_denoised.observations.txt",
        recursive=T, full.names=T), invert=T, value=T
    )
    if (length(infercnv_filename) < 1) {
      infercnv_filename <- list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
      "/Output/InferCNV/supervised_clustering/sample"),
      pattern = "infercnv.12_denoised.observations.txt",
      recursive=T, full.names=T)
    }

    print(infercnv_filename)
    infercnv_filenames <- c(infercnv_filenames, infercnv_filename)
    HMM_filenames <- c(HMM_filenames, 
      list(list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
      "/Output/InferCNV/supervised_clustering/subpop"),
      pattern = "infercnv.14_HMM_predHMMi6.*.hmm_mode-subclusters.repr_intensities.observations.txt",
      recursive=T, full.names=T)))
  }
}
names(infercnv_filenames) <- sample_ids
names(HMM_filenames) <- sample_ids

# load inferCNV output, aggregating data sets that have been split:
infercnv_output <- lapply(infercnv_filenames, function(x) {
  return(as.data.frame(t(read.table(x))))
})

#save.image(paste0(Rdata_dir, "/all_infercnv_and_integrated_data_loaded.RData"))

# find CNV levels and correlation with 5% cancer:
infercnv_value_type <- "sum_of_abs_value"
#infercnv_value_type <- "mean_of_squares"
z=1
infercnv_measures <- lapply(infercnv_output, function(x) {
  print(sample_ids[z])
  # sum total CNV levels for each cell and normalise to gene no:
  # infercnv_levels <- apply(x, 1, function(y) {
  #   y[is.na(y)] <- 0
  #   return(sum(y))
  # })
  # # normalise to number of genes in plot:
  # infercnv_levels <- infercnv_levels/ncol(x)
  infercnv_levels <- apply(x, 1, function(y) {
    y[is.na(y)] <- 0
    return(round(sum(abs(y))/ncol(x), 6))
    #return(mean(y^2))
  })
  
  # normalise to number of genes in plot:
  #infercnv_levels <- infercnv_levels/ncol(x)
  # calculate correlation with top 5% cancer:
  epithelial_cells <- rownames(x)
  epithelial_infercnv <- infercnv_levels[names(infercnv_levels) %in% epithelial_cells]
  epithelial_infercnv <- epithelial_infercnv[order(epithelial_infercnv, decreasing=T)]
  top_infercnv <- head(epithelial_infercnv, length(epithelial_infercnv)*0.05)
  # find average genome-wide CNV predictions across genome:
  top_infercnv_average <- apply(x[names(top_infercnv),], 2, mean)
  # find correlations of each cell's CNVs with top_infercnv_average:
  infercnv_correlations <- apply(x, 1, function(y) {
    if (length(unique(as.numeric(y))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(y), top_infercnv_average)
      # create significant correlation vector only including p < 0.05:
      sig_cor_0.05 <- cor$estimate
      sig_cor_0.05[sig_cor_0.05 < 0.05] <- 0
      cor_result <- data.frame(round(cor$estimate, 6), round(cor$p.value, 6), round(sig_cor_0.05, 6))
    }
    return(cor_result)
  })
  cor_df <- do.call("rbind", infercnv_correlations)
  combined_df <- cbind(infercnv_levels, cor_df)
  if (include_calls) {
    # create columns of normal calls for different CNV and correlation cutoffs:
    infercnv_cutoffs <- c(0.1, 0.2, 0.3)
    for (l in 1:3) {
      cutoff_col_1 <- rep("non_normal", nrow(combined_df))
      cutoff_col_1[combined_df$infercnv_levels < infercnv_cutoffs[l]] <- "normal"
      cutoff_col_2 <- rep("non_normal", nrow(combined_df))
      cutoff_col_2[
        combined_df$infercnv_levels < infercnv_cutoffs[l] & combined_df$sig_cor_0.05 < 0.3
      ] <- "normal"
  
      cutoff_df <- data.frame(cutoff_col_1, cutoff_col_2)
      
      combined_df <- cbind(combined_df, cutoff_df)
  
      colnames(combined_df)[(ncol(combined_df)-1):(ncol(combined_df))] <- c(
        paste0("normal_call_", infercnv_cutoffs[l], "_infercnv_cutoff"), 
        paste0("normal_call_", infercnv_cutoffs[l], "_infercnv_0.3_correlation_cutoffs")
      )
    }
  } 
  z <<- z+1
  return(combined_df)
})
names(infercnv_measures) <- sample_ids
infercnv_measures_df <- do.call("rbind", infercnv_measures)
rownames(infercnv_measures_df) <- gsub("^.*\\.", "", rownames(infercnv_measures_df))
colnames(infercnv_measures_df)[1:4] <- c("infercnv_levels", "infercnv_correlation_top_0.05_cancer", 
  "infercnv_correlation_p_value", "significant_infercnv_correlation_0.05")

# save infercnv measures as text file and RDS:
write.table(infercnv_measures_df, paste0(seurat_path, "/infercnv_measures_", infercnv_value_type, "_df.txt"))
saveRDS(infercnv_measures, paste0(seurat_path, "infercnv_measures_", infercnv_value_type, ".RData"))

# add new columns to integrated_object metadata:
integrated_object@meta.data$infercnv_levels <- rep(NA, nrow(integrated_object@meta.data))
integrated_object@meta.data$infercnv_correlation_top_0.05_cancer <- 
  rep(NA, nrow(integrated_object@meta.data))
integrated_object@meta.data$infercnv_correlation_p_value <- 
  rep(NA, nrow(integrated_object@meta.data))
integrated_object@meta.data$significant_infercnv_correlation_0.05 <- 
  rep(NA, nrow(integrated_object@meta.data))

# find locations of cells in infercnv_measures_df in integrated object metadata:
m <- match(rownames(infercnv_measures_df), rownames(integrated_object@meta.data))

# add infercnv values to seurat object metadata:
integrated_object@meta.data$infercnv_levels[m] <- infercnv_measures_df$infercnv_levels

integrated_object@meta.data$infercnv_correlation_top_0.05_cancer[m] <- 
  infercnv_measures_df$infercnv_correlation_top_0.05_cancer

integrated_object@meta.data$infercnv_correlation_p_value[m] <- 
  infercnv_measures_df$infercnv_correlation_p_value

integrated_object@meta.data$significant_infercnv_correlation_0.05[m] <- 
  infercnv_measures_df$significant_infercnv_correlation_0.05

# check infercnv levels present for all epithelial cells:
epi_cells <- colnames(integrated_object)[grep("pithelial", 
  integrated_object@meta.data$garnett_call_ext_major)]
m <- match(epi_cells, colnames(integrated_object))
all_infercnv_levels <- integrated_object@meta.data$infercnv_levels[m]
names(all_infercnv_levels) <- colnames(integrated_object)[m]

any(is.na(integrated_object@meta.data$infercnv_levels[m]))
which(is.na(integrated_object@meta.data$infercnv_levels[m]))

# plot histograms and density plots of infercnv levels for each sample:
custom_sample_ids <- c("CID4463", "CID4471")
for (i in 1:length(custom_sample_ids)) {
  print(custom_sample_ids[i])
  # define output directory:
  out_dir <- paste0(seurat_path, "seurat_", custom_sample_ids[i], "/Output/Plots/")
  system(paste0("mkdir -p ", out_dir))
  
  # fetch infercnv measures df for sample:
  measure_df <- eval(parse(text=paste0("infercnv_measures$", custom_sample_ids[i])))
  
  # create density plot of infercnv values:
  density_plot <- density(measure_df$infercnv_levels, bw="SJ")
  
  #  prepare data for quad plots:
  infercnv_level_vs_correlation <- 
    subset(measure_df, select = c(infercnv_levels, sig_cor_0.05))
  
  # scale data and visualise:
  scaled_infercnv_level_vs_correlation <- scale(infercnv_level_vs_correlation) %>% as.data.frame()
  
  # scaled_scatter <- ggplot(infercnv_level_vs_correlation,
  #   aes(x=infercnv_levels, y=sig_cor_0.05)) +
  #   geom_point() +
  #   theme_minimal()
  # scaled_scatter
  
  # run silhouette cluster analysis to determine clusters and thresholds:
  pamk_result <- pamk(scaled_infercnv_level_vs_correlation, krange=2:4)
  pamk_result$nc
  silhouette_result <- pam(scaled_infercnv_level_vs_correlation, 
                           pamk_result$nc)
  
  # plot silhouette values:
  sil_values <- as.data.frame(silhouette_result$silinfo$widths)
  #barplot(sil_values$sil_width)
  # sil_values$cell_ids <- rownames(sil_values)
  # sil_values <- sil_values[order(sil_values$cluster),]
  # p <- ggplot(sil_values, aes(x=cell_ids, y=sil_width))
  # p <- p + geom_bar(stat="identity")
  
  infercnv_level_vs_correlation <- merge(
    measure_df,
    data.frame(row.names=names(silhouette_result$clustering),
               cluster=silhouette_result$clustering),
    by="row.names"
  )
  
  # determine normal cluster:
  cluster_split <- split(infercnv_level_vs_correlation, 
                         infercnv_level_vs_correlation$cluster)
  cluster_max_levels <- lapply(cluster_split, function(x) max(x$sig_cor_0.05))
  max_vals <- do.call("c", cluster_max_levels)
  normal_cluster <- which(max_vals == min(max_vals))
  cancer_cluster <- which(max_vals == max(max_vals))
  x_upper_lim_normal <- max(cluster_split[[normal_cluster]]$infercnv_levels)
  x_lower_lim_cancer <- min(cluster_split[[cancer_cluster]]$infercnv_levels)
  x_int <- x_lower_lim_cancer + (x_upper_lim_normal - x_lower_lim_cancer)/2
  y_int <- max(cluster_split[[normal_cluster]]$sig_cor_0.05)
  
  if (pamk_result$nc == 3) {
    
    # create quad plot:
    p <- ggplot(infercnv_level_vs_correlation, 
                aes(x=infercnv_levels, y=sig_cor_0.05, color=as.factor(cluster)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("black", "#b2182b", "#74add1"), 
                                labels=c("Unassigned", "Cancer", "Normal"))
    p <- p + xlab("Infercnv level")
    p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    p
    quad_plot <- p
  } else if (pamk_result$nc == 2) {
    
    # create quad plot:
    p <- ggplot(infercnv_level_vs_correlation, 
                aes(x=infercnv_levels, y=sig_cor_0.05, color=as.factor(cluster)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("#74add1", "#b2182b"), 
                                labels=c("Normal", "Cancer"))
    p <- p + xlab("Infercnv level")
    p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
    p <- p + theme(legend.title = element_blank())
    p <- p + geom_vline(xintercept = x_int)
    p <- p + geom_hline(yintercept = y_int)
    p
    quad_plot <- p
  } else if (pamk_result$nc == 1) {
    # create quad plot:
    p <- ggplot(infercnv_level_vs_correlation, 
                aes(x=infercnv_levels, y=sig_cor_0.05, color=as.factor(cluster)))
    p <- p + geom_point()
    p <- p + scale_color_manual(values=c("#b2182b"), 
                                labels=c("Cancer"))
    p <- p + xlab("Infercnv level")
    p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
    p <- p + theme(legend.title = element_blank())
    p
    quad_plot <- p
  }
  
  png(paste0(out_dir, "infercnv_level_", infercnv_value_type, "_distributions_.png"), width = 860, height = 400)
    par(mfrow=c(1,2))
    hist(measure_df$infercnv_levels, main = NULL, xlab = "InferCNV levels")
    plot(density_plot, main=NA, xlab = "Infercnv value")
  dev.off()
  
  #png(paste0(out_dir, "quad_plot_cluster_silhouette_scores2.png"), width = 430, height = 200)
  #  barplot(sil_values$sil_width)
  #dev.off()
  
  png(paste0(out_dir, "normal_call_quad_plot_", infercnv_value_type, ".png"), width = 430, height = 200)
    print(quad_plot)
  dev.off()

}
  
# create quad plot:
p <- ggplot(infercnv_level_vs_correlation, 
  aes(x=infercnv_levels, y=sig_cor_0.05, color=as.factor(cluster)))
p <- p + geom_point()
p <- p + scale_color_manual(values=c("black", "#b2182b", "#74add1"), 
  labels=c("Unassigned", "Cancer", "Normal"))
p <- p + xlab("Infercnv level")
p <- p + ylab("Corr. with top 5% cancer (p<0.05)")
p <- p + theme(legend.title = element_blank())
p <- p + geom_vline(xintercept = x_int)
p <- p + geom_hline(yintercept = y_int)
p


grid.newpage()
pushViewport(viewport(x = 0.01, y = 0.99,
 width = 0.25, height = 0.25, just = c("left", "top")))
hist(eval(parse(text=paste0("infercnv_measures$", sample_id, 
                            "$infercnv_levels"))))
popViewport()
pushViewport(viewport(x = 0.01, y = 0.99,
                      width = 0.25, height = 0.25, just = c("left", "top")))
plot(density_plot, main = "InferCNV value density", xlab = "Infercnv value")
popViewport()
  
  density_plot <- density(
    eval(parse(text=paste0("infercnv_measures$", sample_id, 
    "$infercnv_levels"))), bw="SJ")
   p

hist(wt, main="Histogram of wt")
boxplot(wt, main="Boxplot of wt")

# # create scatter plot of infercnv levels for each sample:
# for (s in list(infercnv_measures$CID4463, infercnv_measures$CID4471)) {
#   infercnv_level_vs_correlation <- 
#     subset(s, select = c(infercnv_levels, sig_cor_0.05))
#   scatter_plot <- ggplot(infercnv_level_vs_correlation, 
#     aes(x=infercnv_levels, y=sig_cor_0.05)) +
#     geom_point() +
#     theme_minimal()
#   scatter_plot
#   
#   # run k-means clustering analysis on quad plot data:
#   scaled_infercnv_level_vs_correlation <- 
#     scale(infercnv_level_vs_correlation) %>% as.data.frame()
#   
#   for (i in 1:3) {
#     kmeans_clusters <- kmeans(scaled_infercnv_level_vs_correlation, centers = i)
#     
#     clustered_infercnv_measurements <- merge(
#       s, 
#       data.frame(row.names=names(kmeans_clusters$cluster), 
#                  cluster=kmeans_clusters$cluster),
#       by="row.names"
#     )
#     
#     grouped_scatter_plot <- ggplot(clustered_infercnv_measurements, 
#                                    aes(x=infercnv_levels, y=sig_cor_0.05, color=cluster)) +
#       geom_point() +
#       theme_minimal()
#     print(grouped_scatter_plot)
#   }
# }

CID4463_dir <- paste0(seurat_path, 
                      "/seurat_CID4463/Output/InferCNV/Plots/")
png(paste0(CID4463_dir, "infercnv_value_plots.png"))

dev.off()

#HMM_filenames <- HMM_filenames[c(1:7, 9:13, 16:23)]

if (include_HMM) {
  # load inferCNV HMM output, aggregating data sets that have been split:
  HMM_output <- lapply(HMM_filenames, function(x) {
    print(x)
    if (length(x) > 1) {
      output_dfs <- lapply(x, function(y) as.data.frame(t(read.table(y))))
      # remove columns missing from one or more dfs:
      for (d in 1:length(output_dfs)) {
        if (d==1) {
          keep_colnames <- colnames(output_dfs[[d]])
        } else {
          missing_cols <- c(
            keep_colnames[!(keep_colnames %in% colnames(output_dfs[[d]]))],
            colnames(output_dfs[[d]])[
              !(colnames(output_dfs[[d]]) %in% keep_colnames)
            ]
          )
          keep_colnames <- keep_colnames[!(keep_colnames %in% missing_cols)]
        }
      }
      HMM_output <- lapply(output_dfs, function(y) subset(y, select=keep_colnames))
      # join dataframes and remove duplicates
      HMM_output <- do.call("rbind", HMM_output)
      HMM_output <- HMM_output[!duplicated(rownames(HMM_output)),]
    } else {
      HMM_output <- as.data.frame(t(read.table(x)))
    }
    return(HMM_output)
  })
  
  # find CNV levels and correlation with 5% cancer:
  HMM_CNV_measures <- lapply(HMM_output, function(x) {
    old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
    new_scores <- c(2, 1, 0, 1, 2, 3)
    for (s in 1:length(new_scores)) {
      x[x == old_scores[s]] <- new_scores[s]
    }
    # sum total CNV levels for each cell and normalise to gene no:
    HMM_CNV_levels <- apply(x, 1, function(y) {
      y[is.na(y)] <- 0
      return(sum(y))
    })
    # normalise to number of genes in plot:
    HMM_CNV_levels <- round(abs(HMM_CNV_levels)/ncol(x), 3)
    # calculate correlation with top 5% cancer:
    epithelial_cells <- rownames(x)
    epithelial_CNV <- HMM_CNV_levels[names(HMM_CNV_levels) %in% epithelial_cells]
    epithelial_CNV <- epithelial_CNV[order(epithelial_CNV, decreasing=T)]
    top_CNV <- head(epithelial_CNV, length(epithelial_CNV)*0.05)
    # find average genome-wide CNV predictions across genome:
    top_CNV_average <- apply(x[names(top_CNV),], 2, mean)
    # find correlations of each cell's CNVs with top_CNV_CNV_average:
    CNV_correlations <- apply(x, 1, function(y) {
      if (length(unique(as.numeric(y))) == 1) {
        cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
          cor.p.value="no_CNVs_recorded")
      } else {
        cor <- cor.test(as.numeric(y), top_CNV_average)
        # create significant correlation vector only including p < 0.05:
        sig_cor_0.05 <- cor$estimate
        sig_cor_0.05[sig_cor_0.05 < 0.05] <- 0
        cor_result <- data.frame(cor$estimate, cor$p.value, sig_cor_0.05)
      }
      return(cor_result)
    })
    cor_df <- do.call("rbind", CNV_correlations)
    combined_df <- cbind(HMM_CNV_levels, cor_df)
    # create columns of normal calls for different CNV and correlation cutoffs:
    CNV_cutoffs <- c(0.1, 0.2, 0.3)
    for (l in 1:3) {
      cutoff_col_1 <- rep("non_normal", nrow(combined_df))
      cutoff_col_1[combined_df$HMM_CNV_levels < CNV_cutoffs[l]] <- "normal"
      cutoff_col_2 <- rep("non_normal", nrow(combined_df))
      cutoff_col_2[
        combined_df$HMM_CNV_levels < CNV_cutoffs[l] & combined_df$sig_cor_0.05 < 0.3
      ] <- "normal"
  
      cutoff_df <- data.frame(cutoff_col_1, cutoff_col_2)
      
      combined_df <- cbind(combined_df, cutoff_df)
  
      colnames(combined_df)[(ncol(combined_df)-1):(ncol(combined_df))] <- c(
        paste0("normal_call_", CNV_cutoffs[l], "_HMM_CNV_cutoff"), 
        paste0("normal_call_", CNV_cutoffs[l], "_HMM_CNV_0.3_correlation_cutoffs")
      )
    }
    
    return(combined_df)
  })
  HMM_CNV_measures_df <- do.call("rbind", HMM_CNV_measures)
  rownames(HMM_CNV_measures_df) <- gsub("^.*\\.", "", rownames(HMM_CNV_measures_df))
  colnames(HMM_CNV_measures_df) <- c("HMM_CNV_levels", "HMM_CNV_correlation_0.05_cancer", 
    "HMM_CNV_correlation_p_value", "significant_HMM_CNV_correlation_0.05")
  HMM_minimal <- subset(HMM_CNV_measures_df, select=c(HMM_CNV_levels, 
    HMM_CNV_correlation_0.05_cancer, HMM_CNV_correlation_p_value, 
    significant_HMM_CNV_correlation_0.05))
  infercnv_measures_df <- merge(infercnv_measures_df, HMM_minimal, by="row.names")
}

#write.table(infercnv_measures_df, paste0(seurat_path, "/infercnv_measures_", infercnv_value_type, "_df.txt"))
#saveRDS(infercnv_measures, paste0(seurat_path, "infercnv_measures_", infercnv_value_type, ".RData"))



#for (i in 1:length(sample_ids)) {
#  if (i==1) {
#  	GIN_filenames <- list(
#  	  list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
#      "/Output/InferCNV/supervised_clustering/subpop"),
#      pattern = "GIN",
#      recursive=T, full.names=T)
#    )
#    correlation_filenames <- list(
#      list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
#      "/Output/InferCNV/supervised_clustering/subpop"),
#      pattern = "correlation",
#      recursive=T, full.names=T)
#    )
#  } else {
#    GIN_filenames[[i]] <- list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
#      "/Output/InferCNV/supervised_clustering/subpop"),
#      pattern = "GIN",
#      recursive=T, full.names=T)
#    correlation_filenames[[i]] <- list.files(paste0(seurat_path, "seurat_", sample_ids[i], 
#      "/Output/InferCNV/supervised_clustering/subpop"),
#      pattern = "correlation",
#      recursive=T, full.names=T)
#  }
#}
#
#names(seurat_list) <- sample_ids
#
## remove infercnv done on whole sample if exists:
#      output_files <- as.list(output_files[
#        grep("[0-9]/infercnv|Aggregate_group/infercnv", output_files)
#      ])
#      # load dataframes:
#      output_dfs <- lapply(output_files, function(x) as.data.frame(t(read.table(x))))#