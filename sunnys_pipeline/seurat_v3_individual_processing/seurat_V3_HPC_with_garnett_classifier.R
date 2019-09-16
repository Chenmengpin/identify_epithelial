#### SEURAT V3 HPC R SCRIPT WITH DROPLET UTILS FILTERING
#     written by; Sunny Z Wu
#     last modified; 20190220
#
#     # New commands in seurat v3
#     https://satijalab.org/seurat/essential_commands.html
#
### inputs (in order);
      # 01 PROJECT/SAMPLE NAME
      # 02 SPECIES (human/mouse)
      # 03 WORKING DIRECTORY (if submitting job from directory use $(pwd))
      # 04 PATH TO CELLRANGER RAW MATRIX
      # 05 GENES TO PLOT FILE
            # - genes.csv
                # panel,CAFs,immune
                # ACTB,PDGFRB,ITGAM
                # EPCAM,PDGFRA,LY6G6D
                # KRT18,COL1A1,LY6C1
                # KRT5,COL1A2,CD33
                # ESR1,ACTA2,FUT4
                # PGR,MCAM,CD15
                # ERBB2,PDPN,CD68
                # COL1A1,FAP,ADGRE1
                # PDGFRB,CD34,CD274
                # PDGFRA,CXCL12,PDCD1LG2
                # ACTA2,THY1,PDCD1
                # CD34,S100A4,LGALS9
                # PECAM1,,CD80
                # PTPRC,,CD86
                # CD3D,,CD276
                # CD8A,,TNFRSF18
                # FOXP3,,TNFSF18
                # CD19,,
                # JCHAIN,,
                # CD68,,
      # 06 SEURAT PARAMS FILE
          # - params.csv
          # Input params file should be in the following format of a csv file;
                # SEURAT_PROCESSING_PARAMS_FILE,OPTION
                # #PERFORM_EMPTYDROPS_FILTERING?,--
                # Emptydrops_filtering,T
                # Emptydrops_filtering_FDR,0.1
                # #IF_FALSE_PERFORM_MANUAL_FILTERING,--
                # Gene_high_cutoff,Inf
                # Gene_low_cutoff,250
                # UMI_high_cutoff,Inf
                # UMI_low_cutoff,500
                # #MITOCHONDRIAL_CUTOFF,--
                # Mitochondrial_UMI_cutoff,0.2
                # #REGRESS_OUT_VARIATION,--
                # Regress_out_variation_on_nUMIs_and_nGenes,T
                # #PERFORM_JACKSTRAW_PC_CUTOFF,--
                # Jackstraw_PC_cutoff,F
                # Number_of_PCs_to_compute,100
                # jackstraw_pval_cutoff_1,0.0001
                # jackstraw_pval_cutoff_2,0.01
                # jackstraw_pval_cutoff_3,0.05
                # #IF_FALSE_PERFORM_MANUAL_PC_CUTOFFS,--
                # PC_CUTOFF_1,20
                # PC_CUTOFF_2,30
                # PC_CUTOFF_3,50
                # #CLUSTERING_RESOLUTIONS,--
                # RES_1,0.4
                # RES_2,0.8
                # RES_3,1.2
                # #DIFFERENTIAL_GENE_EXPRESSION,--
                # Perform_differential_gene_expression,F
                # minimum_fraction_detected_min.pct,0.5
                # minimum_fraction_diff_min.diff.pct,0.05
                # minimum_threshold_diff_min.thresh.use,0.75
                # #EXPORT_MATRICES,--
                # Export_raw_UMI_matrix,T
                # Export_raw_normalised_matrix,T
### outputs; 
        # - All raw and processed seurat objects at Output/Rdata
        # - Output/Figures
        # - Emptydrops filtering metrics at Output/EmptyDrops
        # - All DGE and metrics at Output/Gene_expression_and_stats
        # - Matrices and metadata at Output/Matrices_and_metadata
        # - Rmarkdown html summmary file at params at Output/Analysis_params_and_Rmd
#
### REQUIRES R v.3.5.0
### QSUB ARGUMENTS
#     qsub 
#     -cwd 
#     -pe smp 32 
#     -l h_vmem=200G 
#     -P TumourProgression 
#     -b y 
#     -j y 
#     -V 
#     -N s_sampleID
#     "R CMD BATCH 
#     --no-save 
#     '--args 
#     sampleID
#     human
#     /working/directory/path/
#     /path_to_raw_matrix/GRCh38
#     gene_set_file.csv
#     seurat_params_file.csv'
#     /path/to/this/script.R" 
#
### QSUB ONE LINER:
# qsub -cwd -pe smp 16 -l h_vmem=300G -b y -j y -V -P TumourProgression "${R} CMD BATCH --no-save '--args CID44041 human $(pwd) /share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4404/count_4404_primary_GRCh38/outs/raw_gene_bc_matrices/GRCh38/ ./seurat_gene_input_file.csv ./seurat_params_file.csv' ./seurat_v3_HPC.R"
#
# LOAD PACKAGES -----------------------------------------------------------

library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)
library(tidyr)
library(eply)
library(ggplot2)

# COMMAND LINE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- date()
temp_args <-
  commandArgs(trailingOnly = T)

# 01 PROJECT/SAMPLE NAME
temp_project_name <- temp_args[1]
# 02 SPECIES (human/mouse)
temp_species_type <- temp_args[2]
# 03 WORKING DIRECTORY (if submitting job from directory use $(pwd))
temp_wd <- temp_args[3]
setwd(temp_wd)

# 04 PATH TO CELLRANGER RAW MATRIX
temp_raw_matrix <- temp_args[4]
# 05 GENES TO PLOT FILE
temp_gene_plot_file <- temp_args[5]
# 06 SEURAT PARAMS FILE
temp_params_file <- temp_args[6]

# READ PARAMS -----------------------------------------------------

temp_params <- 
  read.csv(temp_params_file, 
           row.names = "SEURAT_PROCESSING_PARAMS_FILE")


# Input as matrix.txt file (instead of 10X output)
temp_matrix_txt_file <- 
  as.logical(temp_params["Input_as_matrix_txt_file",]) 

# Use emptydrops function from DropletUtils for filtering?
temp_use_droplet_utils_filtering <-
  as.logical(temp_params["Emptydrops_filtering",]) 

if(temp_use_droplet_utils_filtering){
  library(DropletUtils)
}

# Filtering thresholds
# FDR rate for dropletUtils filtering (recommended is 0.01)
temp_FDR_threshold <- 
  as.numeric(as.character(temp_params["Emptydrops_filtering_FDR",]))

#nGene filtering params
temp_nGene.high_threshold <- 
  as.numeric(as.character((temp_params["Gene_high_cutoff",]))) 
temp_nGene.low_threshold <- 
  as.numeric(as.character((temp_params["Gene_low_cutoff",]))) 

#nUMI filtering params
temp_nUMI.high_threshold <- 
  as.numeric(as.character((temp_params["UMI_high_cutoff",])))  
temp_nUMI.low_threshold <- 
  as.numeric(as.character((temp_params["UMI_low_cutoff",])))  

#percentage of mitochondrial UMIs
temp_percent.mito.high_threshold <- 
  as.numeric(as.character((temp_params["Mitochondrial_UMI_cutoff",]))) 

# Variables to regress out
temp_vars.to.regress <-
  as.logical(temp_params["Regress_out_variation_on_nUMIs_and_nGenes",]) 

temp_vars.to.regress <- {if(temp_vars.to.regress) 
  (c("nFeature_RNA", "nCount_RNA")) else NULL
}

# PCA analysis and t-SNE resolution (reduce for less Jackstraw processing time)
temp_jackstraw_cutoff <- 
  as.logical(temp_params["Jackstraw_PC_cutoff",]) 

temp_PCs_to_compute <-
  as.numeric(as.character((temp_params["Number_of_PCs_to_compute",]))) 
temp_PC_A <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_1",])))  
temp_PC_B <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_2",])))  
temp_PC_C <- 
  as.numeric(as.character((temp_params["jackstraw_pval_cutoff_3",])))  

# MANUAL PC CUTOFFS
if(temp_jackstraw_cutoff == F) {
  temp_PC_A <- 
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_1",])))) 
  temp_PC_B <-
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_2",])))) 
  temp_PC_C <- 
    (1:as.numeric(as.character((temp_params["PC_CUTOFF_3",])))) 
}

# CLUSTERING RESOLUTIONS
temp_res_1 <- 
  as.numeric(as.character((temp_params["RES_1",])))  
temp_res_2 <- 
  as.numeric(as.character((temp_params["RES_2",])))  
temp_res_3 <- 
  as.numeric(as.character((temp_params["RES_3",])))  

#Differential gene expression params
# Do differential expression?
temp_do_differential_expression <- 
  as.logical(temp_params["Perform_differential_gene_expression",]) 
temp_min.pct <- 
  as.numeric(as.character((temp_params["minimum_fraction_detected_min.pct",]))) 
temp_min.diff.pct <- 
  as.numeric(as.character((temp_params["minimum_fraction_diff_min.diff.pct",])))  
temp_thresh.use <- 
  as.numeric(as.character((temp_params["minimum_threshold_diff_min.thresh.use",])))  

# Run Garnett cell type identification
temp_run_garnett <- 
  as.logical(temp_params["Run_Garnett",]) 

if(temp_run_garnett) {
  library(monocle)
  library(garnett)
  library(org.Hs.eg.db)
}


#exporting matrix and associated metadata (T or F)
temp_export_normalised_matrix <- 
  as.logical(temp_params["Export_raw_UMI_matrix",]) 
temp_export_raw_UMI_matrix <- 
  as.logical(temp_params["Export_raw_normalised_matrix",]) 

# SET UP AND FUNCTIONS ------------------------------------------------------------------

#Error processing
options(error=expression(NULL))

#Sub-directory Outputs
dir.create("Output")
dir.create("Output/Gene_expression_and_stats")
dir.create("Output/Figures")
dir.create("Output/Figures/Metrics/")
dir.create("Output/Figures/Gene_markers/")
dir.create("Output/Figures/Gene_markers/TSNE/")
dir.create("Output/Figures/Gene_markers/UMAP/")
dir.create("Output/Figures/TSNE/")
dir.create("Output/Figures/UMAP/")
dir.create("Output/Figures/Heatmaps/")
dir.create("Output/Figures/Cell_type_ID_Garnett/")
dir.create("Output/Figures/Cell_type_ID_Garnett/01_raw_garnett_calls/")
dir.create("Output/Figures/Cell_type_ID_Garnett/02_cluster_call_extended/")
dir.create("Output/Rdata")
dir.create("Output/Matrices_and_metadata")
dir.create("Output/EmptyDrops")
dir.create("Output/Analysis_params_and_Rmd")

# PNG function
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 14, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_small <-
  function(x) {
    png(
      file = (x), 
      width = 5, 
      height = 5, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_big <-
  function(x) {
    png(
      file = (x), 
      width = 22, 
      height = 12, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_longer <-
  function(x) {
    png(
      file = (x), 
      width = 14, 
      height = 12, 
      res = 300, 
      units = 'in'
    )
  }

# first letter up function for mouse gene conversion
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# function to list directories
twee <- function(path = getwd(), level = Inf) {
  
  fad <-
    list.files(path = path, recursive = TRUE,no.. = TRUE, include.dirs = TRUE)
  
  fad_split_up <- strsplit(fad, "/")
  
  too_deep <- lapply(fad_split_up, length) > level
  fad_split_up[too_deep] <- NULL
  
  jfun <- function(x) {
    n <- length(x)
    if(n > 1)
      x[n - 1] <- "|__"
    if(n > 2)
      x[1:(n - 2)] <- "   "
    x <- if(n == 1) c("-- ", x) else c("   ", x)
    x
  }
  fad_subbed_out <- lapply(fad_split_up, jfun)
  
  cat(unlist(lapply(fad_subbed_out, paste, collapse = "")), sep = "\n")
}

# READ10X MATRIX ----------------------------------------------------

if(temp_matrix_txt_file == F) {
temp_seurat_10X.data <-
  Read10X(temp_raw_matrix)
}

if(temp_matrix_txt_file) {
  temp_seurat_10X.data <-
    read.table(temp_raw_matrix, 
               sep = "\t", 
               row.names = 1, 
               header = T)
}

# add unique barcode names
colnames(x = temp_seurat_10X.data) <- 
  paste(temp_project_name, 
        colnames(x = temp_seurat_10X.data), 
        sep = '_')

# remove genome prefix from gene names from multi species alignments
if(temp_species_type == "human") {
  rownames(x = temp_seurat_10X.data) <- gsub("GRCh38_", 
                                             "", 
                                             row.names(temp_seurat_10X.data))
}

if(temp_species_type == "mouse") {
  rownames(x = temp_seurat_10X.data) <- gsub("mm10___", 
                                             "", 
                                             row.names(temp_seurat_10X.data))
}

saveRDS(temp_seurat_10X.data,
        "Output/Rdata/01_Read10X_raw_data.RData")

# DROPLET UTILS -----------------------------------------------------------

# vignette available at https://bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets

if(temp_use_droplet_utils_filtering){
  
set.seed(100)
temp_emptyDrops_out <- emptyDrops(temp_seurat_10X.data,
                                  lower = 250)

write.csv(temp_emptyDrops_out, 
          "Output/EmptyDrops/01_Emptydrops_out.csv")

temp_is.cell <- 
  temp_emptyDrops_out$FDR <= temp_FDR_threshold

sum(temp_is.cell,
    na.rm=TRUE)

temp_emptyDrops_out_subset <- 
  subset(temp_emptyDrops_out, 
         !is.na(FDR))

temp_cell_ids <-
  rownames(temp_emptyDrops_out_subset[temp_emptyDrops_out_subset$FDR <= 0.01,])

write.csv(data.frame(cell_ids = temp_cell_ids),
          "Output/EmptyDrops/02_emptydrops_filtered_cell_ids.csv")

temp_permutation_table <-
  table(Limited=temp_emptyDrops_out$Limited, 
        Significant=temp_is.cell)

write.csv(temp_permutation_table,
          "Output/EmptyDrops/03_permutation_table.csv")

temp_png_function("Output/EmptyDrops/04_Emptydrops_cells_filtered_deviation_plot.png")
plot(temp_emptyDrops_out$Total, 
     -temp_emptyDrops_out$LogProb, 
     col=ifelse(temp_is.cell, 
                "red", 
                "black"),
     xlab="Total UMI count", 
     ylab="-Log Probability")
dev.off()

}

# CREATE SEURAT OBJECT -----------------------------------------------------

seurat_10X <-  CreateSeuratObject(
  counts = temp_seurat_10X.data,
  min.cells = 1, 
  project = temp_project_name,
  min.features = 0
)

# mito content

temp_mito.genes <- grep(
  pattern = (if (temp_species_type == "human")
    "MT-"
    else
      "mt-"),
  x = rownames(x = seurat_10X),
  value = TRUE
)

temp_percent.mito <-
  Matrix::colSums(seurat_10X[temp_mito.genes, ]) / Matrix::colSums(seurat_10X)

seurat_10X@meta.data$percent.mito <- temp_percent.mito

# RAW SUMMARY STATS ---------------------------------------------

# barcode rank plot
temp_barcode_rank_plot <- data.frame(cell_barcode = row.names(seurat_10X@meta.data),
                                     nUMI = seurat_10X@meta.data$nCount_RNA)

temp_barcode_rank_plot <- temp_barcode_rank_plot[order(-temp_barcode_rank_plot$nUMI),]

temp_barcode_rank_plot <- cbind(temp_barcode_rank_plot,
                                barcode_number = c(1:length(temp_barcode_rank_plot$cell_barcode)),
                                sample = temp_project_name)

temp_ggplot <- ggplot(temp_barcode_rank_plot,
                      aes(x=barcode_number,
                          y = nUMI)) + geom_line() + 
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000)) +
  scale_y_continuous(trans = "log10",
                     breaks = c(50, 100, 500, 1000, 2000, 5000, 10000, 50000, 100000))
temp_png_function_small("Output/Figures/Metrics/01_barcode_rank_plot.png")
print(temp_ggplot)
dev.off()


# filter noise for raw plotting
seurat_10X <-  subset(x = seurat_10X, 
                      subset = nCount_RNA > 5)
seurat_10X <-  subset(x = seurat_10X, 
                      subset = nFeature_RNA > 5)
seurat_10X_raw <- seurat_10X
saveRDS(seurat_10X,
        "Output/Rdata/02_seurat_object_unprocessed.RData")

options(digits=1)
temp_raw_data_summary <-
  as.data.frame(cbind((as.matrix(summary(
    seurat_10X@meta.data$nFeature_RNA
  ))),
  (as.matrix(summary(
    seurat_10X@meta.data$nCount_RNA
  ))),
  (as.matrix(
    summary(seurat_10X@meta.data$percent.mito)
  ))))

colnames(temp_raw_data_summary) <-
  c("#GENES",
    "#UMIS",
    "#MITO%")

temp_raw_data_summary$CELLNUMBER <- ncol(seurat_10X)
temp_raw_data_summary$CELLNUMBER[2:6] <- "NA" 

write.csv(temp_raw_data_summary, 
          "Output/Gene_expression_and_stats/01_RAW_METRICS.csv", 
          quote = F, 
          row.names = T)

# FILTERING  ---------------------------------------------------------

if(temp_use_droplet_utils_filtering) {

  # emptydrops filtering
  seurat_10X <- 
    subset(seurat_10X, 
           cells = temp_cell_ids)
  
  seurat_10X <-  subset(x = seurat_10X, 
                        subset = nFeature_RNA > 200)
  
  seurat_10X <-  subset(x = seurat_10X, 
                        subset = percent.mito < temp_percent.mito.high_threshold)
  
}

if(temp_use_droplet_utils_filtering == F){

#filter seurat object
  seurat_10X <-  subset(x = seurat_10X, 
                        subset = nCount_RNA > temp_nUMI.low_threshold
                        & nCount_RNA < temp_nUMI.high_threshold &
                          nFeature_RNA > temp_nGene.low_threshold &
                          nFeature_RNA < temp_nGene.high_threshold &
                          percent.mito < temp_percent.mito.high_threshold)

}

# FILTERED SUMMARY STATS AND PLOTS ----------------------------------------

# append filtered status to raw object for plotting
temp_df <- data.frame(barcode = colnames(seurat_10X),
                      status = "filtered")
temp_raw_barcodes <- colnames(seurat_10X_raw)[!colnames(seurat_10X_raw) %in% colnames(seurat_10X)]

if(!length(temp_raw_barcodes) == 0) {
temp_dfB <- data.frame(barcode = temp_raw_barcodes,
                       status = "filtered_out")
temp_dfM <- rbind(temp_df, temp_dfB)
temp_dfM <- data.frame(row.names = temp_dfM$barcode,
                       status = temp_dfM$status)
temp_dfM_sorted <- 
  temp_dfM[colnames(seurat_10X_raw),,drop=FALSE]

temp_dfM_sorted <- 
  as.data.frame.matrix(temp_dfM_sorted)

seurat_10X_raw <- AddMetaData(seurat_10X_raw,
                              metadata = temp_dfM_sorted)
}

if(length(temp_raw_barcodes) == 0) {
  seurat_10X_raw@meta.data$status <- "filtered"
}

# filtered metrics
temp_filtered_data_summary <-
  as.data.frame(cbind((as.matrix(summary(
    seurat_10X@meta.data$nFeature_RNA
  ))),
  (as.matrix(summary(
    seurat_10X@meta.data$nCount_RNA
  ))),
  (as.matrix(
    summary(seurat_10X@meta.data$percent.mito)
  ))))

colnames(temp_filtered_data_summary) <-
  c("#GENES",
    "#UMIS",
    "#MITO%")

temp_filtered_data_summary$CELLNUMBER <- ncol(seurat_10X)
temp_filtered_data_summary$CELLNUMBER[2:6] <- "NA" 

write.csv(temp_filtered_data_summary, 
          "Output/Gene_expression_and_stats/02_EMPTYDROPS_FILTERED_METRICS.csv", 
          quote = F, 
          row.names = T)

Idents(object = seurat_10X_raw) <- "status"

# Plots
temp_png_function("Output/Figures/Metrics/02_FILTERED_VlnPlot_nGene_nUMI_percent_mito.png")
VlnPlot(
  object = seurat_10X_raw,
  features = c("nFeature_RNA",
               "nCount_RNA",
               "percent.mito"),
  cols = 3,
  split.by = "status"
)
dev.off()

temp_png_function("Output/Figures/Metrics/03_FILTERED_GenePlot_nGene_nUMI.png")
FeatureScatter(object = seurat_10X_raw, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA",
               pt.size = 1,
               group.by = "status")
dev.off()
temp_png_function("Output/Figures/Metrics/04_FILTERED_GenePlot_nUMI_percent_mito.png")
FeatureScatter(object = seurat_10X_raw, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mito",
               pt.size = 1, 
               group.by = "status")
dev.off()
temp_png_function("Output/Figures/Metrics/05_FILTERED_GenePlot_nGene_percent_mito.png")
FeatureScatter(object = seurat_10X_raw, 
               feature1 = "nFeature_RNA", 
               feature2 = "percent.mito",
               pt.size = 1, 
               group.by = "status")
dev.off()

# NORMALISATION, FINDVARIABLEGENES AND SCALING ---------------------------------------------------

#log normalisation
seurat_10X <-
  NormalizeData(
    object = seurat_10X,
    display.progress = F
  )

seurat_10X <- FindVariableFeatures(
  object = seurat_10X,
  do.plot = F
)

seurat_10X <- ScaleData(
  object = seurat_10X,
  vars.to.regress = temp_vars.to.regress,
  display.progress = F
)

# PCA ANALYSIS ------------------------------------------------------------

seurat_10X <- RunPCA(
  object = seurat_10X,
  npcs = 100
)

# JACKSTRAW ANALYSIS (NOT SETUP IN V3 DEVELOPERS VERSION) ------------------------------------------------------

# if(temp_jackstraw_cutoff) {
#   
#   temp_jackstraw_start <- print(date())
#   
#   seurat_10X <- JackStraw(
#     object = seurat_10X,
#     dims = 1:temp_PCs_to_compute,
#     reduction = "pca"
#   )
#   
#   temp_jackstraw_finish <- print(date())
#   
# }

# DIMENSIONAL REDUCTION TSNE & UMAP  --------------------------------------------------------------

# TSNE
for(i in c("A", "B", "C")) {
  
  temp_PC <- get(paste0("temp_PC_",i))
  seurat_10X <-
    RunTSNE(object = seurat_10X,
            dims = temp_PC,
            reduction.key = paste0("TSNE",i,"_"),
            reduction.name = paste0("TSNE",i))
}

# UMAP
for(i in c("A", "B", "C")) {
  
  temp_PC <- get(paste0("temp_PC_",i))
  seurat_10X <- 
    RunUMAP(seurat_10X, 
            dims = temp_PC,
            reduction.key = paste0("UMAP",i,"_"),
            reduction.name = paste0("UMAP",i),
            verbose = F
    )
}


for(i in c("A", "B", "C")) {
  for(dr in c("TSNE", "UMAP")) {
    
  temp_png_function(paste0("Output/Figures/",dr,"/",dr,"_PC_",
                           i,
                           ".png"))
  temp_dimplot <- DimPlot(
    object = seurat_10X,
    label.size = 4,
    pt.size = 2,
    reduction = paste0(dr,i)
  )
  print(temp_dimplot)
  dev.off()
  }
  }


# READ INPUT GENE SETS -----------------------------------------------------------------

# detected genes in seurat object
temp_genes.detected <- rownames(seurat_10X)

# read in gene input file
temp_gene_plot <- read.csv(temp_gene_plot_file)

# read marker gene sets
for(i in c(1:ncol(temp_gene_plot))){
  n <- paste0("temp_markers_",
              i)
  assign(n,
         na.omit(temp_gene_plot[i]))
  
  m <- paste0("temp_markers_set_name_",
              i)
  assign(m,
         colnames(temp_gene_plot[i]))
}

# species conversion 
# mouse
if (temp_species_type == "mouse") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_convert <- 
      get(paste0("temp_markers_",
                 i))
    colnames(temp_genes_to_convert) <- "gene"
    
    n <- paste0("temp_mouse_markers_",
                i)
    assign(n,
           firstup(tolower(temp_genes_to_convert$gene))
    )
  }
}
# human
if (temp_species_type == "human") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_convert <- 
      get(paste0("temp_markers_",
                 i))
    colnames(temp_genes_to_convert) <- "gene"
    
    n <- paste0("temp_human_markers_",
                i)
    assign(n,
           temp_genes_to_convert$gene)
  }
}

# filter gene sets for genes that are detected
# mouse
if (temp_species_type == "mouse") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_filter <- 
      get(paste0("temp_mouse_markers_",
                 i))
    temp_genes_to_filter <- 
      data.frame(genes = temp_genes_to_filter)
    temp_genes_to_filter <- 
      subset(temp_genes_to_filter,
             genes %in% temp_genes.detected)
    n <- paste0("temp_mouse_markers_",
                i)
    assign(n,
           temp_genes_to_filter)
  }
}

# human
if (temp_species_type == "human") {
  for(i in c(1:ncol(temp_gene_plot))) {
    temp_genes_to_filter <- 
      get(paste0("temp_human_markers_",
                 i))
    temp_genes_to_filter <- 
      data.frame(genes = temp_genes_to_filter)
    temp_genes_to_filter <- 
      subset(temp_genes_to_filter,
             genes %in% temp_genes.detected)
    n <- paste0("temp_human_markers_",
                i)
    assign(n,
           temp_genes_to_filter)
  }
}

# EXPRESSION OF GENE MARKERS AND METRICS --------------------------------------------------

for(j in c(1:ncol(temp_gene_plot))) {
    for(i in c("A", "B", "C")) {
    
    temp_markers_set_name <- get(paste0("temp_markers_set_name_",
                                        j))
    
    temp_markers_to_use <- get(paste0("temp_",
                                      temp_species_type,
                                      "_markers_",
                                      j))
    
    temp_markers_to_use <- (temp_markers_to_use$genes)
    temp_markers_to_use <- levels(droplevels(temp_markers_to_use))
    
      for(dr in c("TSNE", "UMAP")) {
      
      temp_png_function_big(paste0("Output/Figures/Gene_markers/",dr,"/FeaturePlot_0",
                                   j, 
                                   "_",
                                   temp_markers_set_name,
                                   "_markers_", 
                                   i,
                                   ".png"))
    
      temp_featureplot <- FeaturePlot(
      object = seurat_10X,
      features = temp_markers_to_use,
      pt.size = 1.5,
      reduction = paste0(dr,i))
      print(temp_featureplot)
      dev.off()
      }
  }
}

for(i in c("A", "B", "C")) {
  temp_number <- ncol(temp_gene_plot)+1
  for(dr in c("TSNE", "UMAP")) {
  temp_png_function(paste0("Output/Figures/Gene_markers/",dr,
                           "/", "FeaturePlot_0", temp_number,"_gene_UMI_mito_", 
                           i,
                           ".png"))
  temp_featureplot <- FeaturePlot(
    object = seurat_10X,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mito"),
    pt.size = 1.5,
    reduction = paste0(dr,i)
  )
  print(temp_featureplot)
  dev.off()
  }
  }

# CLUSTERING DIFFERENTIAL EXPRESSION & HEATMAPS --------------------------------------------------------------

temp_PC_res_dataframe <- 
  data.frame(row.names = row.names(seurat_10X@meta.data))

for(i in c("A", "B", "C")) {
  for(j in c(temp_res_1, temp_res_2, temp_res_3)) {
    
    temp_res_use <- as.numeric(gsub("res.",
                                    "",
                                    j))
    temp_dims.use <- get(paste0("temp_PC_",i))
    
    temp_graph_name <- paste0("PC_", i)
    
    seurat_10X <- 
      FindNeighbors(object = seurat_10X,
                    dims = temp_dims.use,
                    graph.name = temp_graph_name
      )
    
    seurat_10X <-
      FindClusters(
        object = seurat_10X,
        resolution = temp_res_use,
        graph.name = temp_graph_name
      )
    
    temp_colname <- paste0(temp_graph_name, "_res.", temp_res_use)
    
    temp_PC_res_dataframe <- cbind(temp_PC_res_dataframe,
                                   (seurat_10X@meta.data[,temp_colname]))
    
    colnames(temp_PC_res_dataframe)[(ncol(temp_PC_res_dataframe))] <- temp_colname
    
    for(dr in c("TSNE", "UMAP")) {
      temp_png_function(paste0("Output/Figures/",dr,"/",dr,"_PC_",
                               i, "_", "res", temp_res_use,
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X,
        label.size = 4,
        pt.size = 2,
        reduction = paste0(dr,i),
        group.by = temp_colname,
        label = T
      )
      print(temp_dimplot)
      dev.off()
    }
    
    if(temp_do_differential_expression) {
      
      Idents(object = seurat_10X) <- temp_colname
      
      options(digits = 4)
      temp_cluster.allmarkers <- FindAllMarkers(
        only.pos = T,
        object = seurat_10X,
        min.pct = temp_min.pct, 
        logfc.threshold = temp_thresh.use,
        min.diff.pct = temp_min.diff.pct, 
        test.use = 'MAST', 
        print.bar = F
      )
      
      temp_cluster.allmarkers <- arrange(temp_cluster.allmarkers,
                                         (cluster),
                                         desc(avg_logFC))
      
      write.csv(temp_cluster.allmarkers,
                file = paste0("Output/Gene_expression_and_stats/FindAllMarkers_PC_",
                              i,
                              "_res.",
                              temp_res_use, 
                              ".csv"))
      
      temp_genes_for_heatmap <- 
        (temp_cluster.allmarkers %>% group_by(cluster) %>% top_n(10,
                                                                 avg_logFC))
      
      temp_cluster_order <- sort(unique(seurat_10X@active.ident))
      
      temp_png_function_big(paste0("Output/Figures/Heatmaps/DoHeatMap_PC_",
                               i,
                               "_res.",
                               temp_res_use,
                               ".png"))
      print(DoHeatmap(
        object = seurat_10X,
        features = temp_genes_for_heatmap$gene,
        group.by = temp_colname,
        size = 0.3
      ))
      dev.off()
    }
  }
}

seurat_10X <- 
  AddMetaData(seurat_10X,
              metadata = temp_PC_res_dataframe)

# GARNETT CELL CLASSIFIER -------------------------------------------------

if(temp_run_garnett) {
  # load classifier
  temp_classifier <- readRDS("/share/ClusterShare/software/contrib/CTP_single_cell/sunwu_sc_pipeline/garnett_classifier_v1_stromalCCA5.Rdata")
  
  # generate monocle CDS object
  temp_raw_data <- GetAssayData(object = seurat_10X, 
                                slot = "counts")
  temp_raw_data <- as.data.frame(temp_raw_data)
  temp_raw_data <- as(as.matrix(temp_raw_data), 
                      "sparseMatrix")
  
  if (temp_species_type == "mouse") {
    rownames(temp_raw_data) <- toupper(row.names(temp_raw_data))
    temp_raw_data <- temp_raw_data[!duplicated(rownames(temp_raw_data)),]
  }
  
  pd <- new("AnnotatedDataFrame", 
            data = seurat_10X@meta.data)
  
  fData <- data.frame(gene_short_name = row.names(temp_raw_data), 
                      row.names = row.names(temp_raw_data))
  fd <- new("AnnotatedDataFrame", data = fData)
  
  lowerDetectionLimit <- 0
  if(all(temp_raw_data == floor(temp_raw_data))) {
    expressionFamily <- negbinomial.size()
  } else if(any(data < 0)){
    expressionFamily <- uninormal()
  } else {
    expressionFamily <- tobit()
  }
  
  monocle_cds <- newCellDataSet(temp_raw_data,
                                phenoData = pd, 
                                featureData = fd,
                                lowerDetectionLimit=lowerDetectionLimit,
                                expressionFamily=expressionFamily)
  monocle_cds <- estimateSizeFactors(monocle_cds)
  
  # classify cells with garnett
  monocle_cds <- classify_cells(monocle_cds,
                                temp_classifier,
                                db = org.Hs.eg.db,
                                cluster_extend = TRUE,
                                cds_gene_id_type = "SYMBOL")
  
  # append to seurat object
  head(pData(monocle_cds))
  temp_garnett_metadata <- pData(monocle_cds)
  temp_garnett_metadata_sorted <- 
    temp_garnett_metadata[rownames(seurat_10X@meta.data),,drop=FALSE]
  temp_garnett_metadata_sorted <- 
    as.data.frame.matrix(temp_garnett_metadata_sorted)
  temp_garnett_metadata_sorted <- data.frame(row.names = row.names(temp_garnett_metadata_sorted),
                                             garnett_call = temp_garnett_metadata_sorted$cell_type,
                                             garnett_call_ext = temp_garnett_metadata_sorted$cluster_ext_type)
  seurat_10X <- AddMetaData(seurat_10X, 
                            metadata = temp_garnett_metadata_sorted)
  
  # Assigning highest cell type call per cluster at highest calculated resolution 
  temp_res_use <- as.numeric(gsub("res.",
                                  "",
                                  temp_res_1))
  temp_colname <- paste0("PC_C", "_res.", temp_res_use)
  temp_table <- seurat_10X@meta.data
  temp_new_idents <- data.frame(row.names = row.names(temp_table),
                                original_garnett_call = temp_table[,"garnett_call"],
                                old_ident = temp_table[,temp_colname],
                                new_ident = NA)
  temp_idents <- unique(seurat_10X@meta.data[,temp_colname])
  temp_idents <- levels(droplevels(temp_idents))
  
  for(i in c(temp_idents)) {
    
    temp_table_subset <- temp_table[temp_table[,temp_colname] == i,]
    temp_subset_table <- as.data.frame(table(temp_table_subset$garnett_call))
    # drop Unknown 
    temp_subset_table <- temp_subset_table[!temp_subset_table$Var1 == "Unknown",]
    temp_id <- temp_subset_table$Var1[which.max(temp_subset_table$Freq)]
    temp_id <- levels(droplevels(temp_id))
    temp_new_idents$new_ident[temp_new_idents$old_ident == i] <- temp_id
    # keep all cycling cells
    temp_new_idents$new_ident[temp_new_idents$original_garnett_call ==  "T-cells cycling"] <- "T-cells cycling"
    temp_new_idents$new_ident[temp_new_idents$original_garnett_call == "Epithelial cycling"] <- "Epithelial cycling"
  }
  
  temp_new_idents_append <- data.frame(row.names = row.names(temp_table),
                                       garnett_cluster_call = temp_new_idents$new_ident)
  seurat_10X <- AddMetaData(seurat_10X, 
                            metadata = temp_new_idents_append)

  # DR plots cluster called
  for(i in c("A", "B", "C")) {
    for(dr in c("TSNE", "UMAP")) {
    temp_png_function(paste0("Output/Figures/Cell_type_ID_Garnett/",dr,"_PC_",
                             i, "_celltype_garnett", 
                             ".png"))
    temp_dimplot <- DimPlot(
      object = seurat_10X,
      label.size = 4,
      pt.size = 2,
      reduction = paste0(dr,i),
      group.by = "garnett_cluster_call",
      label = T
    )
    print(temp_dimplot)
    dev.off()
    }
  }
  
  # NON-CLUSTER CALLED DR plots
  for(i in c("A", "B", "C")) {
    for(dr in c("TSNE", "UMAP")) {
      temp_png_function(paste0("Output/Figures/Cell_type_ID_Garnett/01_raw_garnett_calls/",dr,"_PC_",
                               i, "_celltype_garnett_raw", 
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X,
        label.size = 4,
        pt.size = 2,
        reduction = paste0(dr,i),
        group.by = "garnett_call",
        label = T
      )
      print(temp_dimplot)
      dev.off()
    }
  }
  
  # NON-CLUSTER CALLED DR plots
  for(i in c("A", "B", "C")) {
    for(dr in c("TSNE", "UMAP")) {
      temp_png_function(paste0("Output/Figures/Cell_type_ID_Garnett/02_cluster_call_extended/",dr,"_PC_",
                               i, "_celltype_garnett_ext", 
                               ".png"))
      temp_dimplot <- DimPlot(
        object = seurat_10X,
        label.size = 4,
        pt.size = 2,
        reduction = paste0(dr,i),
        group.by = "garnett_call_ext",
        label = T
      )
      print(temp_dimplot)
      dev.off()
    }
  }
  
  # cell type proportions plot
  temp_table <- seurat_10X@meta.data
  temp_cellprop_df <-
    data.frame(unclass(table(temp_table$garnett_cluster_call)))
  temp_cellprop_df_1 <- data.frame(celltype = row.names(temp_cellprop_df),
                                   value = temp_cellprop_df[,1])
  temp_cellprop_df_1 <- cbind(temp_cellprop_df_1,
                              proportion = (temp_cellprop_df_1$value/sum(temp_cellprop_df_1$value)))

  temp_png_function("Output/Figures/Cell_type_ID_Garnett/01_Celltype_garnett_proportions.png")
  temp_ggplot <- ggplot(temp_cellprop_df_1, 
                        aes(fill=celltype,
                            x=celltype,
                            y=proportion)) +
    geom_bar(position="dodge", 
             stat="identity") +
    scale_y_continuous(labels = scales::percent, breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) + 
    ylab("Percentage") +
    xlab("Cell Type") +
    theme_bw(base_size = 16) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    guides(fill=guide_legend(title="Cell Type"))  
  print(temp_ggplot)
  dev.off()
  
  # write table
  temp_cellprop_df_1$proportion <- temp_cellprop_df_1$proportion*100
  write.csv(temp_cellprop_df_1, 
            "Output/Gene_expression_and_stats/03_celltype_garnett_proportions.csv",
            quote = F,
            row.names = F)  
  
  # markers for VlnPlot
  temp_genes_list <- list(Epithelial = c("EPCAM", "KRT18", "KRT8", "CD24", "KRT19"),
                          Mature_epithelial = c("ESR1", "FOXA1", "GATA3", "PGR", "ELF3"),
                          Myoepithelial = c("KRT5", "KRT14", "KRT17", "ACTA2", "OXTR", "CNN1"),
                          Endothelial = c("PECAM1", "CD34", "VWF", "FLT1", "ACKR1"),
                          T_cell = c("CD8A", "CD8B", "NKG7", "GZMA", "GZMK","CD4", "CD40LG"), 
                          T_Reg = c("FOXP3", "BATF", "TNFRSF4", "TNFRSF18", "CTLA4"),
                          Myeloid = c("CD68", "CD14", "FCER1G", "FCGR3A", "CD86"),
                          B_cell = c("JCHAIN", "IGKC", "MS4A1", "CD79A", "CD19", "TCL1A"),
                          CAF1 = c("PDGFRB", "DCN", "COL1A1", "COL1A2", "PDGFRA"),
                          CAF2 = c("MCAM", "ACTA2", "MYL9", "TAGLN", "CALD1"),
                          Proliferation = c("MKI67", "UBE2C", "AURKB", "TOP2A", "CENPF", "CDK1")
  )
  
  temp_names <- names(temp_genes_list)
  
  for(i in c(1:11)) {
    
    temp_markers <- temp_genes_list[[i]]
    if (temp_species_type == "mouse") {
      temp_markers <- firstup(tolower(temp_markers))
    }
    temp_markers <- temp_markers[temp_markers %in% temp_genes.detected]
    temp_output_name <- temp_names[i]
    
    temp_png_function_longer(paste0("Output/Figures/Cell_type_ID_Garnett/VlnPlot_classifier_genes_",
                               i,
                               "_",
                               temp_output_name,
                               ".png"))
    
    temp_featureplot <- VlnPlot(
      object = seurat_10X,
      features = temp_markers,
      pt.size = 2,
      group.by = "garnett_cluster_call",
      ncol = 3)
    
    print(temp_featureplot)
    dev.off()
  }
  
}

# SAVE PROCESSED OBJECT ---------------------------------------------------

seurat_10X <- SetAllIdent(seurat_10X, 
                          id = "RNA_snn_res.0.4")
saveRDS(seurat_10X,
        "Output/Rdata/03_seurat_object_processed.RData")

# EXPORTING EXPRESSION MATRIX ------------------------------------------------------

#Raw UMI matrix
if (temp_export_raw_UMI_matrix) {
  temp_matrix <- GetAssayData(object = seurat_10X, 
                              assay = "RNA", 
                              slot = "counts")
  
  temp_matrix.sparse <- Matrix(temp_matrix , sparse = T )
  
  write.table(
    temp_matrix.sparse,
    file = "Output/Matrices_and_metadata/normalized_expression_sparse_matrix.txt",
    sep = "\t",
    col.names = T
  )
}

#Normalised matrix
if (temp_export_normalised_matrix) {
  temp_matrix <- GetAssayData(object = seurat_10X, 
                              assay = "RNA", 
                              slot = "data")
  
  temp_matrix.sparse <- Matrix(temp_matrix, 
                               sparse = T )
  write.table(
    temp_matrix.sparse,
    file = "Output/Matrices_and_metadata/raw_UMI_count_sparse_matrix.txt",
    sep = "\t",
    col.names = T
  )
}

#Associated Metadata
  temp_table <- seurat_10X@meta.data
  temp_table_charactered <- apply(temp_table,2,as.character)
  temp_table_charactered <- cbind(barcode = row.names(temp_table),
                                  temp_table_charactered)
  
  
  write.table(temp_table_charactered,
              file = "Output/Matrices_and_metadata/metadata.txt",
              sep = "\t",
              col.names = T,
              row.names = F,
              quote = F
  )




# R MARKDOWN FILE ---------------------------------------------------------

# cellranger stats
temp_path_cellranger <- gsub("raw_gene_bc_matrices/GRCh38/",
                             "metrics_summary.csv",
                             temp_raw_matrix)
if(temp_species_type == "mouse") {
  temp_path_cellranger <- gsub("raw_gene_bc_matrices/mm10/",
                               "metrics_summary.csv",
                               temp_raw_matrix)
}

if(!file.exists(temp_path_cellranger)) {
  temp_path_cellranger <- gsub("raw_gene_bc_matrices/GRCh38_mm10/",
                               "metrics_summary.csv",
                               temp_raw_matrix)
  
  if(temp_species_type == "mouse") {
    temp_path_cellranger <- gsub("raw_gene_bc_matrices/GRCh38_mm10/",
                                 "metrics_summary.csv",
                                 temp_raw_matrix)
  }}
  
if(file.exists(temp_path_cellranger)) {
  temp_csv <- as.data.frame(read.csv(temp_path_cellranger,
                                     stringsAsFactors=FALSE))
  temp_csv <- as.data.frame(t(temp_csv))
  temp_csv$V1 <- (gsub(",", "",temp_csv$V1))
  rownames(temp_csv) <- gsub("\\."," ", rownames(temp_csv))
  colnames(temp_csv) <- c(temp_project_name)
  
  write.csv(temp_csv, 
            "Output/Analysis_params_and_Rmd/01_cellranger_stats.csv", 
            quote = F, 
            row.names = T)
}

# seurat analysis params
temp_params_save <- 
  read.csv(temp_params_file, row.names = "SEURAT_PROCESSING_PARAMS_FILE")

temp_summary <- file("Output/Analysis_params_and_Rmd/02_analysis_params.txt", 
                     open="wt")
sink(temp_summary, 
     type="output")
print(temp_params_save)
sink(type="output")
close(temp_summary)

# directory structure
temp_summary <- file("Output/Analysis_params_and_Rmd/03_output_directory_structure.txt", 
                     open="wt")
sink(temp_summary, 
     type="output")
twee("./Output")
sink(type="output")
close(temp_summary)

# R session info
temp_summary <- file("Output/Analysis_params_and_Rmd/04_Rsessioninfo.txt", 
                     open="wt")
sink(temp_summary, 
     type="output")
print(sessionInfo())
sink(type="output")
close(temp_summary)

# other metrics
temp_PC_A_length <- length(temp_PC_A)
temp_PC_B_length <- length(temp_PC_B)
temp_PC_C_length <- length(temp_PC_C)
temp_finish_time <- date()
temp_cell_number <- ncol(seurat_10X)
temp_total_genes <- nrow(seurat_10X)

# sink to Rmd file
# copy Rmd template file to output
file.copy("/share/ClusterShare/software/contrib/CTP_single_cell/sunwu_sc_pipeline/seurat_v3_rmarkdown.Rmd", 
          "./Output/")
rmarkdown::render("./Output/seurat_v3_rmarkdown.Rmd", 
                  output_file = paste0(temp_project_name,"_output.html"),
                  params = list(temp_project_name = temp_project_name,
                                temp_start_time = temp_start_time,
                                temp_finish_time = temp_finish_time,
                                temp_species_type = temp_species_type,
                                temp_cell_number = temp_cell_number,
                                temp_total_genes = temp_total_genes,
                                temp_PC_A_length = temp_PC_A_length,
                                temp_PC_B_length = temp_PC_B_length,
                                temp_PC_C_length = temp_PC_C_length,
                                temp_res_1 = temp_res_1,
                                temp_res_2 = temp_res_2,
                                temp_res_3 = temp_res_3,
                                temp_markers_set_name_1 = temp_markers_set_name_1,
                                temp_markers_set_name_2 = temp_markers_set_name_2,
                                temp_markers_set_name_3 = temp_markers_set_name_3,
                                temp_do_differential_expression = temp_do_differential_expression,
                                temp_run_garnett = temp_run_garnett
                                ))

# sink to pdf file (TO DO)

#clean environment
rm(list = ls(pattern = "temp"))
rm(seurat_10X)
