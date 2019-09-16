#### SEURAT V3 HPC R SCRIPT WITH DROPLET UTILS FILTERING
#     written by; Sunny Z Wu
#     last modified; 20190501
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
# qsub -cwd -pe smp 23 -l h_vmem=300G -b y -j y -V -P TumourProgression "${R} CMD BATCH --no-save '--args CID44041 human $(pwd) /share/ScratchGeneral/sunwu/Chromium_10X/count_data_2/CID4404/count_4404_primary_GRCh38/outs/raw_gene_bc_matrices/GRCh38/ ./seurat_gene_input_file.csv ./seurat_V3_HPC_params_file.csv /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/refs/' ./seurat_V3_HPC_with_garnett_classifier_and_infercnv.R"
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
# 04 PATH TO CELLRANGER RAW MATRIX
temp_raw_matrix <- temp_args[4]
# 05 GENES TO PLOT FILE
temp_gene_plot_file <- temp_args[5]
# 06 SEURAT PARAMS FILE
temp_params_file <- temp_args[6]
# 07 INFERCNV REFERENCE DIRECTORY
temp_infercnv_reference_dir <- temp_args[7]


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

# InferCNV params
temp_run_infercnv <- 
  as.logical(temp_params["Run_InferCNV",])
temp_infercnv_resolution <- as.character(temp_params["infercnv_resolution",])
temp_infercnv_cutoff <- 
  as.numeric(as.character(temp_params["infercnv_cutoff",])) # cuts genes with average expression 
                                # across all cells < cutoff
temp_infercnv_window <- 
  as.numeric(as.character(temp_params["infercnv_window",]))   # sliding window value
temp_infercnv_threshold <- 
  as.numeric(as.character(temp_params["infercnv_threshold",]))  # reduces final gene values > threshold
                                  # to threshold and < -threshold to
                                  # -threshold
temp_denoise_value <- 
  as.numeric(as.character(temp_params["infercnv_denoise_value",]))  # all gene values within 
                                    # denoise_value std devs of mean 
                                    # of values will be reduced to 0
temp_run_HMM <- 
  as.logical(temp_params["run_HMM_prediction",])  # Hidden Markov Model-based
                            # CNV calls in separate plot
temp_HMM_prediction_type <- 
  as.character(temp_params["HMM_prediction_type",]) # choose from 3 or 6 classes                        
                            # of CNV calls

print(temp_run_HMM)

if(temp_run_infercnv) {
  lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
  library(infercnv, lib.loc=lib_loc)
  library(HiddenMarkov, lib.loc=lib_loc)
  library(ComplexHeatmap, lib.loc = lib_loc)
  library(circlize, lib.loc = lib_loc)
  library(reshape2)
  library(grid)
  library(RColorBrewer)
}

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
dir.create("Output/InferCNV")
dir.create("Output/InferCNV/input_files")

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

### InferCNV functions ###

prepare_infercnv_metadata <- function(seurat_object, temp_resolution) {
  
  # fetch PC A, res 0.8 clusters
  Idents(seurat_object) <- 
    eval(parse(text = paste0("seurat_object@meta.data$", temp_resolution)))

  # create vector of Garnett annotated clusters
  annotated_idents <- gsub(
    " ", "_", paste0(seurat_object@meta.data$garnett_cluster_call, 
    " ", Idents(seurat_object))
  )
  # only keep annotation with largest no cells per existing cluster
  annotated_ident_numbers <- table(annotated_idents)
  annotated_ident_cluster_numbers <- gsub("^.*_", "", names(annotated_ident_numbers))

  for ( c in levels(Idents(seurat_object)) ) {
    annotated_ident_index <- which(annotated_ident_cluster_numbers == c)
    potential_clusters <- annotated_ident_numbers[annotated_ident_index]
    correct_cluster <- 
      names(potential_clusters[potential_clusters == max(potential_clusters)])
      # label cells in seurat object as the correct cluster
      levels(Idents(seurat_object))[levels(Idents(seurat_object)) == c] <- 
        correct_cluster
  }

  # create infercnv metadata
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
    cell_type = Idents(seurat_object), stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]
  # record number per cell type
  number_per_cell_type <- as.data.frame(table(temp_metadata$cell_type))

  return(list(temp_metadata, number_per_cell_type, seurat_object))
}

create_infercnv_object <- function(raw_counts_path, annotation_file_path, gene_order_path) {
  infercnv_obj = CreateInfercnvObject(
            raw_counts_matrix=raw_counts_path,
            annotations_file=annotation_file_path,
            delim="\t",
            gene_order_file=gene_order_path,
            ref_group_names=normals
  )
  return(infercnv_obj)
}

run_infercnv <- function(infercnv_object, outdir, cutoff, window, threshold, denoise_value,
  runHMM, typeHMM) {
  if (runHMM) {
      infercnv::run(infercnv_object,
            num_threads=25,
            out_dir=outdir,
            cutoff=cutoff,
            window_length=window,
            max_centered_threshold=threshold,
            cluster_by_groups=T,
            plot_steps=F,
            denoise=T,
            sd_amplifier=denoise_value,
            HMM=runHMM,
            HMM_type=typeHMM
      )
  } else {
    infercnv::run(infercnv_object,
            num_threads=25,
            out_dir=outdir,
            cutoff=cutoff,
            window_length=window,
            max_centered_threshold=threshold,
            cluster_by_groups=T,
            plot_steps=F,
            denoise=T,
            sd_amplifier=denoise_value
    )
  }
}

# function to call ggplot default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

create_group_annotation <- function(df, metadata) {
  m <- match(rownames(df), metadata$cell_ids)
  df$groups <- metadata$cell_type[m]
  # sort df and metadata by cluster number
  df <- df[order(
    as.numeric(
        gsub("^.*_", "", df$groups)
    )
  ),]

  df$groups <- gsub("_", " ", df$groups)
  # store groups column in separate data frame and remove from data frame:
  group_annotation_df <- data.frame(group = df$groups, row.names = rownames(df),
      stringsAsFactors = F)

  # define group heatmap colours:
    cluster_nos <- as.numeric(gsub("^.*\\ ", "", unique(group_annotation_df$group)))
    cols <- gg_color_hue(max(cluster_nos+1))[cluster_nos+1]
    names(cols) <- unique(group_annotation_df$group)
    
    # create heatmap of group annot:
    group_annotation <- Heatmap(group_annotation_df, col = cols, 
      name = "group_annotation", 
      width = unit(4, "mm"), 
      show_row_names = F, show_column_names = F,
      heatmap_legend_param = list(title = "Cluster", 
      title_gp = gpar(fontsize = 8, fontface = "bold"), 
      labels_gp = gpar(fontsize = 8)))
    result_list <- list(group_annotation_df, group_annotation)
    names(result_list) <- c("group_annotation_df", "group_annotation")

    return(result_list)
}

# fetch co-ordinates of chromosomal boundaries:
fetch_chromosome_boundaries <- function(df) {
  # load in gene annotations:
  gene_order <- read.table(paste0(temp_infercnv_reference_dir, "infercnv_gene_order.txt"), 
    header=F, as.is=T)
  # subset to only include genes in df:
  genes <- gene_order[gene_order$V1 %in% colnames(df),1:2]
  # remove genes on chrM and chrY from df:
  genes_to_keep <- genes$V1[genes$V2 == "chrM|chrY"]
  df <- df[,!(colnames(df) %in% genes_to_keep)]
  # split and determine lengths of each chromosome:
  chr_split <- split(genes$V2, genes$V2)
  chr_lengths <- unlist(lapply(chr_split, length))
  chr_lengths <- chr_lengths[!(names(chr_lengths) %in% c("chrY", "chrM"))]
  # order chromosomes:
  chr_lengths <- c(
    chr_lengths[
      order(
        as.numeric(
          as.character(
            gsub(
              "chr", "", names(chr_lengths[!(names(chr_lengths) %in% c("chrX"))])
            )
          )
        )
      )
      ], chr_lengths[names(chr_lengths) %in% c("chrX")])
  
  for ( l in 1:length(chr_lengths) ) {
    if (l==1) {
      chr_ends <- c(chr_lengths[l])
    } else {
      chr_ends[l] <- chr_ends[l-1] + chr_lengths[l]
    }
  }
  names(chr_ends) <- names(chr_lengths)
  # find total length of plot:
  total_length <- ncol(df)
  # define vertical line co-ordinates for each chromosomal end:
  end_pos <- chr_ends/total_length
  # define chromosome label co-ordinates for midpoint of each chromosome:
  # find centre of each chromosome:
  for ( i in 1:length(chr_lengths) ) {
    if (i==1) {
      lab_pos <- c(chr_lengths[i]/2)
    } else {
      lab_pos[i] <- chr_ends[i-1] + (chr_lengths[i]/2)
    }
  }
  lab_pos <- lab_pos/chr_ends[length(chr_ends)]
  names(lab_pos) <- names(chr_lengths)
  result_list <- list(lengths = chr_lengths, ends = chr_ends, end_pos = end_pos, 
    lab_pos = lab_pos)

  return(result_list)
}

# function to create nUMI and nGene barplot annotations:
create_QC_annotation <- function(seurat_object, df) {
  
  # load seurat object and create QC metrics df:
  qc_df <- data.frame(
  seurat_object@meta.data$nCount_RNA,
  seurat_object@meta.data$nFeature_RNA,
  row.names = as.character(names(Idents(seurat_object))))
  colnames(qc_df) <-  c("nUMI", "nGene")

  nUMI_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      qc_df$nUMI,
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
      qc_df$nGene,
      gp = gpar(
        col = "#9ECAE1", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )

  result_list <- list(nUMI_annotation, nGene_annotation, qc_df)
  names(result_list) <- c("nUMI_annotation", "nGene_annotation", "qc_df")

  return(result_list)

}

# function to annotate InferCNV heatmap with genome-wide CNV frequencies 
# of PAM50 subtypes defined by METABRIC
annotate_PAM50_CNV <- function(observation_df, cnv_frequencies, temp_subtype, chr_ends,
  chr_lengths) {
    
    cnv <- cnv_frequencies[,c(1:4, 6, which(colnames(cnv_frequencies) %in% c(paste0("Loss", temp_subtype), 
      paste0("Gain", temp_subtype))))]
    
    # keep only protein-coding genes:
    cnv <- cnv[grep("protein_coding", cnv$ATTRIBUTE),]
    
    # append meta_cnv info onto observation_df column to put it in order:
    temp <- as.data.frame(t(observation_df[1,]))
    temp$FEATURE <- rownames(temp)
    cnv_df <- merge(temp, cnv, by="FEATURE", all.x = T)
    
    # remove duplicates from cnv_df:
    cnv_df <- cnv_df[-which(duplicated(cnv_df$FEATURE)),]
    
    # order cnv_df by observation_df order:
    m <- match(colnames(observation_df), cnv_df$FEATURE)
    cnv_df <- cnv_df[m,]
    
    # remove unwanted columns:
    cnv_df <- cnv_df[,grep("CID|IND|PID|PDX|start|stop|ATTRIBUTE", colnames(cnv_df), invert=T)]
    
    # replace NAs with 0:
    ind <- grep("Loss", colnames(cnv_df))
    cnv_df[,ind][is.na(cnv_df[,ind])] <- 0
    ind <- grep("Gain", colnames(cnv_df))
    cnv_df[,ind][is.na(cnv_df[,ind])] <- 0
    
    for (j in 1:length(chr_lengths)) {
      if (j==1) {
        cnv_df$Chr[1:chr_lengths[j]] <- i
      } else {
        cnv_df$Chr[chr_lengths[j-1]:chr_lengths[j]] <- j
      }
    }
    
    # change '23' to 'X' in chromosome column:
    cnv_df$Chr[cnv_df$Chr==23] <- "X"
    
    # make loss values negative:
    ind <- grep("Loss", colnames(cnv_df))
    cnv_df[,ind] <- -cnv_df[,ind]
    
    # create area plot of CNVs:
    # adjust plot df:
    cnv_df$BPcum <- seq(1, nrow(cnv_df))
    
    m_df <- melt(cnv_df, id.vars = c("FEATURE", "Chr", "BPcum"))
    m_df <- m_df[order(m_df$BPcum),]
    
    # prepare df for labelling:
    colnames(m_df) <- c("Gene", "Chr", "Genomic_Region", "CNV_type", "Frequency")
    m_df$CNV_type <- gsub(temp_subtype, "", m_df$CNV_type)
    
    # define colours:
    cols <- gg_color_hue(2)
    
    # create area plot:
    p <- ggplot(m_df, aes(x=Genomic_Region, y=Frequency))
    p <- p + geom_area(aes(color=factor(CNV_type, levels = c("Gain", "Loss"))))
    p <- p + scale_fill_manual(c(cols[1], cols[2]))
    p <- p + scale_x_continuous(label = names( c(chr_ends, "") ), breaks = c(0, chr_ends),
                                expand = c(0,0))
    p <- p + theme_bw()
    p <- p + scale_y_continuous(name=paste0(temp_subtype), breaks = c(-1, 0, 1), 
                                limits = c(-1, 1))
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   text = element_text(size=6),
                   axis.title.y = element_text(size=8, angle=0, vjust = 0.5),
                   panel.border = element_blank(),
                   panel.grid.major.x = element_line(colour = "#383838", size = 0.2),
                   panel.grid.minor.x = element_blank(),
                   legend.position="none")
    return(grid.grabExpr({print(p)}))
}

create_CNV_genes_annotation <- function(heat_df, gene_df) {

    heatmap_CNV_genes <- data.frame(gene = colnames(heat_df))

    # isolate CNV-associated genes in heatmap df:
    gene_df_sub <- gene_df[gene_df$Hugo_symbol %in% heatmap_CNV_genes$gene,]

    # if >30 genes, reduce to 30, gains and losses are prioritised:
    gene_df_sub <- gene_df_sub[1:30,]

    heatmap_CNV_genes$type <- "do_not_include"
    for ( t in c("gain", "loss", "either")) {
      heatmap_CNV_genes$type[heatmap_CNV_genes$gene %in% 
      gene_df_sub$Hugo_symbol[gene_df_sub$type == t]] <- t
    }

    label_subset <- which(heatmap_CNV_genes$type != "do_not_include")
    heatmap_CNV_labels <- as.character(heatmap_CNV_genes$gene[label_subset])

    # expand annotation bands by including a number of genes around GOI depending
    # on plot width:
    all_do_not_include <- which(heatmap_CNV_genes$type != "do_not_include")
    expand_no <- as.integer(ncol(heatmap_df)/1650)
    for (n in 1:expand_no) {
        for ( temp_index in all_do_not_include ) {
        
          heatmap_CNV_genes$type[temp_index + n] <- heatmap_CNV_genes$type[temp_index]
        
          if ( (temp_index - n) > 0 ) {
              heatmap_CNV_genes$type[temp_index - n] <- heatmap_CNV_genes$type[temp_index]
          }
        }
    }
    # create CNV_genes text annotation vector:
    heatmap_CNV_genes <- subset(heatmap_CNV_genes, select = -gene)

    CNV_genes_annot = HeatmapAnnotation(df = heatmap_CNV_genes, show_legend = F,
        col = list(type = c("gain" = "#F2210E", "loss" = "#3691BC", 
        "either"="black", "do_not_include"="#E7E4D3")),
        text = anno_link(label_subset, heatmap_CNV_labels, side = "bottom", 
        labels_gp = gpar(fontsize = 6, rot=2)))

    return(CNV_genes_annot)
}

if (!file.exists("Output/Rdata/03_seurat_object_processed.RData")) {

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

} else {
  seurat_10X <- readRDS("Output/Rdata/03_seurat_object_processed.RData")
}


# INFERCNV ------------------------------------------------------

if (temp_run_infercnv) {

  # prepare infercnv metadata and write to files
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_10X, temp_infercnv_resolution)
  names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
  seurat_10X <- infercnv_metadata$seurat

  # if no epithelial or CAF clusters present, run without normals:
  if ( length(grep("pithelial", infercnv_metadata$number_per_group$Var1)) < 1 & 
    length(grep("CAF", infercnv_metadata$number_per_group$Var1)) < 1) {
    print("No epithelial or CAF clusters detected, running with no 'normal' groups...")
  } else if ( length(grep("pithelial", infercnv_metadata$number_per_group$Var1)) < 1 ) {
    # if no epithelial clusters are detected, re-annotate CAF clusters as epithelial as
    # they are likely being mis-classified, if no CAFs abort run:
    print("No epithelial clusters detected, reclassifying CAFs as epithelial...")
    infercnv_metadata$metadata$cell_type <- as.character(infercnv_metadata$metadata$cell_type)
    CAF_metadata <- infercnv_metadata$metadata$cell_type[
      grep("CAF", infercnv_metadata$metadata$cell_type)
    ]
    new_metadata <- gsub(
      "Epithelial[1,2]", "Epithelial", gsub(
        "CAF", "Epithelial", CAF_metadata
      )
    )
    infercnv_metadata$metadata$cell_type[
      grep("CAF", infercnv_metadata$metadata$cell_type)
    ] <- new_metadata
    levels(Idents(seurat_10X)) <- gsub("CAF[1,2]", "Epithelial", levels(Idents(seurat_10X)))
    
    # touch file as indicator that CAFs were reclassified:
    system("touch Output/InferCNV/WARNING_CAF_clusters_were_reclassified_to_epithelial")
  }
  print(paste0("Idents are: ", levels(Idents(seurat_10X))))
  write.table(infercnv_metadata$metadata, "Output/InferCNV/input_files/metadata.txt", 
    quote=F, sep="\t", col.names=F, row.names=F)
  write.table(infercnv_metadata$number_per_group, 
    "Output/InferCNV/input_files/number_per_group.txt", 
    quote=F, col.names=F, row.names=F, sep="\t")
  saveRDS(seurat_10X, "Output/Rdata/04_seurat_object_annotated.RData")
  
  # define normals as stromal groups, excluding epithelial and unassigned cells:
  if ( length(grep("pithelial", infercnv_metadata$number_per_group$Var1)) < 1 & 
  length(grep("CAF", infercnv_metadata$number_per_group$Var1)) < 1) {
    normals <- NULL
  } else {
    normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
      levels(Idents(seurat_10X)), value=T, invert=T)
  }

  # create raw matrix input file
  print("Creating inferCNV raw counts file...")
  if (!file.exists("Output/InferCNV/input_files/raw_counts_matrix.txt")) {
    count_df <- as.matrix(GetAssayData(seurat_10X , slot = "counts"))
    #count_df <- as.matrix(seurat_10X@assays$integrated)
    write.table(count_df, 
    "Output/InferCNV/input_files/raw_counts_matrix.txt", quote=F,
    sep="\t", col.names=T, row.names=T)
  }

  # generate cluster metric plots for epithelial clusters:
  epithelial_clusters <- grep("pithelial", levels(Idents(seurat_10X)), value=T)
  temp_png_function("Output/InferCNV/metrics_by_epithelial_cluster.png")
    temp_violinplot <- VlnPlot(
      object = seurat_10X,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
      pt.size = 1.5,
      idents = epithelial_clusters
    )
    print(temp_violinplot)
  dev.off()

  print("Creating inferCNV object...")
  raw_path <- paste0(temp_wd, "/Output/InferCNV/input_files/raw_counts_matrix.txt")
  annotation_path <- paste0(temp_wd, "/Output/InferCNV/input_files/metadata.txt")
  gene_path <- paste0(temp_infercnv_reference_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- create_infercnv_object(raw_path, annotation_path, gene_path)

  print("InferCNV object created, running inferCNV...")
  out_dir <- "Output/InferCNV/supervised_clustering/"
  infercnv_output <- run_infercnv(initial_infercnv_object, out_dir, temp_infercnv_cutoff, 
    temp_infercnv_window, temp_infercnv_threshold, temp_denoise_value, temp_run_HMM, 
    temp_HMM_prediction_type)
  
  # load observations for customised heatmap
  if (temp_run_HMM) {
    infercnv_observations <- read.table(paste0(out_dir,  
      "/infercnv.15_denoised.observations.txt"), fill=T)
  } else {
    infercnv_observations <- read.table(paste0(out_dir, 
      "/infercnv.12_denoised.observations.txt"), fill=T)
  }
  print("InferCNV run finished")


  # INFERCNV annotated heatmap ------------------------------------------------------
  # format observations for heatmap
  heatmap_df <- as.data.frame(t(infercnv_observations))
  # remove CAFs from heatmap df and metadata:
  cells_to_remove <- 
    rownames(infercnv_metadata$metadata)[grep("CAF", 
      infercnv_metadata$metadata$cell_type)]
  heatmap_df <- heatmap_df[!(rownames(heatmap_df) %in% cells_to_remove),]
  infercnv_metadata$metadata <- 
    infercnv_metadata$metadata[!(infercnv_metadata$metadata$cell_ids %in% 
    cells_to_remove),]
  # remove 'luminal' from epithelial cells as Garnett labels everything luminal:
  infercnv_metadata$metadata$cell_type <- gsub("Luminal_", "", 
  infercnv_metadata$metadata$cell_type)
  # create group annotation df
  group_annotation <- create_group_annotation(heatmap_df, infercnv_metadata$metadata)
  # ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
  m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
  heatmap_df <- heatmap_df[m,]
  # define chromosome lengths in terms of number of filtered genes
  chr_data <- fetch_chromosome_boundaries(heatmap_df)
  # create heatmap annotation for genome-wide PAM50 subtype CNV frequency
  PAM50_subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")
  METABRIC_CNV_frequencies <- 
  read.table(paste0(temp_infercnv_reference_dir, "infercnv_metabric_cnv.txt"), 
    header=T, as.is=T, fill=T)
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
  CNV_genes <- read.table(paste0(temp_infercnv_reference_dir, 
    "./infercnv_brca_genes_associated_with_CNVs.txt"), header = T, as.is = T)
  # create CNV_genes annotation:
  print("Annotating CNV-associated genes...")
  CNV_genes_annotation <- create_CNV_genes_annotation(heatmap_df, CNV_genes)
  
  # replace column labels with two letter label for positioning:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  print("Generating final heatmap...")
  
  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = colorRamp2(c(min(plot_object), 1, max(plot_object)), 
                     c("#00106B", "white", "#680700"), space = "sRGB"), 
    cluster_columns = F, cluster_rows = F,
    split = group_annotation$group_annotation_df$group,
    show_row_names = F, show_column_names = T,
    column_names_gp = gpar(col = "white"),
    show_row_dend = FALSE,
    bottom_annotation = CNV_genes_annotation, bottom_annotation_height = unit(2, "cm"),
    gap = unit(1, "cm"),
    heatmap_legend_param = list(title = "Modified\nexpression", color_bar = "continuous", 
    grid_height = unit(1.5, "cm"), legend_direction = "horizontal", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 6)),
    use_raster = T, raster_device = c("png")
  )
  # determine co-ordinates of horizontal lines at group borders:
  spl_groups <- split(group_annotation$group_annotation_df$group, 
  group_annotation$group_annotation_df$group)
  spl_groups <- spl_groups[unique(group_annotation$group_annotation_df$group)]
  if ( length(spl_groups) > 1 ) {
    for ( n in 1:(length(spl_groups)-1) ) {
      if (n==1) {
        hlines <- c(length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group))
      } else {
        hlines[n] <- hlines[n-1] + length(spl_groups[[n]])/length(group_annotation$group_annotation_df$group)
      }
    }
    hlines <- 1-hlines
  }

  annotated_heatmap <- grid.grabExpr(
    draw(group_annotation$group_annotation + final_heatmap, heatmap_legend_side = "left")
  )
  
  pdf("Output/InferCNV/supervised_clustering/infercnv_final_heatmap.pdf", 
  height = 8, width = 11)

    grid.newpage()
    
    # plot Normal subtype:
    pushViewport(viewport(x = 0.095, y = 0.090,
                          width = 0.9975-0.094, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[5]])
    popViewport()
    
    # plot Basal subtype:
    pushViewport(viewport(x = 0.103, y = 0.169,
                          width = 0.9975-0.102, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[4]])
    popViewport()
    
    # plot Her2 subtype:
    pushViewport(viewport(x = 0.107, y = 0.251, 
                          width = 0.9975-0.106, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[3]])
    popViewport()
    
    # plot LumB subtype:
    pushViewport(viewport(x = 0.102, y = 0.333, 
                          width = 0.9975-0.101, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[2]])
    popViewport()
    
    # plot LumA subtype:
    pushViewport(viewport(x = 0.102, y = 0.411, 
                          width = 0.9975-0.101, height = 0.08, just = c("left", "top")))
    grid.draw(metabric_plots[[1]])
    popViewport()
    
    # plot heatmap:
    pushViewport(viewport(x = 0.5, y = 0.4, 
                          width = 1, height = 0.6, just = "bottom"))
    grid.draw(annotated_heatmap)
    decorate_heatmap_body("hm", {
  
      for ( e in 1:length(chr_data$end_pos) ) {
        grid.lines(c(chr_data$end_pos[e], chr_data$end_pos[e]), c(0, 1), gp = gpar(lwd = 1, 
          col = "#383838"))
        grid.text(gsub("chr", "", names(chr_data$lab_pos)[e]), chr_data$lab_pos[e], 
          unit(0, "npc") + unit(-2.1, "mm"), gp=gpar(fontsize=8))
      }

      if ( length(spl_groups) > 1 ) {
        for ( m in 1:length(hlines) ) {
          grid.lines(c(0, 1), c(hlines[m], hlines[m]), gp = gpar(lwd = 1, col = "#383838"))
        }
      }
    })
  
    popViewport()

  dev.off()
	
	# convert pdf to png:
	system(paste0("convert -density 150 ", 
		"Output/InferCNV/supervised_clustering/infercnv_final_heatmap.pdf ", 
		"-quality 90 Output/InferCNV/supervised_clustering/infercnv_final_heatmap.png"))
  # removed unneeded files:
  system("rm Output/InferCNV/supervised_clustering/*obj")
  system("rm Output/InferCNV/supervised_clustering/*.dat")
  system("rm -r Output/InferCNV/supervised_clustering/BayesNetOutput*")
  system("rm -r Output/InferCNV/input_files/raw_counts_matrix.txt")
  if (temp_run_HMM) {
    system("rm Output/InferCNV/supervised_clustering/*12*")
    system("rm Output/InferCNV/supervised_clustering/*13*")
  }    
}


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
