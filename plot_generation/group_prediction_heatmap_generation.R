# bash call of script:
#samples=( "normals" )
#annots=( "FALSE" )
#for a in ${annots[@]}; do
#  for s in ${samples[@]}; do
#    mkdir -p /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/logs/$s.annot.$a.group.preds/
#    /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript --vanilla /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/group_prediction_heatmap_generation.R $s $a
#    qsub -wd /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/logs/$s.annot.$a.group.preds/ -pe smp 16 -N $s.annot.$a.group.preds -b y -j y -V -P TumourProgression "${R} CMD BATCH  --no-save '--args $s $a' /share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/scripts/group_prediction_heatmap_generation.R"
#  done;
#done;

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

#custom
#sample_ids <- c("CID4386", "CID43862")
#mets <- c("CID4386", "CID43862")
#names(mets) <- c("Brain met", "Subcutaneous met")
#mets <- "none"
#include_annotations <- TRUE
#subset <- TRUE

subset <- FALSE

temp_args <-
  commandArgs(trailingOnly = T)
temp_sample <- temp_args[1]
print(paste0("temp_sample: ", temp_sample))
include_annotations <- as.logical(temp_args[2])
print(paste0("include annotations: ", include_annotations))

#temp_sample <- "normals"
#include_annotations <- as.logical("FALSE")

if (temp_sample == "normals") {
  sample_ids <- c("CID4520N", "CID4523N", "CID4530N", "IND4", "IND5", "IND6", "IND7")
  mets <- "none"
} else if (temp_sample == "ER") {
  sample_ids <- c("CID3941", "CID3948", "CID4067", "CID4290", "CID4461", "CID4463", 
    "CID4465", "CID4471", "CID4530", "CID4535")
  mets <- "none"
} else if (temp_sample == "HER2") {
  sample_ids <- c("CID3586", "CID4066", "CID4398", 
    "CID3921", "CID3963", "CID45171", "CID45172")
  mets <- "CID45172"
  names(mets) <- "Lymph node met"
} else if (temp_sample == "TNBC") {
  sample_ids <- c("CID4386", "CID43862", "CID43863", "CID44042", 
    "CID44971", "CID44991", "CID44992", "CID4515")
  mets <- c("CID4386", "CID43862", "CID43863", "CID44042", "CID4408", "CID4409", "CID44972")
  names(mets) <- c("Brain met", "Subcutaneous met", "Brain met",
    "Lymph node met", "Lymph node met", "Liver met", "Lymph node met")
} else if (temp_sample == "metaplastic") {
  sample_ids <- c("CID4513", "CID4523")
  mets <- "none"
}

missing_genes_colour <- "white"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
seurat_path <- paste0(project_dir, "results/seurat/")

if (include_annotations) {
  out_path <- paste0(seurat_path, "brca_mini_atlas/", 
    paste(sample_ids, collapse = "_"), "_predictions/")
} else {
  out_path <- paste0(seurat_path, "brca_mini_atlas/", 
    paste(sample_ids, collapse = "_"), "_predictions_no_annotations/")
}
ref_dir <- paste0(project_dir, "/refs/")
system(paste0("mkdir -p ", out_path))


#########################################################################################
### 0. Define functions ###
#########################################################################################

create_group_annotation <- function(df, metadata) {
  m <- match(rownames(df), metadata$cell_ids)
  df$groups <- metadata$cell_type[m]
  df$groups <- gsub("_", " ", df$groups)
  # store groups column in separate data frame and remove from data frame:
  group_annotation_df <- data.frame(group = df$groups, ids = rownames(df),
      stringsAsFactors = F)
  rownames(group_annotation_df) <- group_annotation_df$ids
  group_annotation_df <- data.frame(
    group = group_annotation_df$group, row.names = rownames(group_annotation_df))
  group_annotation_df$group <- as.character(group_annotation_df$group)

  sample_annotation_df <- group_annotation_df
  sample_annotation_df$group <- gsub(" .*$", "", 
  	sample_annotation_df$group)

  site_annotation_df <- sample_annotation_df
  if (mets != "none") {
    for (m in 1:length(mets)) {
      site_annotation_df$group[site_annotation_df$group == mets[m]] <- 
      gsub(
        mets[m], 
        names(mets)[m], 
        site_annotation_df$group[site_annotation_df$group == mets[m]]
      )
    }
    site_annotation_df$group[grep("[m,M]et", site_annotation_df$group, invert = T)] <- "Primary"
  }
  
  # create list of all dfs:
  df_list <- list(group_annotation_df, sample_annotation_df, 
  	site_annotation_df)
  names(df_list) <- c("group_annotation_df", "sample_annotation_df", 
  	"site_annotation_df")

  # determine order of cells:
  all_df <- merge(site_annotation_df, sample_annotation_df, by = "row.names")
  temp_group_df <- group_annotation_df
  temp_group_df$Row.names <- rownames(group_annotation_df)
  all_df <- merge(all_df, temp_group_df, by = "Row.names")

  # order dfs by cell type then sample then site:
  for (i in 1:length(sample_ids)) {
    print(i)
    temp_df <- all_df[
      sample_ids[i] == gsub(" (.*)", "", all_df$group),
    ]
    if (nrow(temp_df) > 1) {
      # order by number:
      temp_df <- temp_df[
        order(
          as.numeric(
            as.character(
              gsub("^.* ", "", temp_df$group)
            )
          )
        ), 
      ]
      # order by cell type:
      temp_order <- unlist(
        lapply(strsplit(temp_df$group, " "), function(x) {
          return(x[2])
        })
      )
      temp_df <- temp_df[order(temp_order),]
  
      if (i==1) {
        result_df <- temp_df
      } else {
        result_df <- rbind(result_df, temp_df)
      }
    }
  }
    
  # order by site:
  all_df <- rbind(result_df[grep("met", all_df$group.x, invert = T),],
    result_df[grep("met", all_df$group.x),])
  cell_order <- all_df$Row.names

  df_list <- lapply(df_list, function(x) {
      res <- data.frame(row.names = cell_order, x[cell_order,])
      colnames(res) <- "group"
      return(res)
  })

  # check dfs are in correct gene order:
  for (r in df_list) {
    if (class(r) == "data.frame") {
      print(paste0("Is df in correct cell order? ", identical(as.character(rownames(r)), 
        as.character(cell_order))))
    }
  }

  # define group heatmap colours:
  col_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
    "#660200", "#918940", "black")
  col_palette <- c(col_palette, col_palette)
  cluster_number <- length(unique(df_list$group_annotation_df$group))
  cluster_cols <- col_palette[1:cluster_number]
  names(cluster_cols) <- unique(df_list$group_annotation_df$group)
  
  # create heatmap of group annot:
  group_annotation <- Heatmap(df_list$group_annotation_df, col = cluster_cols, 
    name = "group_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Cluster", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

   # define group heatmap colours:
   col_palette2 <- c("#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
  	"#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
    "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC")
  sample_number <- length(unique(df_list$sample_annotation_df$group))
  sample_cols <- col_palette2[1:sample_number]
  names(sample_cols) <- unique(df_list$sample_annotation_df$group)

  # create heatmap of sample annot:
  sample_annotation <- Heatmap(df_list$sample_annotation_df, col = sample_cols, 
    name = "sample_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Sample", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  # define group heatmap colours:
  col_palette3 <- c("#b2182b", 
            "#85929E", "#9B59B6", "#74add1",
            "#1b7837", "#b8e186", "#fed976",
            "#e7298a", "#18ffff", "#ef6c00",
            "#A93226", "black","orange",
            "#b8bc53", "#5628ce", "#fa909c",
            "#8ff331","#270e26")
  site_number <- length(unique(df_list$site_annotation_df$group))
  site_cols <- col_palette3[1:site_number]
    names(site_cols) <- unique(df_list$site_annotation_df$group)

  # create heatmap of sample annot:
  site_annotation <- Heatmap(df_list$site_annotation_df, col = site_cols, 
    name = "site_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Sample", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  result_list <- list(df_list$group_annotation_df, group_annotation,
    df_list$sample_annotation_df, sample_annotation, df_list$site_annotation_df, site_annotation)
  names(result_list) <- c("group_annotation_df", "group_annotation",
    "sample_annotation_df", "sample_annotation", "site_annotation_df", "site_annotation")

  return(result_list)
}

# fetch co-ordinates of chromosomal boundaries:
fetch_chromosome_boundaries <- function(df, ref_dir) {
  # load in gene annotations:
  gene_order <- read.table(paste0(ref_dir, "infercnv_gene_order.txt"), 
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

# function to call ggplot default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# function to create GIN level annotation:
create_GIN_annotation <- function(df) {

  # adjust HMM values
  # 0 = complete loss, change to 2
  # 2 = loss of one copy, change to 1
  # 3 = neutral, change to 0
  # 4 = addition of one copy, change to 1
  # 5 = addition of two copies, change to 2
  # 6 = addition of three or more copies, change to 3
  old_scores <- c(0, 2, 3, 4, 5, 6)
  new_scores <- c(2, 1, 0, 1, 2, 3)
  for (s in 1:length(new_scores)) {
    df[df == old_scores[s]] <- new_scores[s]
  }

  # sum total CNV levels for each cell:
  GIN_levels <- apply(df, 1, function(x) {
  	x[is.na(x)] <- 0
  	return(sum(x))
  })
  
  # create GIN annotation:
  GIN_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      GIN_levels,
      gp = gpar(
        col = "#AF548E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )
  result_list <- list(GIN_annotation, GIN_levels)
  names(result_list) <- c("GIN_annotation", "GIN_levels")

  return(result_list)
}

# function to create nUMI and nGene barplot annotations:
create_QC_annotation <- function(seurat_objects_list, df) {
  
  # load each seurat object and create QC metrics df:
  for (j in 1:length(seurat_objects_list)) {
    seurat_10X <- seurat_objects_list[[j]]
    qc_df <- data.frame(
  	seurat_10X@meta.data$nCount_RNA,
  	seurat_10X@meta.data$nFeature_RNA,
  	row.names = as.character(names(Idents(seurat_10X)))
    )
    colnames(qc_df) <-  c("nUMI", "nGene")

    if (j==1) {
      all_qc_df <- qc_df
    } else {
      all_qc_df <- rbind(all_qc_df, qc_df)
    }
  }

  m <- match(
    rownames(heatmap_df), rownames(all_qc_df)
  )
  all_qc_df <- all_qc_df[m,]

  nUMI_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      all_qc_df$nUMI,
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
      all_qc_df$nGene,
      gp = gpar(
        col = "#9ECAE1", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )

  result_list <- list(nUMI_annotation, nGene_annotation, all_qc_df)
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
    expand_no <- as.integer(ncol(heatmap_df)/1300)
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


#########################################################################################
### 1. Load data ###
#########################################################################################

# fetch vector of all genes:
all_genes <- as.character(read.table(paste0(ref_dir, "/infercnv_gene_order.txt"))$V1)

# load samples, metadata and fetch vector of genes contained in samples:
for (i in 1:length(sample_ids)) {
  print(i)
  metadata_dir <- paste0(seurat_path, "brca_mini_atlas/seurat_", sample_ids[i], 
    "/Output/InferCNV/input_files/")
  infercnv_metadata <- read.table(paste0(metadata_dir, "metadata.txt"),
    header=F)
  colnames(infercnv_metadata) <- c("cell_ids", "cell_type")
  # remove CAFs from heatmap df and metadata:
  cells_to_remove <- 
  infercnv_metadata$cell_ids[grep("pithelial", infercnv_metadata$cell_type, invert=T)]
  infercnv_metadata <- infercnv_metadata[!(infercnv_metadata$cell_ids %in% cells_to_remove),]
  print(dim(infercnv_metadata))

  if ( length(grep(sample_ids[i], infercnv_metadata$cell_type)) < 1) {
    infercnv_metadata$cell_type <- paste0(sample_ids[i], "_", infercnv_metadata$cell_type)
  }

  infercnv_output_filename <- list.files(
    paste0(seurat_path, "brca_mini_atlas/seurat_", 
    sample_ids[i], 
    "/Output/InferCNV/supervised_clustering"), 
    pattern = "infercnv.14_HMM_predHMMi6.hmm_mode-samples.repr_intensities.observations.txt",
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

###
#save.image(paste0(out_path, "/temp_image.RData"))
###

# remove 'luminal' from epithelial cells as Garnett labels everything luminal:
metadata_df$cell_type <- gsub("Luminal_", "", metadata_df$cell_type)
# create group annotation df
group_annotation <- create_group_annotation(heatmap_df, metadata_df)
# ensure heatmap_df has same cell order as group_annotation$group_annotation_df:
m <- match(rownames(group_annotation$group_annotation_df), rownames(heatmap_df))
heatmap_df <- heatmap_df[m,]

###
#save.image(paste0(out_path, "temp_image_2_subset.RData"))
###

# define chromosome lengths in terms of number of filtered genes
chr_data <- fetch_chromosome_boundaries(heatmap_df, ref_dir)

# add GIN annotation:
GIN_annotation <- create_GIN_annotation(heatmap_df)

# add nUMI and nGene barplot annotations:
for (i in 1:length(sample_ids)) {
  print(paste0("Loading ", sample_ids[i], " seurat object..."))
  if (i==1) {
    seurat_list <- list(readRDS(paste0(seurat_path, "brca_mini_atlas/", "seurat_", 
    sample_ids[i], "/Output/Rdata/03_seurat_object_processed.RData")))
  } else {
    seurat_list[[i]] <- readRDS(paste0(seurat_path, "brca_mini_atlas/", "seurat_", 
    sample_ids[i], "/Output/Rdata/03_seurat_object_processed.RData"))
  }
}
QC_annotation <- create_QC_annotation(seurat_list, heatmap_df)

if (include_annotations) {

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

  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = colorRamp2(c(0, 2, 3, 4, 5, 6), 
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

  if (mets != "none") {

    ht_list <- group_annotation$site_annotation + group_annotation$sample_annotation + 
    group_annotation$group_annotation + final_heatmap + GIN_annotation$GIN_annotation + 
    QC_annotation$nUMI_annotation + QC_annotation$nGene_annotation
    annotated_heatmap <- grid.grabExpr(
    	draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
    )
  
    # determine where starting co-ordinates for heatmap are based upon longest cluster name
    # (0.00604 units per character):
    longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
    x_coord <- longest_cluster_name*0.00455
    if (missing_genes_colour != "white") {
      pdf(paste0(out_path, "final_infercnv_heatmap_missing_genes_coloured.pdf"), 
          height = 14.5, width = 17)
    } else {
        pdf(paste0(out_path, "final_infercnv_heatmap.pdf"), 
          height = 14.5, width = 17)
    }
  
      grid.newpage()
      # plot Normal subtype:
      pushViewport(viewport(x = x_coord, y = 0.090,
                            width = 0.834+0.012, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[5]])
      popViewport()
      
      # plot Basal subtype:
      pushViewport(viewport(x = x_coord+0.008, y = 0.169,
                            width = 0.834+0.004, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[4]])
      popViewport()
      
      # plot Her2 subtype:
      pushViewport(viewport(x = x_coord+0.012, y = 0.251, 
                            width = 0.834, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[3]])
      popViewport()
      
      # plot LumB subtype:
      pushViewport(viewport(x = x_coord+0.007, y = 0.333, 
                            width = 0.834+0.005, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[2]])
      popViewport()
      
      # plot LumA subtype:
      pushViewport(viewport(x = x_coord+0.007, y = 0.411, 
                            width = 0.834+0.005, height = 0.08, just = c("left", "top")))
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
  
      pushViewport(viewport(x=x_coord + 0.83, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("GIN", rot=55)
      popViewport()
  
      pushViewport(viewport(x=x_coord + 0.855, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("nUMI", rot=55)
      popViewport()
  
      pushViewport(viewport(x=x_coord + 0.865, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("nGene", rot=55)
      popViewport()
  
    dev.off()

  } else {

    ht_list <- group_annotation$sample_annotation + 
      group_annotation$group_annotation + final_heatmap + GIN_annotation$GIN_annotation + 
      QC_annotation$nUMI_annotation + QC_annotation$nGene_annotation
    
    # determine where starting co-ordinates for heatmap are based upon longest cluster name
    # (0.0035 units per character):
    longest_cluster_name <- max(nchar(unique(as.character(group_annotation$group_annotation_df$group))))
    x_coord <- longest_cluster_name*0.0035

    annotated_heatmap <- grid.grabExpr(
    	draw(ht_list, gap = unit(3, "mm"), heatmap_legend_side = "left")
    )
  
    if (missing_genes_colour != "white") {
      pdf(paste0(out_path, "final_infercnv_heatmap_missing_genes_coloured.pdf"), 
          height = 14.5, width = 17)
    } else {
        pdf(paste0(out_path, "final_infercnv_heatmap.pdf"), 
          height = 14.5, width = 17)
    }
  
      grid.newpage()
      # plot Normal subtype:
      pushViewport(viewport(x = x_coord+0.003, y = 0.090,
                            width = 0.82+0.019, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[5]])
      popViewport()
      
      # plot Basal subtype:
      pushViewport(viewport(x = x_coord+0.008, y = 0.169,
                            width = 0.82+0.014, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[4]])
      popViewport()
      
      # plot Her2 subtype:
      pushViewport(viewport(x = x_coord+0.01, y = 0.251, 
                            width = 0.82+0.012, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[3]])
      popViewport()
      
      # plot LumB subtype:
      pushViewport(viewport(x = x_coord+0.007, y = 0.333, 
                            width = 0.82+0.015, height = 0.08, just = c("left", "top")))
      grid.draw(metabric_plots[[2]])
      popViewport()
      
      # plot LumA subtype:
      pushViewport(viewport(x = x_coord+0.007, y = 0.411, 
                            width = 0.82+0.015, height = 0.08, just = c("left", "top")))
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
  
      pushViewport(viewport(x=x_coord + 0.85, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("GIN", rot=55)
      popViewport()
  
      pushViewport(viewport(x=x_coord + 0.87, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("nUMI", rot=55)
      popViewport()
  
      pushViewport(viewport(x=x_coord + 0.885, y=0.405, width = 0.1, height = 0.1, just = "bottom"))
      #grid.draw(lollipop)
      grid.text("nGene", rot=55)
      popViewport()
  
    dev.off()
  }

} else {
  
  # replace column labels with two letter label for positioning:
  plot_object <- heatmap_df
  colnames(plot_object) <- rep("la", ncol(plot_object))
  
  print("Generating final heatmap...")
  
  na_less_vector <- unlist(plot_object)
  na_less_vector <- na_less_vector[!is.na(na_less_vector)]

  # create main CNV heatmap:
  final_heatmap <- Heatmap(
    plot_object, name = paste0("hm"), 
    col = colorRamp2(c(0, 2, 3, 4, 5, 6), 
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

  ht_list <- group_annotation$sample_annotation + group_annotation$group_annotation + 
    final_heatmap + GIN_annotation$GIN_annotation + 
      QC_annotation$nUMI_annotation + QC_annotation$nGene_annotation
  annotated_heatmap <- grid.grabExpr(
    draw(ht_list, heatmap_legend_side = "left")
  )
  if (missing_genes_colour != "white") {
    pdf(paste0(out_path, "final_infercnv_prediction_heatmap_missing_genes_coloured_non_annotated.pdf"), 
        height = 9, width = 14)
  } else {
      pdf(paste0(out_path, "final_infercnv_prediction_heatmap_non_annotated.pdf"), 
        height = 9, width = 14)
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
  
  dev.off()

}

# convert pdf to png:
system(paste0("for p in ", out_path, "*.pdf; do echo $p; f=$(basename $p); echo $f; ",
  "new=$(echo $f | sed 's/.pdf/.png/'); echo $new; ", 
  "convert -density 150 ", out_path, "$f -quality 90 ", out_path, "$new; done"))

print(paste0("Group heatmap created, output in ", out_path))

