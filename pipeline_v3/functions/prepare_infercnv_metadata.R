prepare_infercnv_metadata <- function(seurat_object, combined_sample,
  temp_resolution, subset_data = FALSE, count_df, throw_unneeded, pool_cells) {
  
  # fetch PC A, res 0.8 clusters:
  if (!combined_sample) {
    Idents(seurat_object) <- 
      eval(parse(text = paste0("seurat_object@meta.data$", temp_resolution)))
    annotated_idents <- gsub(
      " ", "_", paste0(seurat_object@meta.data$garnett_call_ext_major, " ", 
      Idents(seurat_object))
    )
    Idents(seurat_object) <- factor(annotated_idents)
#    # only keep annotation with largest no cells per existing cluster
#    annotated_ident_numbers <- table(annotated_idents)
#    annotated_ident_cluster_numbers <- gsub("^.*_CID", "CID", 
#      names(annotated_ident_numbers))
#    annotated_ident_cluster_numbers <- gsub("^.*_IND", "IND", 
#      annotated_ident_cluster_numbers)
#    annotated_ident_cluster_numbers <- gsub("^.*_PDX", "PDX", 
#      annotated_ident_cluster_numbers)
#    annotated_ident_cluster_numbers <- gsub("^.*[a-z]_", "", 
#      annotated_ident_cluster_numbers)
#    annotated_ident_cluster_numbers <- gsub("^.*CAF.*_", "", 
#      annotated_ident_cluster_numbers)
#  
#    for ( c in levels(Idents(seurat_object)) ) {
#      print(c)
#      annotated_ident_index <- which(annotated_ident_cluster_numbers == c)
#      potential_clusters <- annotated_ident_numbers[annotated_ident_index]
#      correct_cluster <- 
#        names(potential_clusters[potential_clusters == 
#          max(potential_clusters)])
#        # label cells in seurat object as the correct cluster
#        levels(Idents(seurat_object))[levels(Idents(seurat_object)) == c] <- 
#          correct_cluster
#    }
  }

  if (pool_cells) {
    temp_idents <- as.character(Idents(seurat_object))
    temp_idents[grep("pithelial", temp_idents)] <- "Epithelial"
    temp_idents[grep("nknown", temp_idents)] <- "Unknown"
    temp_idents[grep("CAF", temp_idents)] <- "CAFs"
    temp_idents[grep("pithelial|nknown|CAF", temp_idents, invert=T)] <- "Stromal"
    Idents(seurat_object) <- factor(temp_idents)
  }

  # create infercnv metadata
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
    cell_type = Idents(seurat_object), stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]

  if (throw_unneeded) {
    temp_metadata <- temp_metadata[grep("[u,U]nknown|[u,U]nassigned|CAFs", 
      temp_metadata$cell_type, invert=T),]
  }
  
  if (subset_data) {
    temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
  }

  # record number per cell type
  number_per_cell_type <- as.data.frame(table(temp_metadata$cell_type))

  return(list(temp_metadata, number_per_cell_type, seurat_object))
}
