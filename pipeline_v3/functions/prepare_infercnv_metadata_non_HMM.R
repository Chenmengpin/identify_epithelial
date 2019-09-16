prepare_infercnv_metadata <- function(seurat_object, combined_sample,
  temp_resolution, subset_data = FALSE, count_df, supervised_clustering) {
  
  # fetch clusters of correct resolution:
  Idents(seurat_object) <- 
    eval(parse(text = paste0("seurat_object@meta.data$", temp_resolution)))
  annotated_idents <- gsub(
    " ", "_", paste0(seurat_object@meta.data$garnett_call_ext_major, " ", 
    Idents(seurat_object))
  )
  # only keep annotation with largest no cells per existing cluster
  annotated_ident_numbers <- table(annotated_idents)
  annotated_ident_cluster_numbers <- gsub("^.*_CID", "CID", 
    names(annotated_ident_numbers))
  annotated_ident_cluster_numbers <- gsub("^.*_IND", "IND", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*_PDX", "PDX", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*[a-z]_", "", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*CAF.*_", "", 
    annotated_ident_cluster_numbers)

  for ( c in levels(Idents(seurat_object)) ) {
    print(c)
    annotated_ident_index <- which(annotated_ident_cluster_numbers == c)
    potential_clusters <- annotated_ident_numbers[annotated_ident_index]
    correct_cluster <- 
      names(potential_clusters[potential_clusters == 
        max(potential_clusters)])
      # label cells in seurat object as the correct cluster
      levels(Idents(seurat_object))[levels(Idents(seurat_object)) == c] <- 
        correct_cluster
  }


  # create infercnv metadata
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
    cell_type = Idents(seurat_object), stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]

  if (subset_data) {
    temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
  }

  # record number per cell type
  number_per_cell_type <- as.data.frame(table(temp_metadata$cell_type))

  return(list(temp_metadata, number_per_cell_type, seurat_object))
}
