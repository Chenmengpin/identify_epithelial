prepare_infercnv_metadata <- function(seurat_object, temp_resolution, 
  temp_temp_half_stromal) {
  
  # fetch correct clusters:
  if (temp_resolution != "none") {
    Idents(seurat_object) <- 
      eval(parse(text = paste0("seurat_object@meta.data$", temp_resolution)))
    garnett_calls <- eval(
      parse(
        text = paste0("seurat_object@meta.data$garnett_seurat_cluster_call_major_",
          temp_resolution)
      )
    )

    annotated_idents <- gsub(
      " ", "_", paste0(garnett_calls, "_", Idents(seurat_object))
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
    annotated_ident_cluster_numbers <- gsub("^.*CAF[1,2]_", "", 
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
  }

  # create infercnv metadata
  temp_metadata <- data.frame(cell_ids = names(Idents(seurat_object)), 
    cell_type = Idents(seurat_object), stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]

  if (temp_half_stromal) {
    stromal_ids <- temp_metadata$cell_ids[
      grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
        temp_metadata$cell_type, invert=T)
    ]
    control_stromal <- sample(stromal_ids, length(stromal_ids)/2, 
      length(stromal_ids))
    cell_type <- as.character(temp_metadata$cell_type)
    cell_type[temp_metadata$cell_ids %in% control_stromal] <-
      "Control_Stromal"
    temp_metadata$cell_type <- factor(cell_type)
  }

  # record number per cell type
  number_per_cell_type <- as.data.frame(table(temp_metadata$cell_type))

  return(list(temp_metadata, number_per_cell_type, seurat_object))
}
