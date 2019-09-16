prepare_infercnv_metadata <- function(seurat_object, temp_resolution, subset_data = FALSE, count_df) {
  
  cell_ids <- names(seurat_object@ident)
  if (temp_resolution != "none") {
    seurat_object@ident <- factor(paste0(gsub("_.*$", "", names(seurat_object@ident)), "_", 
    eval(parse(text=paste0("seurat_object@meta.data$", temp_resolution)))))
    names(seurat_object@ident) <- cell_ids
  }

  # create infercnv metadata
  temp_metadata <- data.frame(cell_ids = names(seurat_object@ident), 
    cell_type = seurat_object@ident, stringsAsFactors = F)
  temp_metadata <- temp_metadata[order(temp_metadata$cell_type),]
  m <- match(colnames(count_df), temp_metadata$cell_id)
  temp_metadata <- temp_metadata[m,]

  if (subset_data) {
    temp_metadata <- temp_metadata[temp_metadata$cell_ids %in% colnames(count_df),]
  }

  # record number per cell type
  number_per_cell_type <- as.data.frame(table(temp_metadata$cell_type))

  return(list(temp_metadata, number_per_cell_type, seurat_object))
}
