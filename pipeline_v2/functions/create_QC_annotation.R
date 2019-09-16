create_QC_annotation <- function(seurat_object, df) {
  
  # load each seurat object and create QC metrics df:
  qc_df <- data.frame(
	seurat_object@meta.data$nUMI,
	seurat_object@meta.data$nGene,
	row.names = as.character(names(seurat_object@ident))
  )
  colnames(qc_df) <-  c("nUMI", "nGene")

  m <- match(
    rownames(heatmap_df), rownames(qc_df)
  )
  qc_df <- qc_df[m,]

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
