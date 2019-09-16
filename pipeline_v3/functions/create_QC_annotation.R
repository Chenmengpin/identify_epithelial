create_QC_annotation <- function(seurat_object, df) {
  
  # create QC metrics df:
  seurat_10X <- seurat_object
  qc_df <- data.frame(
	seurat_10X@meta.data$nCount_RNA,
	seurat_10X@meta.data$nFeature_RNA,
	row.names = as.character(names(Idents(seurat_10X)))
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
