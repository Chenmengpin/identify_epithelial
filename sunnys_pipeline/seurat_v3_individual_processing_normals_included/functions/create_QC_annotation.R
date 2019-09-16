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
