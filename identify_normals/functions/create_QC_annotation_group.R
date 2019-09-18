create_QC_annotation_group <- function(df, seurat_object) {
  
  # create QC metrics df:
  seurat_10X <- seurat_object[[p]]
  qc_df <- data.frame(
	seurat_10X@meta.data$nCount_RNA,
	seurat_10X@meta.data$nFeature_RNA,
	row.names = as.character(names(Idents(seurat_10X)))
  )
  colnames(qc_df) <-  c("nUMI", "nGene")

   m <- match(
    rownames(df), rownames(qc_df)
  )
  qc_df <- qc_df[m,]

  p <<- p+1

  return(qc_df)

}
