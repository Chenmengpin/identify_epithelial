create_QC_annotation <- function(qc_df, df) {
  
  # add zero QC values for SNP array CNVs
  array_rownames <- rownames(df)[grep("array", rownames(df))]
  if (length(array_rownames) > 0) {
    array_qc <- data.frame(
      row.names = array_rownames, 
      nUMI = rep(0, length(array_rownames)),
      nGene = rep(0, length(array_rownames))
    )
    qc_df <- rbind(array_qc, qc_df)
  }
  
  # order the two dfs:
  m <- match(
    rownames(df), rownames(qc_df)
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
