create_GIN_annotation <- function(df, type_annotation_df) {

    # adjust HMM values
  # 0 = complete loss, change to 2
  # 0.5 = loss of one copy, change to 1
  # 1 = neutral, change to 0
  # 1.5 = addition of one copy, change to 1
  # 2 = addition of two copies, change to 2
  # 3 = addition of three or more copies, change to 3
  old_scores <- c(0, 0.5, 1, 1.5, 2, 3)
  new_scores <- c(2, 1, 0, 1, 2, 3)
  for (s in 1:length(new_scores)) {
    df[df == old_scores[s]] <- new_scores[s]
  }

  # sum total CNV levels for each cell:
  GIN_levels <- apply(df, 1, function(x) {
    x[is.na(x)] <- 0
    return(sum(x))
  })

  # normalise to number of genes in plot:
  GIN_levels <- GIN_levels/ncol(df)
  
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

  cancer_GIN <- GIN_levels[type_annotation_df$group == "Metastasis" | 
    type_annotation_df$group == "PDX"]
  min_cancer_GIN <- paste0(round(min(cancer_GIN), 3), " CNV/gene")
  max_cancer_GIN <- paste0(round(max(cancer_GIN), 3), " CNV/gene")

  normal_GIN <- GIN_levels[type_annotation_df$group == "Normal"]
  min_normal_GIN <- paste0(round(min(normal_GIN), 3), " CNV/gene")
  max_normal_GIN <- paste0(round(max(normal_GIN), 3), " CNV/gene")

  result_list <- list(GIN_annotation, GIN_levels, min_normal_GIN,
    max_normal_GIN, min_cancer_GIN, max_cancer_GIN)
  names(result_list) <- c("GIN_annotation", "GIN_levels", "min_normal_GIN",
    "max_normal_GIN", "min_cancer_GIN", "max_cancer_GIN")

  return(result_list)
}
