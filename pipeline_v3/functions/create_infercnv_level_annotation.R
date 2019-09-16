create_infercnv_level_annotation <- function(df) {

  library(scales)
  # sum total CNV levels for each cell:
  infercnv_level <- apply(df, 1, function(x) {
    x[is.na(x)] <- 0
    return(sum(x))
  })

  # normalise to number of genes in plot:
  infercnv_level <- round(rescale(infercnv_level/ncol(df), c(1,100)), 0)
  
  # create infercnv_level annotation:
  infercnv_annotation <- rowAnnotation(
    correlation_annotation = anno_barplot(
      infercnv_level,
      gp = gpar(
        col = "#AF548E", 
        width = unit(4, "cm")
      ), 
      border = FALSE, 
      which = "row", 
      axis = F
    )
  )

  # create infercnv_level data frame:
  infercnv_level <- data.frame(cell_id=rownames(df), infercnv_level=infercnv_level)

  result_list <- list(infercnv_annotation, infercnv_level)
  names(result_list) <- c("infercnv_annotation", "infercnv_level")

  return(result_list)
}
