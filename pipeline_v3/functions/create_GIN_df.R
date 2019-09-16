create_GIN_df <- function(df) {

  if (is.null(df)) {
    return(NULL)
  } else {
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
    GIN_levels <- round(GIN_levels/ncol(df), 3)
  
    # create GIN levels data frame:
    GIN_levels <- data.frame(cell_id=rownames(df), GIN=GIN_levels)
    
    return(GIN_levels)
  }
}
