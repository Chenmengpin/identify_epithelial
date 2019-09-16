create_group_annotation <- function(df, metadata) {
  m <- match(rownames(df), metadata$cell_ids)
  df$groups <- metadata$cell_type[m]
  df$groups <- gsub("_", " ", df$groups)
  # store groups column in separate data frame and remove from data frame:
  group_annotation_df <- data.frame(group = df$groups, ids = rownames(df),
      stringsAsFactors = F)
  rownames(group_annotation_df) <- group_annotation_df$ids
  group_annotation_df <- data.frame(
    group = group_annotation_df$group, row.names = rownames(group_annotation_df))
  group_annotation_df$group <- as.character(group_annotation_df$group)

  sample_annotation_df <- group_annotation_df
  sample_annotation_df$group <- gsub(" .*$", "", sample_annotation_df$group)

  # create list of all dfs:
  df_list <- list(group_annotation_df, sample_annotation_df)
  names(df_list) <- c("group_annotation_df", "sample_annotation_df")

  # determine order of cells:
  all_df <- merge(sample_annotation_df, group_annotation_df, by = "row.names")

  # order dfs by cell type:
  all_df <- all_df[order(all_df$group.y),]
  cell_order <- all_df$Row.names

  df_list <- lapply(df_list, function(x) {
      res <- data.frame(row.names = cell_order, x[cell_order,])
      colnames(res) <- "group"
      return(res)
  })

  # check dfs are in correct gene order:
  for (r in df_list) {
    if (class(r) == "data.frame") {
      print(paste0("Is df in correct cell order? ", identical(as.character(rownames(r)), 
        as.character(cell_order))))
    }
  }

  # define group heatmap colours:
  col_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
    "#660200", "#918940", "black")
  col_palette <- c(col_palette, col_palette)
  cluster_number <- length(unique(df_list$group_annotation_df$group))
  cluster_cols <- col_palette[1:cluster_number]
  names(cluster_cols) <- unique(df_list$group_annotation_df$group)
  
  # create heatmap of group annot:
  group_annotation <- Heatmap(df_list$group_annotation_df, col = cluster_cols, 
    name = "group_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Cluster", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  # define group heatmap colours:
  col_palette2 <- c("#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
    "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
    "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC")
  sample_number <- length(unique(df_list$sample_annotation_df$group))
  sample_cols <- col_palette2[1:sample_number]
  names(sample_cols) <- unique(df_list$sample_annotation_df$group)

  # create heatmap of sample annot:
  sample_annotation <- Heatmap(df_list$sample_annotation_df, col = sample_cols, 
    name = "sample_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Sample", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  result_list <- list(df_list$group_annotation_df, group_annotation,
    df_list$sample_annotation_df, sample_annotation)
  names(result_list) <- c("group_annotation_df", "group_annotation",
    "sample_annotation_df", "sample_annotation")

  return(result_list)
}
