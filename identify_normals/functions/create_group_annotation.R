create_group_annotation <- function(df, metadata) {
  m <- match(rownames(df), metadata$cell_ids)
  m <- m[!(is.na(m))]
  df$groups <- metadata$cell_type[m]
  df$groups <- gsub("_", " ", df$groups)
  # store groups column in separate data frame and remove from data frame:
  group_annotation_df <- data.frame(group = df$groups, ids = rownames(df),
      stringsAsFactors = F)
  rownames(group_annotation_df) <- group_annotation_df$ids
  group_annotation_df <- data.frame(
    group = group_annotation_df$group, row.names = rownames(group_annotation_df))
  group_annotation_df$group <- as.character(group_annotation_df$group)

  # define group heatmap colours:
  extra_colours <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  col_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
    "#660200", "#918940", "black", "#9ECAE1", "#3F007D", "#67000D", "#FDD0A2", "#08519C",
    "#DBB335", "#5AA050", "#807DBA", "#1A6000", "#F16913", "#FD8D3C", "#DEEBF7", 
    "#7F2704", "#DADAEB", "#FC9272", "#BCBDDC", extra_colours, "#b2182b", "#85929E", 
    "#9B59B6", "#74add1","#1b7837", "#b8e186", "#fed976","#e7298a", "#18ffff", "#ef6c00",
    "#A93226", "black","orange", "#b8bc53", "#5628ce", "#fa909c", "#8ff331","#270e26")
  cluster_number <- length(unique(group_annotation_df$group))
  cluster_cols <- col_palette[1:cluster_number]
  names(cluster_cols) <- unique(group_annotation_df$group)
  
  # create heatmap of group annot:
  group_annotation <- Heatmap(group_annotation_df, col = cluster_cols, 
    name = "group_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Cluster", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  result_list <- list(group_annotation_df, group_annotation)
  names(result_list) <- c("df", "annotation")

  return(result_list)
}
