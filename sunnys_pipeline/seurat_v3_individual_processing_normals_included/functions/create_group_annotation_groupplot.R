create_group_annotation_groupplot <- function(df, metadata) {
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
  sample_annotation_df$group <- gsub("_.*$", "", 
    rownames(sample_annotation_df))

  type_annotation_df <- sample_annotation_df

  type_annotation_df$group[grep("PDX", type_annotation_df$group)] <- "PDX"
  type_annotation_df$group[grep("CID", type_annotation_df$group)] <- "Metastasis"
  type_annotation_df$group[grep("[0-9]N", type_annotation_df$group)] <- "Normal"

  # create list of all dfs:
  df_list <- list(group_annotation_df, sample_annotation_df, 
    type_annotation_df)
  names(df_list) <- c("group_annotation_df", "sample_annotation_df", 
    "type_annotation_df")

  # determine order of cells:
  all_df <- merge(type_annotation_df, sample_annotation_df, by = "row.names")
  temp_group_df <- group_annotation_df
  temp_group_df$Row.names <- rownames(group_annotation_df)
  all_df <- merge(all_df, temp_group_df, by = "Row.names")

  # order dfs by cell type then sample then type:
  for (i in 1:length(sample_ids)) {
    print(i)
    temp_df <- all_df[
      sample_ids[i] == gsub("_.*$", "", all_df$Row.names),
    ]
    if (nrow(temp_df) > 1) {
      # order by number:
      temp_df <- temp_df[
        order(
          as.numeric(
            as.character(
              gsub("^.* ", "", temp_df$group)
            )
          )
        ), 
      ]
      # order by cell type:
      temp_order <- unlist(
        lapply(strsplit(temp_df$group, " "), function(x) {
          return(x[2])
        })
      )
      temp_df <- temp_df[order(temp_order),]
  
      if (i==1) {
        result_df <- temp_df
      } else {
        result_df <- rbind(result_df, temp_df)
      }
    }
  }
    
  # order by type:
  all_df <- rbind(result_df[grep("Normal", result_df$group.x),],
    result_df[grep("Metastasis", result_df$group.x),])
  all_df <- rbind(all_df,result_df[grep("PDX", result_df$group.x),])
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
    labels_gp = gpar(fontsize = 8)),
    show_heatmap_legend =F)

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

  # define group heatmap colours:
  col_palette3 <- c("#b2182b", 
            "#85929E", "#9B59B6", "#74add1",
            "#1b7837", "#b8e186", "#fed976",
            "#e7298a", "#18ffff", "#ef6c00",
            "#A93226", "black","orange",
            "#b8bc53", "#5628ce", "#fa909c",
            "#8ff331","#270e26")
  type_number <- length(unique(df_list$type_annotation_df$group))
  type_cols <- col_palette3[1:type_number]
    names(type_cols) <- unique(df_list$type_annotation_df$group)

  # create heatmap of sample annot:
  type_annotation <- Heatmap(df_list$type_annotation_df, col = type_cols, 
    name = "type_annotation", 
    width = unit(4, "mm"), 
    show_row_names = F, show_column_names = F,
    heatmap_legend_param = list(title = "Sample", 
    title_gp = gpar(fontsize = 8, fontface = "bold"), 
    labels_gp = gpar(fontsize = 8)))

  result_list <- list(df_list$group_annotation_df, group_annotation,
    df_list$sample_annotation_df, sample_annotation, df_list$type_annotation_df, type_annotation)
  names(result_list) <- c("group_annotation_df", "group_annotation",
    "sample_annotation_df", "sample_annotation", "type_annotation_df", "type_annotation")

  return(result_list)
}
