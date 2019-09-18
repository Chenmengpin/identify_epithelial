create_normal_call_annotation <- function(df, metadata, GIN_levels, 
  cor_thresh, GIN_thresh) {
  # determine top 5% of tumour cells:
  epithelial_cells <- rownames(df)
  epithelial_GIN <- GIN_levels[GIN_levels$cell_id %in% epithelial_cells,]
  epithelial_GIN <- epithelial_GIN[order(epithelial_GIN$GIN, decreasing=T),]
  top_GIN <- head(epithelial_GIN, nrow(epithelial_GIN)*0.05)
  # find average genome-wide CNV predictions across genome:
  top_GIN_CNV_average <- apply(df[top_GIN$cell_id,], 2, mean)
  # find correlations of each cell's CNVs with top_GIN_CNV_average:
  GIN_correlations <- apply(df, 1, function(x) {
    if (length(unique(as.numeric(x))) == 1) {
      cor_result <- data.frame(cor.estimate="no_CNVs_recorded", 
        cor.p.value="no_CNVs_recorded")
    } else {
      cor <- cor.test(as.numeric(x), top_GIN_CNV_average, method = "kendall")
      cor_result <- data.frame(cor$estimate, cor$p.value)
    }
    return(cor_result)
  })
  cor_df <- do.call("rbind", GIN_correlations)
  combined_df <- cbind(cor_df, GIN_levels$GIN)
  sig_cors <- combined_df[combined_df$cor.p.value < 0.05,]
  sig_cors <- data.frame(cell_id=rownames(sig_cors), cor=sig_cors$cor.estimate)
  normal_cors <- sig_cors$cell_id[sig_cors$cor < 0.3]
  normal_GIN <- epithelial_GIN$cell_id[epithelial_GIN$GIN < 0.2]
  normal_calls <- normal_cors[normal_cors %in% normal_GIN]
  normal_call_df <- data.frame(row.names = epithelial_cells, 
    call = rep("other", length(epithelial_cells)), stringsAsFactors=F)
  normal_call_df$call[rownames(normal_call_df) %in% normal_calls] <- "normal"
  normal_call_annotation <- Heatmap(normal_call_df, 
    col = c("other" = "#E7E4D3", "normal" = "#1B7837"), 
    name = "normal_call_annotation", width = unit(6, "mm"), 
    show_row_names = F, show_column_names = F, 
    heatmap_legend_param = list(title = "Normal epithelial calls", title_gp = gpar(fontsize = 8, 
    fontface = "bold"), labels_gp = gpar(fontsize = 6)))
  res <- list(normal_call_df, normal_call_annotation, cor_df)
  names(res) <- c("normal_call_df", "normal_call_annotation", "correlation_df")
  return(res)
}