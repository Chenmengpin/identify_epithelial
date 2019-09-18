function (infercnv_obj, p_val, hclust_method, partition_method,
    restrict_to_DE_genes = FALSE)
{
    flog.info(sprintf("define_signif_tumor_subclusters(p_val=%g",
        p_val))
    tumor_groups <- infercnv_obj@observation_grouped_cell_indices
    res = list()
    normal_expr_data = infercnv_obj@expr.data[, unlist(infercnv_obj@reference_grouped_cell_indices)]
    for (tumor_group in names(tumor_groups)) {
        flog.info(sprintf("define_signif_tumor_subclusters(), tumor: %s",
            tumor_group))
        tumor_group_idx <- tumor_groups[[tumor_group]]
        tumor_expr_data <- infercnv_obj@expr.data[, tumor_group_idx]
        if (restrict_to_DE_genes) {
            p_vals <- .find_DE_stat_significance(normal_expr_data,
                tumor_expr_data)
            DE_gene_idx = which(p_vals < p_val)
            tumor_expr_data = tumor_expr_data[DE_gene_idx, ,
                drop = FALSE]
        }
        tumor_subcluster_info <- .single_tumor_subclustering(tumor_group,
            tumor_group_idx, tumor_expr_data, p_val, hclust_method,
            partition_method)
        res$hc[[tumor_group]] <- tumor_subcluster_info$hc
        res$subclusters[[tumor_group]] <- tumor_subcluster_info$subclusters
    }
    infercnv_obj@tumor_subclusters <- res
    if (!is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- define_signif_tumor_subclusters(infercnv_obj@.hspike,
            p_val, hclust_method, partition_method, restrict_to_DE_genes)
    }
    return(infercnv_obj)
}