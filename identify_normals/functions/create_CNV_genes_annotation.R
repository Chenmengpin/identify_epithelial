create_CNV_genes_annotation <- function(heat_df, gene_df) {

    heatmap_CNV_genes <- data.frame(gene = colnames(heat_df))

    # isolate CNV-associated genes in heatmap df:
    gene_df_sub <- gene_df[gene_df$Hugo_symbol %in% heatmap_CNV_genes$gene,]

    # if >30 genes, reduce to 30, gains and losses are prioritised:
    gene_df_sub <- gene_df_sub[1:30,]

    heatmap_CNV_genes$type <- "do_not_include"
    for ( t in c("gain", "loss", "either")) {
      heatmap_CNV_genes$type[heatmap_CNV_genes$gene %in% 
      gene_df_sub$Hugo_symbol[gene_df_sub$type == t]] <- t
    }

    label_subset <- which(heatmap_CNV_genes$type != "do_not_include")
    heatmap_CNV_labels <- as.character(heatmap_CNV_genes$gene[label_subset])

    # expand annotation bands by including a number of genes around GOI depending
    # on plot width:
    all_do_not_include <- which(heatmap_CNV_genes$type != "do_not_include")
    expand_no <- as.integer(ncol(heatmap_df)/1300)
    for (n in 1:expand_no) {
        for ( temp_index in all_do_not_include ) {
        
          heatmap_CNV_genes$type[temp_index + n] <- heatmap_CNV_genes$type[temp_index]
        
          if ( (temp_index - n) > 0 ) {
              heatmap_CNV_genes$type[temp_index - n] <- heatmap_CNV_genes$type[temp_index]
          }
        }
    }
    # create CNV_genes text annotation vector:
    heatmap_CNV_genes <- subset(heatmap_CNV_genes, select = -gene)

    CNV_genes_annot = HeatmapAnnotation(df = heatmap_CNV_genes, show_legend = F,
        col = list(type = c("gain" = "#F2210E", "loss" = "#3691BC", 
        "either"="black", "do_not_include"="#E7E4D3")),
        text = anno_link(label_subset, heatmap_CNV_labels, side = "bottom", 
        labels_gp = gpar(fontsize = 6, rot=2)))

    return(CNV_genes_annot)
}
