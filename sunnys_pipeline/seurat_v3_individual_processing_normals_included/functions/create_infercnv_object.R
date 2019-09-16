create_infercnv_object <- function(raw_counts_path, annotation_file_path, 
  gene_order_path, normals) {
  infercnv_obj = CreateInfercnvObject(
            raw_counts_matrix=raw_counts_path,
            annotations_file=annotation_file_path,
            delim="\t",
            gene_order_file=gene_order_path,
            ref_group_names=normals
  )
  return(infercnv_obj)
}
