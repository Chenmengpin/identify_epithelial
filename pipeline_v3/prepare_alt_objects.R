# prepare and combine 'alt' datasets which do not have enough stromal cells with stromal cells
# from datasets which do

library(Seurat)
library(DropletUtils)

project_name <- "identify_epithelial"
subproject_name <- "normal_threshold"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")

alt_ids <- c("IND4", "IND5", "IND6", "IND7")
stromal_ids <- c("CID44991", "CID44992", "CID4471", "CID4530")

annotate_identities <- function(seurat_object) {
  annotated_idents <- gsub(
    " ", "_", paste0(seurat_object@meta.data$garnett_cluster_call, " ", 
    Idents(seurat_object))
  )
  # only keep annotation with largest no cells per existing cluster
  annotated_ident_numbers <- table(annotated_idents)
  annotated_ident_cluster_numbers <- gsub("^.*_CID", "CID", 
    names(annotated_ident_numbers))
  annotated_ident_cluster_numbers <- gsub("^.*_IND", "IND", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*_PDX", "PDX", 
    annotated_ident_cluster_numbers)

  for ( c in levels(Idents(seurat_object)) ) {
    annotated_ident_index <- which(annotated_ident_cluster_numbers == c)
    potential_clusters <- annotated_ident_numbers[annotated_ident_index]
    correct_cluster <- 
      names(potential_clusters[potential_clusters == 
        max(potential_clusters)])
      # label cells in seurat object as the correct cluster
      levels(Idents(seurat_object))[levels(Idents(seurat_object)) == c] <- 
        correct_cluster
  }
  # change "Epithelial_Luminal" to "Epithelial":
  levels(Idents(seurat_object)) <- gsub("Epithelial_Luminal", "Epithelial", 
    levels(Idents(seurat_object)))
  return(seurat_object)
}

for (i in 1:length(alt_ids)) {

  print(paste0("Processing primary ", stromal_ids[i], " and ALT ", 
  	ids[i], "..."))

  input_dir <- paste0(seurat_path, "seurat_", ids[i], 
  	"/Output/InferCNV/input_files/") 
  system(paste0("mkdir -p ", input_dir))

  alt <- readRDS(paste0(seurat_path, "/seurat_", ids[i], 
  	"/Output/Rdata/03_seurat_object_processed.RData"))
  Idents(alt) <- paste0(ids[i], "_", alt@meta.data$PC_A_res.0.8)
  alt@meta.data$garnett_cluster_call <- "Epithelial"
  
  prim <- readRDS(paste0(seurat_path, "/seurat_", stromal_ids[i], 
  	"/Output/Rdata/03_seurat_object_processed.RData"))
  Idents(prim) <- paste0(stromal_ids[i], "_", prim@meta.data$PC_A_res.0.8)

  # downsample ALT to meet mean nUMI/cell of primary:
  # fetch mean nUMI counts:
  mean_prim_count <- mean(apply(GetAssayData(prim , slot = "counts"), 2, 
  	sum))
  mean_alt_count <- mean(apply(GetAssayData(alt , slot = "counts"), 2, 
  	sum))
  print(paste0("Mean nUMI/cell in primary: ", mean_prim_count))
  print(paste0("Mean nUMI/cell in ALT: ", mean_alt_count))

  if (mean_prim_count < mean_alt_count) {
    # downsample samples to mean_prim_count:
    scaling_factor <- mean_prim_count/mean_alt_count
    new_raw_matrix <- downsampleMatrix(GetAssayData(alt, slot = "counts"), 
      scaling_factor)
    alt@assays$RNA@counts <- new_raw_matrix
    print(paste0("Downsampled mean nUMI/cell in ALT: ", 
      mean(apply(alt@assays$RNA@counts, 2, sum))))
  } else if (mean_alt_count < mean_prim_count) {
    # downsample samples to mean_alt_count:
    scaling_factor <- mean_alt_count/mean_prim_count
    new_raw_matrix <- downsampleMatrix(GetAssayData(prim, slot = "counts"), 
      scaling_factor)
    prim@assays$RNA@counts <- new_raw_matrix
    print(paste0("Downsampled mean nUMI/cell in primary: ", 
      mean(apply(prim@assays$RNA@counts, 2, sum))))
  }

  prim_annotated <- annotate_identities(prim)
  alt_annotated <- annotate_identities(alt)

  # fetch cell ids from clusters to keep:
  clusters_to_remove <-  grep("pithelial", levels(Idents(prim_annotated)), value=T)
  cells_to_keep <- names(Idents(prim_annotated))[!(Idents(prim_annotated) %in% 
  	clusters_to_remove)]

  # subset prim_annotated:
  prim_annotated <- SubsetData(prim_annotated, cells=cells_to_keep)

  # merge objects:
  joined_obj <- merge(prim_annotated, alt_annotated, merge.data=FALSE)

  saveRDS(joined_obj, paste0(seurat_path, "/seurat_", ids[i], 
  	"/Output/Rdata/04_seurat_object_combined.RData"))

}
