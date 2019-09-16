#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R

library(Seurat, lib.loc="/share/ScratchGeneral/jamtor/R/3.5.0")

#args = commandArgs(trailingOnly=TRUE)
#anchor_no <- as.numeric(args[1])
#dims_to_consider <- as.numeric(args[2])

sample_ids <- c("CID4520N", "IND4", "IND5", "IND6", "IND7")

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
in_dir <- paste0(project_dir, "results/seurat/brca_mini_atlas_030719/")
seurat_filename <- "03_seurat_object_processed.RData"
out_dir <- paste0(project_dir, "results/seurat/brca_mini_atlas_030719/seurat_normals/Output/Rdata/")
#, paste(sample_ids, collapse = "_"))
cluster_resolution <- "PC_A_res.1"


#########################################################################################
### 0. Load data ###
#########################################################################################

# load samples:
for (i in 1:length(sample_ids)) {

	seurat_10X <- readRDS(paste0(in_dir, "seurat_", sample_ids[i], "/Output/Rdata/",
		seurat_filename))
	if (i==1) {
		seurat_list <- list(seurat_10X)
	} else {
		seurat_list[[i]] <- seurat_10X
	}
}
names(seurat_list) <- sample_ids


#########################################################################################
### 1. Combine seurat objects ###
#########################################################################################

# exchange annotation in ident slots with that in 
# meta.data$PC_A_res.1
i=1
annotated_seurat_list <- lapply(seurat_list, function(x) {
  # change cluster resolution:
  Idents(x) <- paste0(sample_ids[i], "_", 
  	eval(parse(text=paste0("x@meta.data$", cluster_resolution))))
  # annotate with garnett annotations:
  annotated_idents <- gsub(
    " ", "_", paste0(
      eval(parse(text=paste0("x@meta.data$garnett_seurat_cluster_call_major_",
      cluster_resolution))), 
    " ", Idents(x))
  )
  # only keep annotation with largest no cells per existing cluster
  annotated_ident_numbers <- table(annotated_idents)
  annotated_ident_cluster_numbers <- gsub("^.*_CID", "CID", 
    names(annotated_ident_numbers))
  annotated_ident_cluster_numbers <- gsub("^.*_IND", "IND", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*_PDX", "PDX", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*[a-z]_", "", 
    annotated_ident_cluster_numbers)
  annotated_ident_cluster_numbers <- gsub("^.*CAF.*_", "", 
    annotated_ident_cluster_numbers)
  
  for ( c in levels(Idents(x)) ) {
    print(c)
    annotated_ident_index <- which(annotated_ident_cluster_numbers == c)
    potential_clusters <- annotated_ident_numbers[annotated_ident_index]
    correct_cluster <- 
      names(potential_clusters[potential_clusters == 
        max(potential_clusters)])
      # label cells in seurat object as the correct cluster
      levels(Idents(x))[levels(Idents(x)) == c] <- 
        correct_cluster
  }
  i <<- i+1
  return(x)
})

annotated_seurat_list_minus_first <- 
	annotated_seurat_list[2:length(annotated_seurat_list)]
merged_seurat <- merge(annotated_seurat_list[[1]], 
	annotated_seurat_list_minus_first, do.normalise = F)

system(paste0("mkdir -p ", out_dir))
saveRDS(merged_seurat, paste0(out_dir, "04_seurat_object_combined.RData"))


#########################################################################################
### 2. Integrate data using Seurat v3 CCA ###
#########################################################################################

# find min gene number:
min_gene_no <- min(unlist(lapply(annotated_seurat_list, 
	function(x) nrow(GetAssayData(x)))))

# integrate data:
variable_genes_for_CCA <- min_gene_no
num_cc_compute <- 20

print(date())

new_seurat_list <- lapply(seurat_list, function(x) {
	FindVariableFeatures(
      object = x,
      do.plot = F,
      nfeatures=variable_genes_for_CCA
    )
})

reference_list <- new_seurat_list[sample_ids]

seurat_anchors <- FindIntegrationAnchors(object.list = reference_list, dims = 1:num_cc_compute, 
	anchor.features = variable_genes_for_CCA)
print(date())
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:num_cc_compute)

print(date())

out_dir <- paste0(out_path, "_CCA/")
system(paste0("mkdir -p ", out_dir))
saveRDS(seurat_integrated, paste0(out_dir, "/05_seurat_object_integrated.RData"))
