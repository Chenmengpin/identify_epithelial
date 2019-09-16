#! /share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/Rscript

library(Seurat)
lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(infercnv, lib.loc=lib_loc)
library(HiddenMarkov, lib.loc=lib_loc)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
print(args)
sample_name <- args[1]
ncores=args[2]
h_vmem <- args[3]

#sample_name <- "CID3963"
#ncores <- "3"
#h_vmem <- "5M"
project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_030719"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")
sample_dir <- paste0(seurat_path, "/seurat_", sample_name, "/Output/")
setwd(sample_dir)
out_path <- paste0("InferCNV/supervised_clustering/subset_test/", ncores, 
	"_cores_", h_vmem, "_memory/")
func_dir <- paste0(project_dir, "/scripts/infercnv_subset_test_functions/")
ref_dir <- paste0(project_dir, "/refs/")


################################################################################
### 0. Define functions ###
################################################################################

prepare_infercnv_metadata <- dget(paste0(func_dir, "prepare_infercnv_metadata.R"))
create_infercnv_object <- dget(paste0(func_dir, "create_infercnv_object.R"))
run_infercnv <- dget(paste0(func_dir, "run_infercnv.R"))
run_infercnv <- possibly(run_infercnv, otherwise = "error_occurred")


################################################################################
### 1. Load seurat objects and subset ###
################################################################################

seurat_10X <- readRDS(paste0(seurat_path, "seurat_", sample_name, 
  "/Output/Rdata/03_seurat_object_processed.RData"))
Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.1

# generate vector of cell numbers to try:
cell_numbers <- c(100, 200, 400, 800, 1600, 3200, 6400)
#for (i in 1:9) {
#  if (i==1) {
#  	cell_numbers <- 1
#  } else if (i==2) {
#  	cell_numbers <- c(cell_numbers, 2)
#  } else if (i<4) {
#  	cell_numbers <- c(cell_numbers, (cell_numbers[i-1])^2)
#  } else {
#  	cell_numbers <- c(cell_numbers, (cell_numbers[i-1])*2)
#  }
#}

min_UMIs <- 1000
min_Genes <- 1000
# select good quality epithelial and stromal cells to include:
print(paste0("Original no cells: ", length(Idents(seurat_10X))))
quality_cells <- 
  names(Idents(seurat_10X))[seurat_10X@meta.data$nCount > min_UMIs & 
  seurat_10X@meta.data$nFeature > min_Genes]
sub_object1 <- subset(seurat_10X, cells=quality_cells)
potential_epithelial_cells <- names(Idents(sub_object1))[grep("pithelial", 
	sub_object1@meta.data$garnett_seurat_cluster_call_subset_PC_A_res.1)]
print(paste0("Number of cells in epithelial/myoepithelial clusters: ", 
  length(potential_epithelial_cells)))
EPCAM_expression <- FetchData(sub_object1, "EPCAM")
epithelial_cells <- rownames(EPCAM_expression)[
  EPCAM_expression$EPCAM > 0
]
KRT18_expression <- FetchData(sub_object1, "KRT18")
epithelial_cells <- c(epithelial_cells, rownames(KRT18_expression)[
  KRT18_expression$KRT18 > 0
])
KRT8_expression <- FetchData(sub_object1, "KRT8")
epithelial_cells <- c(epithelial_cells, rownames(KRT8_expression)[
  KRT8_expression$KRT18 > 0
])
KRT5_expression <- FetchData(sub_object1, "KRT5")
epithelial_cells <- c(epithelial_cells, rownames(KRT5_expression)[
  KRT5_expression$KRT5 > 0
])
KRT14_expression <- FetchData(sub_object1, "KRT14")
epithelial_cells <- unique(c(epithelial_cells, rownames(KRT14_expression)[
  KRT14_expression$KRT14 > 0
]))

epithelial_cells <- 
	potential_epithelial_cells[potential_epithelial_cells %in% epithelial_cells]
epithelial_cells <- unique(epithelial_cells)
potential_stromal_cells <- names(Idents(sub_object1))[grep("pithelial", 
	sub_object1@meta.data$garnett_cluster_call, invert=T)]
print(paste0("Number of cells in stromal clusters: ", 
  length(potential_stromal_cells)))
stromal_cells <- rownames(EPCAM_expression)[
  EPCAM_expression$EPCAM == 0
]
for (marker in c("KRT18", "KRT8", "KRT5", "KRT14", "ACTA2", "PDGFRA",
	"PDGFRB", "CD34", "MCAM")) {
	exp <- FetchData(sub_object1, marker)
	non_marker <- rownames(exp)[
    exp[,1] == 0
  ]
  stromal_cells <- stromal_cells[stromal_cells %in% non_marker]
}

# subset cells:
for (j in cell_numbers) {

  print(paste0("InferCNV beginning with ", j, " cells per group..."))

  input_dir <- paste0(out_path, "subset_", as.character(j), "_cells/input_files/")
  system(paste0("mkdir -p ", sample_dir, input_dir))

  # subset seurat object:
  cells_to_keep <- c(sample(epithelial_cells, cell_numbers[j]), 
  	sample(stromal_cells, cell_numbers[j]))
  if (length(epithelial_cells <= j)) {
  	 seurat_object <- subset(seurat_10X, cells = cells_to_keep)
  } else {
  	seurat_object <- seurat_10X
  }
  print(paste0("Total cell number = ", length(Idents(seurat_object))))


  ################################################################################
  ### 2. Generate input matrix and metadata files ###
  ################################################################################
  
  # create raw matrix input file
  
  print("Creating inferCNV raw counts file...")
  count_df <- as.matrix(GetAssayData(seurat_object , slot = "counts"))
  write.table(count_df, paste0(input_dir, "input_matrix.txt"), quote=F,
  sep="\t", col.names=T, row.names=T)

  # prepare infercnv metadata and write to files
  print("Creating inferCNV metadata file...")
  infercnv_metadata <- prepare_infercnv_metadata(seurat_object, FALSE,
    "PC_A_res.1", subset_data=FALSE, count_df)
  names(infercnv_metadata) <- c("metadata", "number_per_group", "seurat")
  seurat_object <- infercnv_metadata$seurat
  
  print(paste0("Idents are: ", levels(Idents(seurat_object))))
  write.table(infercnv_metadata$metadata, paste0(input_dir, "metadata.txt"), 
    quote=F, sep="\t", col.names=F, row.names=F)
  write.table(infercnv_metadata$number_per_group, paste0(input_dir, 
    "number_per_group.txt"), quote=F, col.names=F, row.names=F, sep="\t")
  
  normals <- grep("[e,E]pithelial|[m,M]yoepithelial|CAF|[u,U]nassigned", 
      unique(Idents(seurat_object)[names(Idents(seurat_object)) %in% colnames(count_df)]), value=T, 
      invert=T)
  
  print(paste0("Normals are: ", normals))
  

  ################################################################################
  ### 3. Run InferCNV ###
  ################################################################################

  print("Creating inferCNV object...")
  raw_path <- paste0(input_dir, "input_matrix.txt")
  annotation_path <- paste0(input_dir, "metadata.txt")
  gene_path <- paste0(ref_dir, "infercnv_gene_order.txt")
  initial_infercnv_object <- create_infercnv_object(raw_path, annotation_path,
    gene_path, normals)
  
  print("InferCNV object created, running inferCNV...")
  infercnv_output <- run_infercnv(initial_infercnv_object, numcores-1, out_dir, 0.1, 
    101, 3, 1.3, TRUE, "i6", subpop_mode=TRUE)
  print("InferCNV finished")
 
}







