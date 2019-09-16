library(Seurat)
project_name <- "identify_epithelial"
subproject_name <- "brca_mini_atlas_130819"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
seurat_path <- paste0(project_dir, "results/seurat/", subproject_name, "/")

integrated_object <- readRDS(
  "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Jun2019/02_integration/output/CCA_CCA23Jun2019/Output/Rdata/03_seurat_CCA_aligned_processed.Rdata")
)
integrated_epithelial <- names(Idents(integrated_object))[
  grep("pithelial", integrated_object@meta.data$garnett_call_ext_major)
]

sample_ids <- unique(gsub("\\_.*", "", integrated_epithelial))
sample_ids <- sample_ids[sample_ids != "CID44971CHUNKS2"]

for (s in sample_ids) {
  sample_epithelial <- grep(s, integrated_epithelial, value=T)
  saveRDS(sample_epithelial, 
  	paste0(seurat_path, "seurat_", s, "/integrated_epithelial.RData"))
}