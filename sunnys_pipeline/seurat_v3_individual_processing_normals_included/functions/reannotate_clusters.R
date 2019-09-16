library(Seurat)

sample_name <- "seurat_IND4"

seurat_10X <- readRDS(paste0("/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/normal_threshold/", 
	sample_name, "/Output/Rdata/03_seurat_object_processed.RData"))

Idents(seurat_10X) <- seurat_10X@meta.data$PC_A_res.0.8

#IND4:
cell_types <- c("Epithelial", "Myoepithelial", "Myoepithelial", "Epithelial", 
	"Epithelial", "Epithelial", "Unknown")

#IND5:
#cell_types <- c("Myoepithelial", "Myoepithelial", "Epithelial", "Myoepithelial",
#	"Epithelial", "Myoepithelial", "Epithelial", "Myoepithelial", "Epithelial")

#IND6:
#cell_types <- c("Myoepithelial", "Myoepithelial", "Myoepithelial", "Epithelial", 
#	"Myoepithelial", "Myoepithelial", "Epithelial", "Myoepithelial", "Myoepithelial", 
#	"Unknown", "Unknown")

#IND7:
#cell_types <- rep("Epithelial", length(levels(Idents(seurat_10X) )))



names(cell_types) <- levels(Idents(seurat_10X))

for ( i in levels(Idents(seurat_10X)) ) {
	garnett <- as.character(seurat_10X@meta.data$garnett_cluster_call)
	print(i)
	print(table(garnett[Idents(seurat_10X) == i]))
	garnett[Idents(seurat_10X) == i] <- cell_types[i]
	seurat_10X@meta.data$garnett_cluster_call <- factor(garnett)
	print(table(seurat_10X@meta.data$garnett_cluster_call[Idents(seurat_10X) == i]))
}

#for ( i in 1:(length(levels(Idents(seurat_10X)))-1) ) {
#	garnett <- as.character(seurat_10X@meta.data$garnett_cluster_call)
#	print(levels(Idents(seurat_10X))[i])
#	print(table(garnett[Idents(seurat_10X) == as.character(i)]))
#	garnett[Idents(seurat_10X) == as.character(i)] <- cell_types[i]
#	seurat_10X@meta.data$garnett_cluster_call <- factor(garnett)
#	print(table(seurat_10X@meta.data$garnett_cluster_call[Idents(seurat_10X) == as.character(i)]))
#}

saveRDS(seurat_10X, paste0("/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/normal_threshold/", 
	sample_name, "/Output/Rdata/03_seurat_object_processed.RData"))
