#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R

library(Seurat, lib.loc="/share/ScratchGeneral/jamtor/R/3.5.0")

#args = commandArgs(trailingOnly=TRUE)
#anchor_no <- as.numeric(args[1])
#dims_to_consider <- as.numeric(args[2])

sample_ids <- c("CID44041", "CID4495", "CID44971", "CID44972", "CID44991", 
  "CID44992", "CID4515", "CID43862", "CID43863" )

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
in_dir <- paste0(project_dir, "results/seurat/brca_mini_atlas_030719/")
seurat_filename <- "04_seurat_object_annotated.RData"


#########################################################################################
### 0. Load data ###
#########################################################################################

# load samples:
for (i in 1:length(sample_ids)) {

	print(paste0("Loading ", sample_ids[i], "..."))
	seurat_10X <- readRDS(paste0(in_dir, "seurat_", sample_ids[i], 
	  "/Output/Rdata/", seurat_filename))
#	cell_ids <- names(Idents(seurat_10X))
#	Idents(seurat_10X) <- paste0(sample_ids[i], "_", 
#		seurat_10X@meta.data$PC_A_res.1)

	if (i==1) {
		seurat_list <- list(seurat_10X)
	} else {
		seurat_list[[i]] <- seurat_10X
	}

	system(paste0("mkdir -p ", in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/"))
}
names(seurat_list) <- sample_ids

epithelial_clusters <- lapply(seurat_list, function(x) {
	return(grep("pithelial", levels(Idents(x)), value=T))
})

i=1
lapply(seurat_list, function(x) {
	system(paste0("mkdir -p ", in_dir, "seurat_", sample_ids[i], "/Output/Plots/"))
	png(
    paste0(in_dir, "seurat_", sample_ids[i], "/Output/Plots/epithelial_marker_violin_plot.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_violinplot <- VlnPlot(
	      object = x,
	      features = c("EPCAM", "KRT18", "ESR1", "KRT5", "KRT14", "ELF5", "PGR", 
	      	"ERBB2", "MKI67"),
	      pt.size = 1.5,
	      idents = epithelial_clusters[[i]]
	    )
	    
	    print(temp_violinplot)
    dev.off()

	i <<- i+1
})

i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Plots/epithelial_marker_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c("EPCAM", "KRT18", "ESR1", "KRT5", "KRT14", "ELF5", "PGR", 
	      	"ERBB2", "MKI67"),
	      pt.size = 0.75,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})




i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/custom_panel_markers1_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c("CCL5", "LAG3", "CD3E", "CXCL10", "IFNGR1", "CD27", "NKG7", 
	      	"CD4", "DKK2", "IL12B", "VSIR", "CD274", "PSMB10", "CD40", "EPCAM", "IL15" 
	    ),
	      pt.size = 1,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})

i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/custom_panel_markers2_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c("PDCD1", "CD276", "PDCD1LG2", "CD40LG", "FAS", "IL6", "PECAM1", "CD8A", 
	      	"STAT1", "CD44", "FOXP3", "ITGAM", "PTEN", "CMKLR1", "TIGIT", "CD47"
	    ),
	      pt.size = 1,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})

i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/custom_panel_markers3_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c( 
	      	"GZMB", "ITGAV", "PTPRC", "CXCL9", "AKT1", "CD68", "HAVCR2", "ITGAX", 
	      	"STAT2", "CXCR6", "ARG1", "CD74", "HIF1A", "ITGB2", "STAT3", "HLA-DQA1"
	    ),
	      pt.size = 1,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})

i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/custom_panel_markers4_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c( 
	      	"B2M", "CD86", "ICAM1", "ITGB8", "TBX21", "HLA-DRB1", "BATF3", "CSF1R", 
	      	"ICOSLG", "LY6E", "TNF", "HLA-E", "BCL2", "CTLA4", "IFNAR1", "MKI67"
	    ),
	      pt.size = 1,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})

i=1
lapply(seurat_list, function(x) {
  png(
    paste0(in_dir, "seurat_", sample_ids[i], 
    	"/Output/Figures/Gene_markers/TSNE/custom_panel_markers5_tSNE.png"), 
    width = 14, height = 8, res=300,
    	units = 'in')
		temp_featureplot <- FeaturePlot(
	      object = x,
	      features = c( 
	      	"TNFRSF9", "IDO1", "CCND1", "CTNNB1", "IFNG", "MS4A1", "VEGFA", "UBB", 
	      	 "OAZ1", "SDHA", "POLR2A", "RAB7A"
	    ),
	      pt.size = 1,
	      reduction = "TSNEA",
	      order = T
	    )
	    
	    print(temp_featureplot)
    dev.off()

	i <<- i+1
})


i=1
lapply(seurat_list, function(x) {
	png(
    paste0(
	  "/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/results/seurat/seurat_", 
	  	sample_ids[i], "/epithelial_metrics_violin_plot.png"), width = 14, height = 8, res=300,
    	units = 'in')
		temp_violinplot <- VlnPlot(
	      object = x,
	      features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
	      pt.size = 1.5,
	      idents = epithelial_clusters[[i]]
	    )
	    
	    print(temp_violinplot)
    dev.off()

	i <<- i+1
})



annots <- list(
	as.character(c(0, 1, 3, 5, 8, 13)),
	as.character(c(2, 3, 7, 9, 10, 14, 19)),
	as.character(c(2, 3, 4, 6, 11, 12, 14, 15, 20))
)

for (j in 1:length(seurat_list)) {
	idents <- levels(Idents(seurat_list[[j]]))
	idents <- gsub(paste0(sample_ids[j], "_"), "", idents)
	levels(Idents(seurat_list[[j]]))[idents %in% annots[[j]]] <- paste0("Epithelial_", 
		idents[idents %in% annots[[j]]])
	levels(Idents(seurat_list[[j]]))[!(idents %in% annots[[j]])] <- paste0("Non-epithelial_", 
		idents[!(idents %in% annots[[j]])])

	saveRDS(seurat_list[[j]], paste0(out_path, sample_ids[j], "/04_seurat_object_annotated.RData"))
}





