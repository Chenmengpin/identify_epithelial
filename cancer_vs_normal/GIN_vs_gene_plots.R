#/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R

library(ggplot2)
library(Seurat, lib.loc="/share/ScratchGeneral/jamtor/R/3.5.0")

#args = commandArgs(trailingOnly=TRUE)
#anchor_no <- as.numeric(args[1])
#dims_to_consider <- as.numeric(args[2])

sample_ids <- c("CID3586", "CID3921", "CID3941", "CID3948", "CID3963", 
  "CID4066", "CID4067", "CID4290A", "CID43862", "CID43863", 
  "CID4398", "CID44041", "CID4461", "CID4463", "CID4465", 
  "CID4471", "CID4495", "CID44971", "CID44972", "CID44991", 
  "CID44992", "CID4515", "CID45171", "CID45172", 
  "CID4520N", "CID4523", "CID4523N", "CID4530", "CID4530N", 
  "CID4535")

ER_pos <- c("CID3948", "CID3941", "CID4067", 
  "CID4290A", "CID4461", "CID4463", 
  "CID4465", "CID4471", "CID4530", 
  "CID4535")

HER2_pos <- c("CID3586", "CID4066", "CID4398", "CID3921", 
  "CID3963", "CID45171", "CID45172")

TNBC <- c("CID44041", "CID4495", "CID44971", 
	"CID44991", "CID44992", "CID4515", 
	"CID4386", "CID43862", "CID43863", 
	"CID44042", "CID4408", "CID4409", 
	"CID44972")

all_subsets <- list(ER_pos, HER2_pos, TNBC)

GOI <- c("KRT14", "MKI67", "ESR1", "ELF5", "CCND1")

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/identify_epithelial/")
in_dir <- paste0(project_dir, "results/seurat/brca_mini_atlas_030719/")
seurat_filename <- "03_seurat_object_processed.RData"


#########################################################################################
### 1. Plot samples individually ###
#########################################################################################

# load GIN scores:
GIN_scores <- read.table(paste0(in_dir, "/brca_atlas_epithelial_GIN_scores.txt"),
    header=T, as.is=T)

for (i in 1:length(sample_ids)) {

  print(paste0("Loading ", sample_ids[i], "..."))
  seurat_10X <- readRDS(paste0(in_dir, "seurat_", sample_ids[i], 
    "/Output/Rdata/", seurat_filename))
    
  if (i==1) {
    seurat_list <- list(seurat_10X)
  } else {
    seurat_list[[i]] <- seurat_10X
  }
}
names(seurat_list) <- sample_ids

j=1
lapply(seurat_list, function(x) {
  print(j)
  out_dir <- paste0(in_dir, "seurat_", sample_ids[j], "/Output/Plots/")
  system(paste0("mkdir -p ", out_dir))

  GIN <- GIN_scores[grep(sample_ids[j], GIN_scores$cell_id),]
  RNA_matrix <- x@assays$RNA

  for (k in 1:length(GOI)) {
    if (k==1) {
      gene_expression_df <- data.frame(cell_id = colnames(RNA_matrix),
        expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
    } else {
      gene_expression_df <- rbind(
        gene_expression_df,
        data.frame(cell_id = colnames(RNA_matrix),
        expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
      )
    }
  }
  colnames(gene_expression_df) <- c("cell_id", "expression", "GOI")
    
  genes_vs_GIN <- merge(gene_expression_df, GIN, by="cell_id")
    
  p <- ggplot(genes_vs_GIN, aes(GIN, expression, colour = GOI))
  p <- p + geom_smooth(method=lm)
  p <- p + geom_point()
  
  pdf(paste0(out_dir, sample_ids[j], "_epithelial_genes_vs_GIN.pdf"))
    print(p)
  dev.off()
  j <<- j+1

})


#########################################################################################
### 2. Plot samples by subset ###
#########################################################################################

all_subsets <- list(ER_pos, HER2_pos, TNBC)
names(all_subsets) <- c("ER_pos", "HER2_pos", "TNBC")

for (i in 1:length(all_subsets)) {
  for (j in 1:length(all_subsets[[i]])) {
    print(paste0("Loading ", all_subsets[[i]][j], "..."))
    seurat_10X <- readRDS(paste0(in_dir, "seurat_", all_subsets[[i]][j], 
      "/Output/Rdata/", seurat_filename))

    GIN <- GIN_scores[gsub("_.*$", "", GIN_scores$cell_id) == all_subsets[[i]][j],]
    RNA_matrix <- seurat_10X@assays$RNA

    for (k in 1:length(GOI)) {
      if (k==1) {
        gene_expression_df <- data.frame(cell_id = colnames(RNA_matrix),
          expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
      } else {
        gene_expression_df <- rbind(
          gene_expression_df,
          data.frame(cell_id = colnames(RNA_matrix),
          expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
        )
      }
    }

    colnames(gene_expression_df) <- c("cell_id", "expression", "GOI")
    gene_expression_df$id <- all_subsets[[i]][j]

    if (j==1) {
      genes_vs_GIN <- merge(gene_expression_df, GIN, by="cell_id")
    } else {
      genes_vs_GIN[[j]] <- rbind(
        genes_vs_GIN, merge(gene_expression_df, GIN, by="cell_id")
      )
    }
  }

}



  for (i in 1:length(sub)) {

    print(paste0("Loading ", sub[i], "..."))
    seurat_10X <- readRDS(paste0(in_dir, "seurat_", sub[i], 
      "/Output/Rdata/", seurat_filename))
      
    if (i==1) {
      seurat_list <- list(seurat_10X)
    } else {
      seurat_list[[i]] <- seurat_10X
    }
  }
  names(seurat_list) <- sub
  
  lapply(seurat_list, function(seurat) {
  
    GIN <- GIN_scores[grep(x[j], GIN_scores$cell_id),]
    RNA_matrix <- y@assays$RNA
  
    for (k in 1:length(GOI)) {
      if (k==1) {
        gene_expression_df <- data.frame(cell_id = colnames(RNA_matrix),
          expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
      } else {
        gene_expression_df <- rbind(
          gene_expression_df,
          data.frame(cell_id = colnames(RNA_matrix),
          expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
        )
      }
    }
    colnames(gene_expression_df) <- c("cell_id", "expression", "GOI")

    genes_vs_GIN <- merge(gene_expression_df, GIN, by="cell_id")
    return(genes_vs_GIN)
  })
      


    
      
    p <- ggplot(genes_vs_GIN, aes(GIN, expression, colour = GOI))
    p <- p + geom_smooth(method=lm)
    p <- p + geom_point()
    
    pdf(paste0(out_dir, sample_ids[j], "_epithelial_genes_vs_GIN.pdf"))
      print(p)
    dev.off()
    j <<- j+1
  
  })
})



#########################################################################################
### 3. Plot one line per subset ###
#########################################################################################

for (n in 1:length(all_subsets)) {

  print(paste0("n=", n))
  subset_object <- seurat_list[names(seurat_list) %in% all_subsets[[n]]]
  subset_ids <- names(subset_object)

  m=1
  subset_genes_vs_GIN_temp <- lapply(subset_object, function(x) {
  	print(paste0("m=", m))
    GIN <- GIN_scores[grep(subset_ids[m], GIN_scores$cell_id),]
    RNA_matrix <- x@assays$RNA
  
    for (k in 1:length(GOI)) {
      print(paste0("k=", k))
      if (length(grep(GOI[k], rownames(RNA_matrix))) != 0) {
      	if (k==1) {
          gene_expression_df <- data.frame(cell_id = colnames(RNA_matrix),
            expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
        } else {
          gene_expression_df <- rbind(
            gene_expression_df,
          	data.frame(cell_id = colnames(RNA_matrix),
            expression = as.numeric(RNA_matrix[GOI[k],]), GOI[k])
          )
        }
      }
    }
    colnames(gene_expression_df) <- c("cell_id", "expression", "GOI")
      
    genes_vs_GIN <- merge(gene_expression_df, GIN, by="cell_id")

    m <<- m+1
    return(genes_vs_GIN)
  })

  subset_genes_vs_GIN <- do.call("rbind", subset_genes_vs_GIN_temp)

  p <- ggplot(subset_genes_vs_GIN, aes(GIN, expression, colour = GOI))
  p <- p + geom_smooth(method=lm)
  
  pdf(paste0(in_dir, names(all_subsets)[n], "_epithelial_genes_vs_GIN.pdf"))
    print(p)
  dev.off()

}
