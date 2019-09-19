lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(rtracklayer, lib.loc=lib_loc)
library(GenomicRanges)

project_name <- "identify_epithelial"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")

all_data <- read.table(
  "/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/raw_files/axiom_array/CN_v1_cnv_smooth_signal.cn", 
  header=T, sep="\t"
)

brca_data <- all_data[
  ,grep("Pros|Blood|P[0-9]m|4411|4413|Control", 
  colnames(all_data),invert=T)
]
brca_gr <- GRanges(
  seqnames = Rle(paste0("chr", brca_data$Chromosome)),
  ranges = IRanges(start=brca_data$Position, end = brca_data$Position),
  strand = "*",
  probe = brca_data$ProbeSet,
  symbol = NA
)

#brca_gr_sub <- brca_gr[c(1,2,66,666,6666, 66666, 666666)]

gene_coords <- read.table(
  paste0(ref_dir, "infercnv_gene_order.txt"), header=F, sep="\t"
)
gene_gr <- GRanges(
  seqnames = Rle(gene_coords$V2),
  ranges = IRanges(start=gene_coords$V3, end = gene_coords$V4),
  strand = "*",
  symbol = gene_coords$V1
)

print(paste0("Original number of array probes = ", length(brca_gr)))

# find overlaps between 2 genomic ranges to find which probes correspond to which genes:
olaps <- findOverlaps(brca_gr, gene_gr)
# remove all genes except for first hit for all probes:
olaps <- olaps[!duplicated(queryHits(olaps)),]
# add corresponding gene symbols to each probe entry in brca_gr:
brca_gr$symbol[queryHits(olaps)] <- as.character(gene_gr$symbol[subjectHits(olaps)])
# remove all probes which do not correspond to gene symbol:
brca_gr <- brca_gr[!is.na(brca_gr$symbol)]

print(paste0("Final number of array probes with gene information = ", length(brca_gr)))

### 349003 probes not assigned to any gene symbols! ###

