lib_loc <- "/share/ScratchGeneral/jamtor/R/3.5dev/"
library(TxDb.Hsapiens.UCSC.hg19.knownGene, lib.loc=lib_loc)
library(GenomicFeatures, lib.loc=lib_loc)
library(OrganismDbi, lib.loc=lib_loc)
library(Homo.sapiens, lib.loc=lib_loc)
library(gwascat, lib.loc=lib_loc)
library(liftOver, lib.loc=lib_loc)
library(rtracklayer, lib.loc=lib_loc)
library(GenomicRanges)

project_name <- "identify_epithelial"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/single_cell/", project_name, "/")
ref_dir <- paste0(project_dir, "/refs/")

# load data:
all_data <- read.table(
  "/share/ScratchGeneral/jamtor/projects/single_cell/identify_epithelial/raw_files/axiom_array/CN_v1_cnv_smooth_signal.cn", 
  header=T, sep="\t"
)

# isolate only relevant fields and rename cols:
brca_data <- all_data[
  ,grep("Pros|Blood|P[0-9]m|4411|4413|Control", 
  colnames(all_data),invert=T)
]
colnames(brca_data) <- gsub("^.*KCC_P_", "PDX", colnames(brca_data))
colnames(brca_data) <- gsub("^.*A[0-9][0-9]_", "CID", colnames(brca_data))
colnames(brca_data) <- gsub("\\.CEL", "", colnames(brca_data))

# create granges object for data:
brca_hg19_gr <- GRanges(
  seqnames = Rle(paste0("chr", brca_data$Chromosome)),
  ranges = IRanges(start=brca_data$Position, end = brca_data$Position),
  strand = "*"
)
elementMetadata(brca_hg19_gr) <- subset(brca_data, select=-c(Chromosome, Position))
brca_hg19_gr$symbol <- NA

# retrieve hg19 to hg38 liftover object:
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

seqlevelsStyle(brca_hg19_gr) = "UCSC"
brca_gr_list <- liftOver(brca_hg19_gr, ch)
brca_gr <- unlist(brca_gr_list)

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
array_CNVs <- subset(
  as.data.frame(elementMetadata(brca_gr)), select = -ProbeSet
)

# find mean CNV value for duplicated genes:
split_CNVs <- split(array_CNVs, array_CNVs$symbol)
#table(unlist(lapply(split_CNVs, nrow)))
final_array_CNVs <- do.call(
  "rbind", lapply(split_CNVs, function(x) {
    x <- subset(x, select=-symbol)
    return(apply(x, 2, mean))
  })
)


print(paste0("Final number of array probes with gene information = ", length(brca_gr)))

# save df:
write.table(final_array_CNVs, paste0(ref_dir, "all_array_CNVs.txt"),
  row.names=T, col.names=T, sep="\t", quote=F)


### 349003 probes not assigned to any gene symbols! ###

