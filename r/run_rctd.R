rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(plyr)
library(dplyr)
library(EnsDb.Mmusculus.v79)

# Setup paths
data_dir <- "data"
ref_matrix.path <- paste0(data_dir, "/GSE157079_P0_adult_counts.rds")
ref_clusters.path <- paste0(data_dir, "/GSE157079_P0_adult_clusters.txt")
barcodes.path <- paste0(data_dir, "/barcodes.tsv.gz")
features.path <- paste0(data_dir, "/features.tsv.gz")
matrix.path <- paste0(data_dir, "/matrix.mtx.gz")

# Load reference table
ref_counts <- readRDS(ref_matrix.path)
ref_meta <- read.table(ref_clusters.path,sep="\t", header=TRUE)

# Get gene symbols and convert to gene ids and match for found
ref_gene_names <- rownames(ref_counts)
ref_gene_ids <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ref_gene_names, keytype="SYMBOL", columns=c("SYMBOL","GENEID"))
matched_names <- match(ref_gene_names, ref_gene_ids$SYMBOL)

# Drop rows that arnt found
rownames(ref_counts) <- ref_gene_ids$GENEID[matched_names]
drop <- which(is.na(rownames(ref_counts)))
ref_counts <- ref_counts[-drop,]

# Get cell type mappings and create reference object
cell_types <- ref_meta$clusters; names(cell_types) <- ref_meta$barcodes
cell_types <- as.factor(cell_types)
reference <- Reference(ref_counts, cell_types)

# Load count matrix and annotate with gene ids
count_mat <- readMM(file = matrix.path)
features.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcodes.names = read.delim(barcodes.path, header = FALSE, stringsAsFactors = FALSE)
colnames(count_mat) = barcodes.names$V1
rownames(count_mat) = features.names$V1

# Create coordinate and nUMI inputs for puck
coords <- read.csv(file.path(data_dir,"211004_02.csv"))
coords_filt <- select(coords, c(1,3,4))
rownames(coords_filt) <- coords_filt[,1]; coords_filt[,1] <- NULL
nUMI <- colSums(count_mat)

# Create puck and RCTD objects and run RCTD
puck <- SpatialRNA(coords_filt, count_mat, nUMI)
rctd <- create.RCTD(puck, reference, max_cores = 8)
rctd <- run.RCTD(rctd, doublet_mode = 'doublet')

# Save rctd object
saveRDS(rctd, file = "data/rctd_results.rds")
