# Analysis Instructions

## Init Container

Startup the rstudio container in an envoironment of your choice.

## Analysis

- Analysis requires a pre-requisite cell atlas model for RCTD to run correctly. We used GSE157079 for this analysis, specifically P0 adult counts and clusters
- Assuming spatial data has been process via a 10x or 10x like file structure, the files for each sample should be arranged as barcode, features and matrix `.tsv.gz` files in `_dge` suffixed folders. Adjust the folder names in the setup paths sections of `run_rctd_sample.R` to reflect where the source data resides
- Run `run_rctd_sample.R`. This will create an RCTD rds object
- Adjust the setup paths in `puck_analysis.R`
- Run `puck_analysis.R`, This will generate QC plots at the puck/bead level for all samples and dump and data and plots in the corresponding output directories
- Adjust the setup paths in `gene_analysis.R`
- Run `gene_analysis.R`, This will generate a seurat object for all samples, generate QC plots at the gene level for all samples and dump and data and plots in the corresponding output directories
- To run differential analysis between pucks, edit `gen_deseq2.R` for the required comparison and run. This script will output the deseq2 results for further analysis
