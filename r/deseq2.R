rm(list = ls())
system("type R")

# Load libraries
library(spacexr)
library(Matrix)
library(dplyr)
library(tibble)
library(patchwork)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(DESeq2)

# Params
project.name <- "2022_03_novaseq_PM22082"
output.name <- "sham_vs_0.5gy_GCR"
#output.name <- "sham_vs_5gy_GRB"
#output.name <- "sham_vs_1gy_SPE"


# Init
target_conditions <- c("untreated","untreated","untreated", "treated", "treated", "treated")

#                           1-F41-CTRL   2-F50-CTRL   3-F44-CTRL   
samples_ids_sham_ctrl <- c("211013_23", "211125_11", "211013_03")

#                     1-F61-GRB    2-F68-GRB    3-F65-GRB   
samples_ids_grb <- c("211125_22", "211125_04", "211125_21")

#                     1-F01-GCR    2-F08-GCR    3-F03-GCR   
samples_ids_gcr <- c("211013_21", "211125_20", "211004_12")

#                     1-F21-SPE    2-F27-SPE    3-F22-SPE   
samples_ids_spe <- c("211004_02", "211125_18", "211013_17")



target_samples <- c(samples_ids_sham_ctrl, samples_ids_gcr)

# Tissue types
# tissue_types <- c("DCT", "Early PT", "Endo", "IC", "LOH", "Macro", "Neutro", "NP", "PC", "PCT", "Podo", "Proliferating", "PST", "Stroma 1")
tissue_types <- c("LOH")

# Setup paths
data_dir <- "2022_space_kidney/data"
project_dir <- paste0(data_dir, "/", project.name)
deseq_output_dir <- paste0(data_dir, "/", project.name, "_gen_202303/deseq")

# Loop tissue type
for(j in 1:length(tissue_types)) {
  t_type <- tissue_types[j]
  
  name = paste0(output.name, "_", t_type)
  deseq_output_path <- paste0(deseq_output_dir, "/", name, ".rds")
  
  df_gene_sums <- readRDS(deseq_output_path)
  
  rownames(df_gene_sums) <- df_gene_sums$gene_id
  df_gene_sums$gene_id <- NULL
  
  df_sample_data = data.frame(target_samples, target_conditions)
  rownames(df_sample_data) <- df_sample_data$target_samples
  df_sample_data$target_samples <- NULL
  colnames(df_sample_data) <- c("condition")
  
  dds <- DESeqDataSetFromMatrix(countData = df_gene_sums, colData = df_sample_data, design = ~ 0 + condition)
  dds <- DESeq(dds)
  
  resultsNames(dds) # lists the coefficients
  res <- results(
    dds,
    name="condition_untreated_vs_treated",
    lfcThreshold = 0,
    altHypothesis = "greaterAbs",
    independentFiltering = TRUE,
    alpha = 0.1,
    pAdjustMethod = "BH",
    minmu = 0.5,
    contrast = c("condition", c("untreated", "treated"))
  )
  resOrdered <- res[order(res$padj, res$pvalue, decreasing = FALSE),]
  head(resOrdered)
  
  # or to shrink log fold changes association with condition:
  #res_lfc <- lfcShrink(dds, coef="condition_untreated_vs_treated", type="apeglm")
  
#  pdf(paste0(deseq_output_dir, "/", name, "_ma.pdf"))
#  plot <- plotMA(res, ylim=c(-2,2))
#  print(plot)
#  dev.off()
  
 # pdf(paste0(deseq_output_dir, "/", name, "_ma_lfc.pdf"))
#  plot <- plotMA(res_lfc, ylim=c(-2,2))
#  print(plot)
#  dev.off()
  
  write.csv(resOrdered,paste0(deseq_output_dir, "/", name, "_data.csv"), row.names = TRUE)
}

