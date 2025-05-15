# File Name: rna-seq_pipeline_v3.R
# Author: Kenan Tan
# Date: 2024-01-08
# Description: Comprehensive RNA-seq analysis pipeline for differential expression, visualization, enrichment, and clustering

# ------------------------
# Load necessary libraries
# ------------------------
library(DESeq2)
library(tidyverse)
library(PCAtools)
library(pheatmap)
library(clusterProfiler)
library(KEGGREST)
library(org.At.tair.db)
library(VennDiagram)
library(RColorBrewer)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)
library(cluster)

# ------------------------
# Set working directory
# ------------------------
setwd("U:/code/r_Pipeline/workspace/")

# ------------------------
# Load input data
# ------------------------
count_matrix <- read.csv("vrs4_count_matrix.txt", row.names = 1, sep = "\t")
coldata <- read.csv("ExpDesignMetaData.txt", sep = "\t")
tpm_matrix <- read.csv("vrs4_tpm_matrix.txt", row.names = 1, sep = "\t")

# ------------------------
# Experiment setting
# ------------------------
the_species <- "Barley_v2"
comparisons <- list(
  c("VR_DR", "BW_DR"), c("VR_TM", "BW_TM"),
  c("VR_AP", "BW_AP"), c("VR_SP", "BW_SP")
)

# ------------------------
# Define custom functions
# ------------------------
source("U:/code/r_function/diff_gene_list.R")
source("U:/code/r_function/generate_volcano_plot.R")
source("U:/code/r_function/get_GO_enrich.R")
source("U:/code/r_function/All_DEGs_id.R")
source("U:/code/r_function/convert_genes.R")

# ------------------------
# Basic DESeq2 analysis
# ------------------------
count_matrix <- count_matrix[rowSums(count_matrix) > 10, ]
count_matrix <- round(count_matrix)
tpm_matrix <- tpm_matrix[rowSums(tpm_matrix) > 1, ]

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

# ------------------------
# PCA plot
# ------------------------
vst_data <- assay(vst(dds))
rownames(coldata) <- coldata$sample
p <- pca(vst_data, metadata = coldata, removeVar = 0.1)
write.csv(p$rotated, "pca_df.csv")
screeplot(p)
biplot(p, x = 'PC1', y = 'PC2')

# ------------------------
# TPM heatmap
# ------------------------
pheatmap(tpm_matrix, scale = "row", show_rownames = FALSE, kmeans_k = 10)

# ------------------------
# DEGs identification and volcano plots
# ------------------------
dir.create("DEGs", showWarnings = FALSE)
for (comp in comparisons) {
  comp_name <- paste0(comp[1], "_vs_", comp[2])
  diff_gene_list(logfc = 1, p_adj = 0.05, comparison = comp, name = comp_name)
}

dir.create("volcano_plot", showWarnings = FALSE)
for (comp in comparisons) {
  comp_name <- paste0(comp[1], "_vs_", comp[2])
  generate_volcano_plot(results(dds, contrast = c("condition", comp)), comp_name)
}

# ------------------------
# GO enrichment analysis
# ------------------------
dir.create("GO", showWarnings = FALSE)
for (i in dir("DEGs", pattern = "_up.csv")) {
  get_GO_enrich(i, j = "U:/code/r_Pipeline/workspace/gene2go_clusterprofiler.txt")
}
for (i in dir("DEGs", pattern = "_down.csv")) {
  get_GO_enrich(i, j = "U:/code/r_Pipeline/workspace/gene2go_clusterprofiler.txt")
}

# ------------------------
# Select all DEGs for clustering
# ------------------------
DEGs_name <- All_DEGs_id("DEGs", output_file = "Intersection_all_DEGs_id.txt")
DEGs_tpm_matrix <- tpm_matrix[unlist(DEGs_name), ]
write.csv(DEGs_tpm_matrix, "DEGs_tpm_matrix.csv")

# ------------------------
# Venn diagram (example)
# ------------------------
venn_data <- list(
  FD_M29 = read.csv("DEGs/FD_TP_vs_M29_TP_all.csv")[, 1],
  FD_M30 = read.csv("DEGs/FD_TP_vs_M30_TP_all.csv")[, 1],
  FD_M32 = read.csv("DEGs/FD_TP_vs_M32_TP_all.csv")[, 1]
)
venn.diagram(
  x = venn_data,
  imagetype = "tiff",
  filename = "TP_all.tiff",
  fill = brewer.pal(3, "Set1")
)

# ------------------------
# UpSet plot (example)
# ------------------------
upset_data <- list(
  A = read.csv("DEGs/resW2.5_Spike_all")[, 1],
  B = read.csv("DEGs/resW3_Spike_all")[, 1],
  C = read.csv("DEGs/resW3.5_Spike_all")[, 1],
  D = read.csv("DEGs/resW2.5_Node_all")[, 1],
  E = read.csv("DEGs/resW3_Node_all")[, 1],
  F = read.csv("DEGs/resW3.5_Node_all")[, 1]
)
max_len <- max(sapply(upset_data, length))
upset_df <- as.data.frame(lapply(upset_data, function(x) c(x, rep(NA, max_len - length(x)))))
upset(fromList(upset_df), sets = names(upset_df), order.by = "freq", keep.order = TRUE)

# ------------------------
# Clustering
# ------------------------
# Elbow method
sse <- sapply(1:10, function(k) kmeans(DEGs_tpm_matrix, centers = k)$tot.withinss)
plot(1:10, sse, type = "b", pch = 19, col = "blue", main = "Elbow Method", xlab = "K", ylab = "SSE")

# Silhouette method
asw <- sapply(2:10, function(k) pam(DEGs_tpm_matrix, k)$silinfo$avg.width)
k.best <- which.max(asw) + 1
cat("Optimal K from silhouette:", k.best, "\n")

# K-means clustering
normalized_matrix <- t(scale(t(DEGs_tpm_matrix)))
pam_cluster <- pam(normalized_matrix, k.best)
pam_clust <- data.frame(pam_K = pam_cluster$cluster)
write.csv(pam_clust, "DEGs_heatmap_pam_cluster_info.csv")

# Heatmap
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(
  normalized_matrix,
  col = col_fun,
  cluster_columns = FALSE,
  column_split = data.frame(rep(c("1TP", "2HD"), each = 12)),
  row_split = pam_clust,
  show_row_names = FALSE,
  border = TRUE,
  show_heatmap_legend = FALSE
)
