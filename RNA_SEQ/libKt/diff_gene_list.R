# Differential gene extraction function
diff_gene_list <- function(logfc, comparison, name, p_adj) {
  res <- results(dds, contrast = c("condition", comparison))
  diff <- subset(res, padj < p_adj & abs(log2FoldChange) > logfc)
  diff_up <- subset(diff, log2FoldChange > logfc)
  diff_down <- subset(diff, log2FoldChange < -logfc)
  
  # Save results
  write.csv(diff, file = paste0("DEGs/",name, "_all.csv"))
  write.csv(diff_up, file = paste0("DEGs/",name, "_up.csv"))
  write.csv(diff_down, file = paste0("DEGs/",name, "_down.csv"))
}
