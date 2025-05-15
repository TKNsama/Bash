# 优化后的 GO 富集分析函数
get_GO_enrich <- function(i, j) {
  message("Starting GO enrichment analysis for: ", i)
  
  # Step 1: 读取GO注释文件
  GO_file <- tryCatch({
    read_delim(j, col_names = TRUE)
  }, error = function(e) {
    message("Failed to read GO annotation file: ", e$message)
    return(NULL)
  })
  if (is.null(GO_file)) return(NULL)
  
  term2gene <- GO_file[, c("GO", "GeneID")]
  term2name <- unique(GO_file[, c("GO", "GO_term_name")])
  
  # Step 2: 读取基因列表
  gene_file_path <- paste0("U:/code/r_Pipeline/workspace/DEGs/", i)
  go_enrichment_gene <- tryCatch({
    df <- read_delim(gene_file_path, col_names = TRUE)[, 1]
    as.factor(df[[1]])
  }, error = function(e) {
    message("Failed to read gene list: ", e$message)
    return(NULL)
  })
  if (is.null(go_enrichment_gene)) return(NULL)
  
  # Step 3: 进行GO富集分析
  GO_enrich_analysis <- tryCatch({
    enricher(
      gene = go_enrichment_gene,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    message("GO enrichment analysis failed: ", e$message)
    return(NULL)
  })
  if (is.null(GO_enrich_analysis)) return(NULL)
  
  GO_enrich_analysis_data <- as.data.frame(GO_enrich_analysis)
  
  # Step 4: 创建输出目录
  out_dir <- paste0("GO/", i)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Step 5: Dotplot
  tryCatch({
    pdf(file.path(out_dir, "GOenrich_dotplot.pdf"), width = 10, height = 10)
    print(dotplot(GO_enrich_analysis, showCategory = 10))
    dev.off()
  }, error = function(e) {
    message("Dotplot error: ", e$message)
    try(dev.off(), silent = TRUE)
  })
  
  # Step 6: Barplot
  tryCatch({
    pdf(file.path(out_dir, "GOenrich_barplot.pdf"), width = 10, height = 10)
    print(barplot(GO_enrich_analysis, showCategory = 10))
    dev.off()
  }, error = function(e) {
    message("Barplot error: ", e$message)
    try(dev.off(), silent = TRUE)
  })
  
  # Step 7: 保存GO富集结果
  tryCatch({
    write.table(GO_enrich_analysis_data, file.path(out_dir, "GO_enrichment_results.txt"),
                row.names = FALSE, col.names = TRUE, sep = "\t")
  }, error = function(e) {
    message("Failed to write GO results: ", e$message)
  })
  
  # Step 8: GO分类条形图（BP, CC, MF）
  tryCatch({
    GO_classified <- GO_enrich_analysis_data[!is.na(GO_enrich_analysis_data$Description), ]
    GO_1_2_3 <- rbind(
      head(GO_classified[GO_classified$ONTOLOGY == "BP", ][order(GO_classified$p.adjust), ], 10),
      head(GO_classified[GO_classified$ONTOLOGY == "CC", ][order(GO_classified$p.adjust), ], 10),
      head(GO_classified[GO_classified$ONTOLOGY == "MF", ][order(GO_classified$p.adjust), ], 10)
    )
    
    pdf(file.path(out_dir, "GOenrich_second_class.pdf"), width = 15, height = 10)
    print(
      ggplot(GO_1_2_3, aes(x = reorder(Description, Count), y = Count, fill = -log10(p.adjust))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_gradient(low = "blue", high = "red") +
        facet_grid(ONTOLOGY ~ ., scale = "free") +
        labs(x = "Term", y = "Counts", fill = "FDR(-log10(P.adjust))") +
        theme(axis.title = element_text(size = 15), axis.text.y = element_text(size = 13))
    )
    dev.off()
  }, error = function(e) {
    message("Second classification plot error: ", e$message)
    try(dev.off(), silent = TRUE)
  })
  
  message("Finished GO enrichment for: ", i)
}
