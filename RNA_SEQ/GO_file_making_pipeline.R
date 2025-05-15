# File Name: Make GO file for panBarlex.R
# Author: Kenan Tan
# Date: 2025-05-08
# Description: Description of the script
library(tidyverse)
setwd("U:/code/r_Pipeline/workspace/")
# 读取原始注释文件（包含列名，tab分隔）
df <- read.table("res/bowman_annotations.tsv", sep = "\t", header = FALSE, quote = "\"", stringsAsFactors = FALSE)

# 第2列是 Gene ID，第5列是 GO terms（可能为空）
gene_ids <- df$V2
go_terms <- df$V5

# 删除空GO项的行
valid_rows <- which(go_terms != "")
gene_ids <- gene_ids[valid_rows]
go_terms <- go_terms[valid_rows]

# 按逗号拆分 GO term
gene_go_list <- strsplit(go_terms, ",")

# 展开为两列格式
gene_id_expanded <- rep(gene_ids, sapply(gene_go_list, length))
go_term_expanded <- unlist(gene_go_list)

# 构建初步数据框
go_df <- data.frame(GeneID = gene_id_expanded, GO = go_term_expanded, stringsAsFactors = FALSE)

# --- 添加 GO 描述 和 类别 ---
library(GO.db)

# 添加描述
go_df$GO_term_name <- Term(go_df$GO)

# 添加类别（ONTOLOGY: BP/MF/CC）
go_df$ONTOLOGY <- Ontology(go_df$GO)

# 查看结果
head(go_df)

# 写入最终文件
write.table(go_df, file = "gene2go_clusterprofiler.txt", sep = "\t", quote = FALSE, row.names = FALSE)

