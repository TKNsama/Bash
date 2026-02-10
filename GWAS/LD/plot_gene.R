.libPaths("/filer-5/agruppen/PBP/tan/R_LIBS/Library")

library(tidyverse)
library(stringr)

plot_gwas_gene_ld_single <- function(
  assoc_file,
  ld_file,
  gene_file,
  gene_format = c("bed", "gff"),
  chr,
  start,
  end,
  lead_snp_id,
  out_pdf,
  pdf_width = 8,
  pdf_height = 6
) {

  gene_format <- match.arg(gene_format)

  ### 1. GWAS
  gwas <- read_tsv(assoc_file, show_col_types = FALSE) %>%
    rename(SNP = rs) %>%
    mutate(
      ps = as.numeric(gsub(",", "", ps)),
      P  = as.numeric(p_wald),
      P  = ifelse(is.na(P) | P <= 0, 1e-300, P),
      logP = -log10(P)
    ) %>%
    filter(chr == !!chr, ps >= start, ps <= end)

  if (!lead_snp_id %in% gwas$SNP)
    stop("lead_snp_id not found.")

  ### 2. LD relative to lead SNP
  ld <- read_table(ld_file, col_types = cols())

  ld_lead <- ld %>%
    filter(
      (SNP_A == lead_snp_id & BP_B >= start & BP_B <= end) |
      (SNP_B == lead_snp_id & BP_A >= start & BP_A <= end)
    ) %>%
    mutate(
      SNP = ifelse(SNP_A == lead_snp_id, SNP_B, SNP_A)
    ) %>%
    select(SNP, R2)

  gwas <- gwas %>% left_join(ld_lead, by = "SNP")

  ### 3. genes
  if (gene_format == "bed") {
    genes <- read_tsv(gene_file, col_names = FALSE) %>%
      setNames(c("chr", "start", "end", "gene")) %>%
      filter(chr == !!chr, start <= end, end >= start)
  }

  if (gene_format == "gff") {
    genes <- read_tsv(
      gene_file, comment = "#",
      col_names = FALSE, show_col_types = FALSE
    ) %>%
      setNames(c(
        "chr","source","type","start","end",
        "score","strand","phase","attribute"
      )) %>%
      filter(type == "gene", chr == !!chr) %>%
      mutate(
        gene = str_extract(attribute, "ID=[^;]+") %>%
               str_remove("ID=")
      )
  }

  ### 4. 给 gene 分层（避免重叠）
  genes <- genes %>%
    arrange(start) %>%
    mutate(
      y = -1 - (row_number() %% 2) * 0.4
    )

  ### 5. 画在同一张图
p <- ggplot() +

  ## SNPs
  geom_point(
    data = gwas,
    aes(ps, logP, color = R2),
    size = 1.2
  ) +

  ## lead SNP（修复点）
  geom_point(
    data = subset(gwas, SNP == lead_snp_id),
    aes(ps, logP),
    color = "red",
    size = 3
  ) +

  ## gene blocks
  geom_rect(
    data = genes,
    aes(
      xmin = start,
      xmax = end,
      ymin = y,
      ymax = y + 0.3
    ),
    fill = "steelblue"
  ) +

  ## gene labels
  geom_text(
    data = genes,
    aes(
      x = (start + end) / 2,
      y = y + 0.35,
      label = gene
    ),
    size = 3,
    angle = 45,
    hjust = 0,
    check_overlap = TRUE
  ) +

  scale_color_gradient(
    low = "grey80",
    high = "red",
    na.value = "grey90",
    name = expression(r^2)
  ) +

  scale_x_continuous(limits = c(start, end)) +
  theme_bw()

  ### 6. 输出
  pdf(out_pdf, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()

  message("Saved PDF: ", out_pdf)
}

plot_gwas_gene_ld_single(
  assoc_file  = "output/likehood_1_2/1_filter.tsv",
  ld_file     = "LD_chr6A.ld",
  gene_file   = "/filer-5/agruppen/PBP/tan/wheat/cs_IAAS/CS-IAAS_v1.1_HC.gff3",
  gene_format = "gff",
  chr         = "6A",
  start       = 560000000,
  end         = 600000000,
  lead_snp_id = "6A_581244838",
  out_pdf     = "Chr6A_560_600Mb_GWAS_GENE_LD.pdf"
)


