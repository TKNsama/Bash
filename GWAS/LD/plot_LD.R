.libPaths("/filer-5/agruppen/PBP/tan/R_LIBS/Library")

library(tidyverse)
library(patchwork)
library(stringr)

plot_gwas_gene_ld <- function(
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
  pdf_height = 8
) {

  gene_format <- match.arg(gene_format)

  # ===============================
  # 1. GWAS
  # ===============================
  gwas <- read_tsv(assoc_file, show_col_types = FALSE) %>%
    rename(SNP = rs) %>%
    mutate(
      ps = as.numeric(gsub(",", "", ps)),
      P  = as.numeric(p_wald),
      P  = ifelse(is.na(P) | P <= 0, 1e-300, P),
      logP = -log10(P)
    ) %>%
    filter(chr == !!chr, ps >= start, ps <= end)

  if (nrow(gwas) == 0)
    stop("No SNPs in this region.")

  if (!lead_snp_id %in% gwas$SNP)
    stop("lead_snp_id not found in GWAS.")

  lead_snp <- gwas %>% filter(SNP == lead_snp_id)

  # ===============================
  # 2. LD（相对于 lead SNP）
  # ===============================
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

  gwas <- gwas %>%
    left_join(ld_lead, by = "SNP")

  # ===============================
  # 3. Gene annotation
  # ===============================
  if (gene_format == "bed") {
    genes <- read_tsv(gene_file, col_names = FALSE, show_col_types = FALSE) %>%
      setNames(c("chr", "start", "end", "gene"))
  }

  if (gene_format == "gff") {
    genes <- read_tsv(
      gene_file,
      comment = "#",
      col_names = FALSE,
      show_col_types = FALSE
    ) %>%
      setNames(c(
        "chr","source","type","start","end",
        "score","strand","phase","attribute"
      )) %>%
      filter(type == "gene") %>%
      mutate(
        gene = str_extract(attribute, "ID=[^;]+") %>%
               str_remove("ID=")
      )
  }

  # 统一 chr 名（防止 Chr6A / chr6A）
  genes <- genes %>%
    mutate(chr_clean = str_remove(chr, "^(Chr|chr)")) %>%
    filter(
      chr_clean == !!chr,
      start <= end,
      end >= start
    )

  if (nrow(genes) == 0)
    warning("No genes found in this region.")
  # 给 gene 分层（避免重叠）
  genes <- genes %>%
    arrange(start) %>%
    mutate(track = row_number() %% 3)
  # ===============================
  # 4. GWAS + LD 图
  # ===============================
  p_gwas <- ggplot(gwas, aes(ps, logP)) +
    geom_point(aes(color = R2), size = 1.2) +
    geom_point(
      data = lead_snp,
      color = "red",
      size = 3
    ) +
    scale_color_gradient(
      low = "grey80",
      high = "red",
      na.value = "grey90",
      name = expression(r^2)
    ) +
    scale_x_continuous(limits = c(start, end)) +
    labs(
      x = NULL,
      y = expression(-log[10](p)),
      title = paste0(chr, ":", start/1e6, "-", end/1e6, " Mb")
    ) +
    theme_bw()

  # ===============================
  # 5. Gene track 图
  # ===============================
  p_gene <- ggplot(genes) +
    geom_rect(
      aes(
        xmin = start,
        xmax = end,
        ymin = track,
        ymax = track + 0.6
      ),
      fill = "steelblue"
    ) +
    geom_text(
      aes(
        x = (start + end) / 2,
        y = track + 0.7,
        label = gene
      ),
      size = 3,
      angle = 45,
      hjust = 0,
      check_overlap = TRUE
    ) +
    scale_x_continuous(limits = c(start, end)) +
    scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
    labs(
      x = paste0("Chr ", chr, " position"),
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank()
    )

  # ===============================
  # 6. 组合并输出
  # ===============================
  pdf(out_pdf, width = pdf_width, height = pdf_height)
  print(
    p_gwas /
    p_gene +
    plot_layout(heights = c(3, 1))
  )
  dev.off()

  message("Saved PDF: ", out_pdf)
}
plot_gwas_gene_ld(
  assoc_file  = "output/likehood_1_2/1_res.assoc.txt",
  ld_file     = "LD_chr1B.ld",
  gene_file   = "/filer-5/agruppen/PBP/tan/wheat/cs_IAAS/CS-IAAS_v1.1_HC.gff3",
  gene_format = "gff",
  chr         = "1B",
  start       = 22000000,
  end         = 24000000,
  lead_snp_id = "1B_22859488",
  out_pdf     = "Chr1B_22m_24m_GWAS_GENE_LD.pdf"
)

