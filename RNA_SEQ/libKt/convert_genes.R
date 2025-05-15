# Convert gene IDs to Arabidopsis gene IDs
convert_genes <- function(gene_list, species) {
  orthofinder_data <- read.delim("data/Orthogroups.tsv", header = TRUE, sep = "\t")
  
  species_columns <- list(
    Arabidopsis = "Arabidopsis_thaliana_v10.0_representative_gene_model_updated",
    Barley_v2 = "Hordeum_vulgare_Morex_HC_v2.0",
    Barley_v3 = "Hordeum_vulgare_Morex_HC_v3.0",
    Rice = "Oryza_sativa_HC_v7.0",
    Maize = "Zea_mays_B73_NAM_HC_v5.0",
    Wheat_v1.1 = "Triticum_aestivum_CS_HC_v1.1_filtered",
    Wheat_v2.1 = "Triticum_aestivum_CS_HC_v2.1_filtered"
  )
  
  if (!species %in% names(species_columns)) {
    stop("Species not found in the orthofinder data mapping. Please check the species name.")
  }
  
  species_column <- species_columns[[species]]
  
  orthofinder_data <- orthofinder_data %>%
    mutate(
      Arabidopsis_genes = str_split(Arabidopsis_thaliana_v10.0_representative_gene_model_updated, ","),
      species_genes = str_split(!!sym(species_column), ",") # 动态列名
    )
  
  matched_groups <- orthofinder_data %>%
    filter(map_lgl(species_genes, ~ any(. %in% gene_list)))
  
  if(nrow(matched_groups) == 0) {
    message("No matching groups found for the provided gene list.")
  }
  
  unique_arabidopsis_genes <- matched_groups %>%
    pull(Arabidopsis_genes) %>%
    unlist() %>%
    unique()
  
  unique_arabidopsis_genes <- gsub("\\..*", "", unique_arabidopsis_genes)
  
  unique_arabidopsis_genes <- unique(unique_arabidopsis_genes)
  
  if(length(unique_arabidopsis_genes) == 0) {
    message("No matching Arabidopsis genes found.")
  }
  
  return(unique_arabidopsis_genes)
}
