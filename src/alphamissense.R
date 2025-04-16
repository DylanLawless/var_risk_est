# # Load required libraries
# library(dplyr)
# library(stringr)
# library(ggplot2)
# 
# manual version ----
# https://alphafold.ebi.ac.uk/files/AF-P15918-F1-hg38.csv
# 
# wget https://alphafold.ebi.ac.uk/files/AF-P15918-F1-hg38.csv > ../data/output/AF-P15918-F1-hg38.csv
# 
# gene <- "RAG1"
# 
# # rag1
# system('wget https://alphafold.ebi.ac.uk/files/AF-P15918-F1-hg38.csv -O ../output/AF-P15918-F1-hg38.csv')
# system('wget https://alphafold.ebi.ac.uk/files/AF-P15918-F1-aa-substitutions.csv -O ../output/AF-P15918-F1-aa-substitutions.csv')
# 
# # import
# af_hg38 <- read.table("../output/AF-P15918-F1-hg38.csv", sep = ",", header = TRUE)
# af <- read.table("../output/AF-P15918-F1-aa-substitutions.csv", sep = ",", header = TRUE)
# df <- readRDS(file = "../output/VarRiskEst_PanelAppRex_ID_398_gene_variants.Rds")
# df <- df |> filter(genename == gene)
# 
# # Now, parse HGVSc_VEP to extract REF and ALT:
# df_parsed <- df %>%
#   mutate(
#     REF = str_match(HGVSc_VEP, "c\\.\\d+([ACGT])>([ACGT])")[,2],
#     ALT = str_match(HGVSc_VEP, "c\\.\\d+([ACGT])>([ACGT])")[,3]
#   )
# 
# df_parsed$POS <- df_parsed$`pos(1-based)`
# 
# af_hg38 <- af_hg38 |> select(protein_variant, CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id)
# 
# af_merge <- merge(af, af_hg38, by = c("protein_variant"),   all.x = TRUE)
# 
# # Parse the protein_variant column into original, resnum, and alt
# af_merge_parsed <- af_merge %>%
#   mutate(
#     original = str_sub(protein_variant, 1, 1),
#     resnum   = as.numeric(str_extract(protein_variant, "[0-9]+")),
#     alt      = str_sub(protein_variant, -1, -1)
#   ) |> unique()
# 
# df_af <- merge(af_merge_parsed , df_parsed, by = c("POS", "REF", "ALT"), all.x =  TRUE)
# 
# # Transform occurrence_prob so that 0 becomes 0 and 1 becomes 6
# df_af <- df_af %>%
#   mutate(prob_trans = log10(occurrence_prob + 1e-5) + 5)
# 
# p1 <- ggplot(df_af, aes(x = resnum, y = prob_trans, fill = am_pathogenicity)) +
#   geom_col(position = "identity") +
#   geom_point(aes(y = prob_trans), shape = 21, size = 3, colour = "black", position = "identity") +
#   scale_fill_gradientn(
#     colours = c("navy", "lightblue", "red"),
#     limits  = c(0, 1),
#     na.value = "grey80"
#   ) +
#   labs(
#     x = "Residue sequence number",
#     y = "log10(probability)*\nadjusted",
#     fill = "Pathogenicity"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_blank())
# 
# p1
# 
# # p2: complete matrix as before
# p2 <- ggplot(df_af, aes(x = resnum, y = alt, fill = am_pathogenicity)) +
#   geom_tile() +
#   scale_fill_gradientn(
#     colours = c("navy", "lightblue", "red"),
#     limits  = c(0, 1),
#     na.value = "grey80"
#   ) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   labs(
#     x = "Residue sequence number",
#     y = "Alternative amino acid",
#     fill = "Pathogenicity"
#   ) +
#   theme_minimal(base_size = 12) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# library(patchwork)
# p_alphamissense <-
#   p1 / p2 +
#   plot_layout(guides = 'collect', axis = "collect") +
#   plot_annotation(tag_levels = 'A')
# 
# p_alphamissense
# ggsave("../images/p_alphamissense.png", plot = p_alphamissense, width = 8, height = 5)
# 



# function ----


library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)

process_gene_data <- function(gene, up_id,
                              panel_file = "../output/VarRiskEst_PanelAppRex_ID_398_gene_variants.Rds") {
  ## Step 1: Download the AlphaFold files if they do not exist
  file_hg38 <- sprintf("../output/AF-%s-F1-hg38.csv", up_id)
  file_aa   <- sprintf("../output/AF-%s-F1-aa-substitutions.csv", up_id)
  
  if (!file.exists(file_hg38))
    system(sprintf("wget https://alphafold.ebi.ac.uk/files/AF-%s-F1-hg38.csv -O %s", up_id, file_hg38))
  if (!file.exists(file_aa))
    system(sprintf("wget https://alphafold.ebi.ac.uk/files/AF-%s-F1-aa-substitutions.csv -O %s", up_id, file_aa))
  
  ## Step 2: Read the AlphaFold data files
  af_hg38 <- read.table(file_hg38, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  af <- read.table(file_aa, sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  ## Step 3: Read panel variants and filter by gene
  df <- readRDS(panel_file) %>% filter(genename == gene)
  
  ## Step 4: Parse HGVSc_VEP to extract REF and ALT; get POS from the 1-based position column
  df_parsed <- df %>% 
    mutate(
      REF = str_match(HGVSc_VEP, "c\\.\\d+([ACGT])>([ACGT])")[,2],
      ALT = str_match(HGVSc_VEP, "c\\.\\d+([ACGT])>([ACGT])")[,3]
    )
  df_parsed$POS <- df_parsed$`pos(1-based)`
  
  ## Step 5: Restrict af_hg38 columns to the needed ones
  af_hg38 <- af_hg38 %>% 
    select(protein_variant, CHROM, POS, REF, ALT, genome, uniprot_id, transcript_id)
  
  ## Step 6: Merge the aa substitutions file (af) with af_hg38 info using protein_variant
  af_merge <- merge(af, af_hg38, by = "protein_variant", all.x = TRUE)
  
  ## Step 7: Parse protein_variant into original, resnum and alt
  af_merge_parsed <- af_merge %>%
    mutate(
      original = str_sub(protein_variant, 1, 1),
      resnum   = as.numeric(str_extract(protein_variant, "[0-9]+")),
      alt      = str_sub(protein_variant, -1, -1)
    ) %>% unique()
  
  ## Step 8: Merge the parsed AlphaFold data with panel data using common POS, REF, ALT
  df_af <- merge(af_merge_parsed, df_parsed, by = c("POS", "REF", "ALT"), all.x = TRUE)
  
  ## Step 9: Transform occurrence_prob (so that 0 maps to 0 and 1 to 6) 
  df_af <- df_af %>% mutate(prob_trans = log10(occurrence_prob + 1e-5) + 5)
  
  df_af
}

plot_alphamissense <- function(df_af, output_path) {
  p1 <- ggplot(df_af, aes(x = resnum, y = prob_trans, fill = am_pathogenicity)) +
    geom_col(position = "identity") +
    geom_point(aes(y = prob_trans), alpha= 0.9, shape = 21, size = 3, colour = "black", position = "identity") +
    scale_fill_gradientn(
      colours = c("navy", "lightblue", "red"),
      limits  = c(0, 1),
      breaks  = c(0, 0.5, 1),
      labels  = c("Benign", "Unknown", "Pathogenic"),
      na.value = "grey80"
    ) +
    labs(
      x = "Residue sequence number",
      y = "log10(probability)*\nadjusted",
      fill = "Pathogenicity",
      title = gene
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank())
  
  p2 <- ggplot(df_af, aes(x = resnum, y = alt, fill = am_pathogenicity)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c("navy", "lightblue", "red"),
      limits  = c(0, 1),
      breaks  = c(0, 0.5, 1),
      labels  = c("Benign", "Unknown", "Pathogenic"),
      na.value = "grey80"
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "Residue sequence number",
      y = "Alternative\namino acid",
      fill = "Pathogenicity"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  
  p_alphamissense <- p1 / p2 +
    plot_layout(guides = "collect", axis = "collect", heights = c(1,1.5)) +
    plot_annotation(tag_levels = "A")
  
  ggsave(output_path, plot = p_alphamissense, width = 10, height = 5)
  print(p_alphamissense)
}

## Example usage for RAG1:
gene <- "RAG1"
df_af_gene <- process_gene_data(gene, "P15918")
p_gene <- plot_alphamissense(df_af_gene, paste0("../images/p_alphamissense_", gene, ".pdf"))

## Example usage for nfkb1:
gene <- "NFKB1"
df_af_nfkb1 <- process_gene_data(gene, "P19838")
p_gene <- plot_alphamissense(df_af_nfkb1, paste0("../images/p_alphamissense_", gene, ".pdf"))
# cat("Summary of AlphaMissense-derived pathogenicity scores (am_pathogenicity):\n")
# print(summary(df_af_nfkb1$am_pathogenicity))
# cat("\nSummary of observed probability of disease case (occurrence_prob):\n")
# print(summary(df_af_nfkb1$occurrence_prob))
# ks_result_nfkb1 <- ks.test(df_af_nfkb1$am_pathogenicity, df_af_nfkb1$occurrence_prob)
# cat("\nKolmogorov–Smirnov test result for NFKB1:\n")
# print(ks_result_nfkb1)

## Example usage for cftr:
gene <- "CFTR"
df_af_cftr <- process_gene_data(gene, "P13569")
p_gene <- plot_alphamissense(df_af_cftr, paste0("../images/p_alphamissense_", gene, ".pdf"))
# cat("Summary of AlphaMissense-derived pathogenicity scores (am_pathogenicity):\n")
# print(summary(df_af_nfkb1$am_pathogenicity))
# cat("\nSummary of observed probability of disease case (occurrence_prob):\n")
# print(summary(df_af_nfkb1$occurrence_prob))
# ks_result_nfkb1 <- ks.test(df_af_nfkb1$am_pathogenicity, df_af_nfkb1$occurrence_prob)
# cat("\nKolmogorov–Smirnov test result for NFKB1:\n")
# print(ks_result_nfkb1)



# statistical test ----
df_af_nfkb1
df_af_cftr

plot_occurrence_by_class <- function(df_af, gene) {
  # Perform the Kruskal-Wallis test by clinical classification
  kruskal_result <- kruskal.test(occurrence_prob ~ am_class, data = df_af)
  annot_text <- paste0(gene, ": ", sprintf("Kruskal-Wallis chi-squared = %.3f,\ndf = %d, p = %.2e - Not signif.",
                        kruskal_result$statistic, kruskal_result$parameter, kruskal_result$p.value))
  print(kruskal_result)
  
  p_box <- ggplot(df_af, aes(x = factor(am_class, levels = c("LBen", "Amb", "LPath")), y = prob_trans, fill = am_class)) +
    geom_boxplot(outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c("LBen" = "navy",
                                 "Amb"  = "lightblue",
                                 "LPath"= "red")) +
    scale_x_discrete(labels = c("LBen" = "Benign", "Amb" = "Unknown", "LPath" = "Pathogenic")) +
    labs(subtitle = annot_text,
         x = "Classification",
         y = "Observed Disease\nProbability") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p_box)
}

# Example usage:
# For NFKB1:
p1 <- plot_occurrence_by_class(df_af_nfkb1, "NFKB1")
# For CFTR:
p2 <- plot_occurrence_by_class(df_af_cftr, "CFTR")

p_alpha_kw <- (p1 + p2) +
  plot_layout(guides = "collect", axis = "collect") +
  plot_annotation(tag_levels = "A")

p_alpha_kw

ggsave(file = paste0("../images/p_alphamissense_kw", ".pdf"), plot = p_alpha_kw, width = 10, height = 6)

