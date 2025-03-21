library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2)





# clinvardata <- "~/Desktop/clinvar/variant_summary_head.txt"
# clinvardata <- "~/Desktop/clinvar/variant_summary.txt"
# output <- "~/Desktop/clinvar/"
clinvar_sum <- "~/web/iei_genetics/data/"
clinvar_out <- "~/web/iei_genetics/output/"

gene_variant_summary_clean <- readRDS(file = paste0(clinvar_sum, "gene_variant_summary.Rds"))


# Define a scoring map for individual classification terms.
# Higher positive numbers indicate a more concerning (negative) classification,
# while lower (or negative) numbers indicate a less concerning (or "good") classification.
score_map <- c(
  "Pathogenic" = 5,
  "Likely pathogenic" = 4,
  "Pathogenic, low penetrance" = 3,
  "likely pathogenic, low penetrance" = 2,
  "Benign" = -5,
  "Likely benign" = -4,
  "drug response" = 0,
  "Uncertain significance" = 0,
  "Conflicting classifications of pathogenicity" = 2,
  "risk factor" = 1,
  "association" = 1,
  "no classification for the single variant" = 0,
  "no classifications from unflagged records" = 0,
  "Affects" = 0,
  "other" = 0,
  "not provided" = 0,
  "protective" = -3,
  "likely risk allele" = 1,
  "uncertain risk allele" = 0
)

# Read the header line, remove the leading '#' and split into column names
header_line <- read_lines(clinvardata, n_max = 1)
col_names <- strsplit(sub("^#", "", header_line), "\t")[[1]]

# Function to process each chunk: count ClinicalSignificance per GeneSymbol
process_chunk <- function(chunk, pos) {
  chunk %>%
    group_by(GeneSymbol, ClinicalSignificance) %>%
    summarise(n = n(), .groups = "drop")
}

# List to store results from each chunk
accumulated <- list()

# Callback function for processing each chunk
chunk_callback <- function(x, pos) {
  accumulated[[length(accumulated) + 1]] <<- process_chunk(x, pos)
}

# Read the file in chunks, skipping the header line
read_delim_chunked(
  clinvardata,
  delim = "\t",
  col_names = col_names,
  skip = 1,
  chunk_size = 100000,
  callback = DataFrameCallback$new(chunk_callback)
)

# Combine chunk results and tally counts per GeneSymbol and ClinicalSignificance
final_results <- bind_rows(accumulated) %>%
  group_by(GeneSymbol, ClinicalSignificance) %>%
  summarise(total = sum(n), .groups = "drop")

# print(final_results)

final_results$ClinicalSignificance |> unique()

# Get unique ClinicalSignificance strings from the dataset
unique_clin <- unique(final_results$ClinicalSignificance)

# Replace ";" with "/" and split each string into individual terms
all_terms <- unique(unlist(str_split(str_replace_all(unique_clin, ";", "/"), "/")))
all_terms <- str_trim(all_terms)

# Convert both sets to lower case for case-insensitive matching
mapped_terms <- tolower(names(score_map))
found_terms <- tolower(all_terms)

# Identify terms that are not in the scoring map
unmatched <- setdiff(found_terms, mapped_terms)

if (length(unmatched) > 0) {
  message("The following classification terms are not accounted for in the score_map:")
  print(unmatched)
} else {
  message("All classification terms are accounted for in the score_map.")
}


# Per gene tally score ----
# Function to compute a rank score for a given ClinicalSignificance string.
# It splits the string by both "/" and ";" delimiters and averages the scores.
rank_clin_sig <- function(clin_str) {
  # Replace ";" with "/" so we have one common delimiter.
  unified <- str_replace_all(clin_str, ";", "/")
  # Split by "/"
  terms <- unlist(str_split(unified, "/"))
  # Trim whitespace from each term
  terms <- str_trim(terms)
  # Look up each term in score_map (case-insensitive) and assign a score.
  # If a term isn't found, assign a score of 0.
  scores <- sapply(terms, function(term) {
    key <- names(score_map)[str_to_lower(names(score_map)) == str_to_lower(term)]
    if (length(key) > 0) {
      return(score_map[[key]])
    } else {
      return(0)
    }
  })
  # Return the average score as the rank for the string.
  mean(scores)
}

# Extract unique ClinicalSignificance terms from the final results.
unique_terms <- unique(final_results$ClinicalSignificance)

# Compute the rank for each unique term.
ranked_terms <- data.frame(
  ClinicalSignificance = unique_terms,
  Rank = sapply(unique_terms, rank_clin_sig),
  stringsAsFactors = FALSE
)

print(ranked_terms)


# plot ranks ----
# First, join the final_results with ranked_terms to add the computed rank per ClinicalSignificance entry.
# Assume final_results and ranked_terms are already defined from previous steps.
df_with_rank <- final_results %>%
  left_join(ranked_terms, by = "ClinicalSignificance")

# Now, for each gene (unique GeneSymbol), compute:
#  - total variant count (sum of 'total')
#  - weighted average rank (weighted by the variant count)
gene_summary <- df_with_rank %>%
  group_by(GeneSymbol) %>%
  summarise(
    total_variants = sum(total),
    avg_rank = sum(Rank * total) / sum(total)
  ) %>%
  ungroup()

# First plot: Plot average rank (x-axis) vs. count of variants per gene
# Use shape 21 with black outline and fill based on avg_rank.
p1 <- ggplot(gene_summary, aes(x = avg_rank, y = total_variants, fill = avg_rank)) +
  geom_point(shape = 21, size = 3, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    x = "Average Rank Score",
    y = "Count of Variants per Gene",
    fill = "Avg Rank"
  ) +
  theme_bw()

print(p1)

# Second plot: Normalize the average rank (zero-center) and plot again.
mean_rank <- mean(gene_summary$avg_rank)
gene_summary <- gene_summary %>% mutate(norm_avg_rank = avg_rank - mean_rank)

p2 <- ggplot(gene_summary, aes(x = norm_avg_rank, y = total_variants, fill = norm_avg_rank)) +
  geom_point(shape = 21, size = 3, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    x = "Normalized Average Rank Score (zero centered)",
    y = "Count of Variants per Gene",
    fill = "Norm Avg Rank"
  ) +
  theme_bw()

print(p2)

patch1 <- p1 / p2

# Create histogram of avg_rank: Number of genes per average rank bin
p1_hist <- ggplot(gene_summary, aes(x = avg_rank)) +
  geom_histogram(bins = 30, color = "black", aes(fill = ..x..)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    x = "Average Rank Score",
    y = "Number of Genes",
    fill = "Avg Rank"
  ) +
  theme_bw()

# Normalize avg_rank
mean_rank <- mean(gene_summary$avg_rank)
gene_summary <- gene_summary %>% mutate(norm_avg_rank = avg_rank - mean_rank)

# Create histogram of norm_avg_rank: Number of genes per normalized rank bin
p2_hist <- ggplot(gene_summary, aes(x = norm_avg_rank)) +
  geom_histogram(bins = 30, color = "black", aes(fill = ..x..)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    x = "Normalized Average Rank Score (zero centered)",
    y = "Number of Genes",
    fill = "Norm Avg Rank"
  ) +
  theme_bw()

# Combine the two plots vertically
patch2 <- p1_hist / p2_hist

# Per gene pathogenic count ----
# Create a new column categorizing variants as "Pathogenic" if their computed rank equals 5,
# and "Other" otherwise.
df_with_rank <- df_with_rank %>%
  mutate(variant_class = if_else(round(Rank, 1) == 5, "Pathogenic", "Other"))

# Aggregate counts per gene and classification
gene_variant_summary <- df_with_rank %>%
  group_by(GeneSymbol, variant_class) %>%
  summarise(variant_count = sum(total), .groups = "drop") 

# split clean ----
# This script explores and cleans the GeneSymbol field in gene_variant_summary by splitting rows on ";" and removing duplicates, then aggregating variant counts by gene and classification.

library(dplyr)
library(stringr)
library(tidyr)

rows_with_semicolon <- gene_variant_summary %>% filter(str_detect(GeneSymbol, ";"))
rows_without_semicolon <- gene_variant_summary %>% filter(!str_detect(GeneSymbol, ";"))
cat("Number of rows with semicolon:", nrow(rows_with_semicolon), "\n")
cat("Number of rows without semicolon:", nrow(rows_without_semicolon), "\n")
split_genes <- rows_with_semicolon %>% separate_rows(GeneSymbol, sep = ";") %>% mutate(GeneSymbol = str_trim(GeneSymbol))
unique_no_semicolon <- unique(rows_without_semicolon$GeneSymbol)
unique_split <- unique(split_genes$GeneSymbol)
common_genes <- intersect(unique_no_semicolon, unique_split)
cat("Number of unique gene symbols without semicolon:", length(unique_no_semicolon), "\n")
cat("Number of unique gene symbols from split rows:", length(unique_split), "\n")
cat("Number of common gene symbols:", length(common_genes), "\n")
print(common_genes)
print(unique_split)

gene_variant_summary_clean <- gene_variant_summary %>%
  separate_rows(GeneSymbol, sep = ";") %>%
  mutate(GeneSymbol = str_trim(GeneSymbol)) %>%
  distinct() %>%
  group_by(GeneSymbol, variant_class) %>%
  summarise(variant_count = sum(variant_count), .groups = "drop")

head(gene_variant_summary_clean)

selected_genes <- c("CD247", "CD3D", "CD3E", "CORO1A", "IL2RG", "IL7R", "ITPKB", "JAK3", "LAT", 
                    "LCP2", "PAX1", "PTPRC", "ADA", "AK2", "DCLRE1C", "LIG4", "NHEJ1", "NUDCD3", 
                    "PRKDC", "RAC2", "RAG1", "RAG2", "B2M", "BCL10", "CARD11", "CD3G", "CD40", 
                    "CD40LG", "CD8A", "CHUK", "CIITA", "COPG1", "DOCK2", "DOCK8", "FCHO1", 
                    "FOXI3", "ICOS", "ICOSLG", "IKBKB", "IKZF1")

# Plot as a stacked bar chart:
# - x-axis: GeneSymbol (for many genes, consider filtering or faceting)
# - y-axis: Count of variants per gene
# - Fill: "Pathogenic" (red) vs "Other" (blue)

p_stacked <- gene_variant_summary_clean |> 
  filter(GeneSymbol %in% selected_genes) |>
  ggplot( aes(x = as.factor(GeneSymbol), y = variant_count, fill = variant_class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Other" = "blue", "Pathogenic" = "red")) +
  labs(
    x = "GeneSymbol",
    y = "Variant Count",
    fill = "Classification"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p_stacked)

ggsave(patch1, file = paste0(output, "p_gene_summary_patch1.png"), height = 5, width = 6)
ggsave(p1, file = paste0(output, "p_gene_summary_p1.png"), height = 4, width = 6)
ggsave(patch2, file = paste0(output, "p_gene_summary_hist_patch2.png"), height = 5, width = 6)
ggsave(p1_hist, file = paste0(output, "p_gene_summary_p1_hist.png"), height = 4, width = 6)
ggsave(p_stacked, file = paste0(output, "p_gene_summary_stacked.png"), height = 5, width = 6)

ggsave(patch1, file = paste0(clinvar_out, "p_gene_summary_patch1.png"), height = 5, width = 6)
ggsave(p1, file = paste0(clinvar_out, "p_gene_summary_p1.png"), height = 4, width = 6)
ggsave(patch2, file = paste0(clinvar_out, "p_gene_summary_hist_patch2.png"), height = 5, width = 6)
ggsave(p1_hist, file = paste0(clinvar_out, "p_gene_summary_p1_hist.png"), height = 4, width = 6)
ggsave(p_stacked, file = paste0(clinvar_out, "p_gene_summary_stacked.png"), height = 5, width = 6)

saveRDS(gene_variant_summary_clean, file = paste0(output, "gene_variant_summary.Rds"))
saveRDS(gene_variant_summary_clean, file = paste0(clinvar_sum, "gene_variant_summary.Rds"))
