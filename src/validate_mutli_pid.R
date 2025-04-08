# test ----
we are missing gens from the full run.
running v2


library(readr)
library(dplyr)

# Define a simple processing function for each chunk.
# For testing, we simply add a new column to indicate the chunk's starting position.
process_db_chunk <- function(chunk, pos) {
  chunk <- mutate(chunk, chunk_pos = pos)
  return(chunk)
}

# Initialize an empty list to accumulate chunks.
accumulated <- list()

# Define a callback function that appends the processed chunk to the global list.
chunk_callback <- function(x, pos) {
  accumulated[[length(accumulated) + 1]] <<- process_db_chunk(x, pos)
}

# File to read (using a single file for testing).
data_file <- "../output/scid_datasets.tsv"
message("Processing file: ", data_file)

# Read the file in chunks using read_delim_chunked.
read_delim_chunked(
  file = data_file,
  delim = "\t",
  col_names = TRUE,    # Use header from the file.
  chunk_size = 2,      # Set a small chunk size for testing.
  callback = DataFrameCallback$new(chunk_callback)
)

# Combine all accumulated chunks into one data frame.
result <- bind_rows(accumulated)
print(result)










library(tidyr)

# Convert the dataset from wide to long by pivoting the genetic defect columns
df_long <- pivot_longer(
  df,
  cols = c("SCID", "IL2RG", "RAG1", "Artemis"),
  names_to = "genetic_defect",
  values_to = "count"
)

df_long <- df_long |> filter(!Country == "Israel")

library(dplyr)
library(ggplot2)
library(tidyr)

# Summarise counts per Country and genetic defect, and calculate prevalence
df_summary <- df_long %>%
  group_by(Country, genetic_defect, Population) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(prevalence = total_count / Population)

# Create a bar plot to compare the four genetic defects across populations
p_comparison <- ggplot(df_summary, aes(x = Country, y = prevalence, fill = genetic_defect)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = "Prevalence of Genetic Defects Across Populations",
       x = "Country",
       y = "Prevalence (Count / Population)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

print(p_comparison)
# ggsave("../images/genetic_defects_comparison.png", plot = p_comparison, width = 8, height = 6)






# x ----

library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(stringr)
library(patchwork)
library(ggplot2);theme_set(theme_bw())

PANEL <- 398

# PanelAppRex ----
# Rds format
path_data <- "~/web/PanelAppRex/data"
path_PanelAppData_genes_combined_Rds <- paste0(path_data, "/path_PanelAppData_genes_combined_Rds")
df_core <- readRDS(file= path_PanelAppData_genes_combined_Rds)
colnames(df_core)[colnames(df_core) == 'entity_name'] <- 'Genetic defect'  
df_core <- df_core |> filter(panel_id == PANEL)

# VarRiskEst ----
# source("panelapprex_import.R")

# Load gene tally and variant data for the specified panel using the PANEL variable
varRisEst_gene <- readRDS(file = paste0("~/web/var_risk_est/data/full_run/output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_tally.Rds"))
varRisEst_var  <- readRDS(file = paste0("~/web/var_risk_est/data/full_run/output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_variants.Rds"))



varRisEst_large  <- readRDS(file = paste0("~/web/var_risk_est/data/full_run/var_risk_est_large/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_variants_large.Rds"))




score_map <- c(
  "Pathogenic" = 5,
  "Likely pathogenic" = 4,
  "Pathogenic, low penetrance" = 3,
  "likely pathogenic, low penetrance" = 3,
  "Conflicting classifications of pathogenicity" = 2,
  "risk factor" = 1,
  "association" = 1,
  "likely risk allele" = 1,
  "drug response" = 0,
  "Uncertain significance" = 0,
  "no classification for the single variant" = 0,
  "no classifications from unflagged records" = 0,
  "Affects" = 0,
  "other" = 0,
  "not provided" = 0,
  "uncertain risk allele" = 0,
  "protective" = -3,
  "Likely benign" = -4,
  "Benign" = -5
)

score_df <- tibble(
  classification = names(score_map),
  score = as.numeric(score_map)
) %>%
  # mutate(classification_wrapped = str_wrap(classification, width = 20))
  mutate(
    classification_wrapped = str_wrap(
      str_trunc(classification, width = 37, side = "right", ellipsis = "..."), 
      width = 20
    )
  )

score_df$score <- score_df$score + 0.05
p_score <- ggplot(score_df, aes(x = score, y = reorder(classification_wrapped, score), fill = score )) +
  geom_col(color = "black") +
  scale_fill_gradient2(low = "navy", mid = "lightblue", high = "red", midpoint = 0) +
  labs(
    x = "Score",
    y = "Classification"
  ) +
  guides(fill = "none") 

p_score

get_score <- function(clinvar) {
  # standardise the term: replace underscores with spaces
  clinvar <- str_replace_all(clinvar, "_", " ")
  # split on "/" or "|" (or both)
  terms <- str_split(clinvar, "[/|]", simplify = TRUE)
  terms <- str_trim(terms)
  # get scores for each term; if term not found, assume 0
  scores <- sapply(terms, function(x) if (x %in% names(score_map)) score_map[[x]] else 0)
  mean(scores)
}

# Calculate a score for each row
varRisEst_gene_scored <- varRisEst_gene %>%
  mutate(score = map_dbl(clinvar_clnsig, get_score))

# Compute total count per gene and reorder gene names by descending total count
gene_totals <- varRisEst_gene_scored %>%
  group_by(genename) %>%
  summarise(total_count = sum(count))

varRisEst_gene_scored <- varRisEst_gene_scored %>%
  left_join(gene_totals, by = "genename") %>%
  mutate(genename = reorder(genename, -total_count))

# Order stacking within each gene by increasing score
varRisEst_gene_scored <- varRisEst_gene_scored %>%
  group_by(genename) %>%
  arrange(score) %>%
  mutate(score_factor = factor(score, levels = unique(score))) %>%
  ungroup()

# Create horizontal stacked bar chart with segments ordered by score
p1 <- ggplot(varRisEst_gene_scored, aes(x = count, y = genename, fill = score, ,group = score_factor)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_gradient2(low = "navy", mid = "lightblue", high = "red", midpoint = 0) +
  scale_y_discrete(breaks = function(x) x[seq(1, length(x), 20)]) +
  labs(
    # x = "Classifications per gene (count)",
    y = "Gene name\n(every 20th)",
    fill = "Score"
  )

p1

# Create a numeric index for gene names
# varRisEst_gene_scored <- varRisEst_gene_scored %>%
# mutate(gene_index = as.numeric(genename))

p2 <- ggplot(varRisEst_gene_scored, aes(x = count, y = genename, fill = score, group = score_factor)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_gradient2(low = "navy", mid = "lightblue", high = "red", midpoint = 0) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_discrete(breaks = function(x) x[seq(1, length(x), 20)]) +
  labs(
    # x = " Classifications per gene (%)",
    y = "Gene name\n(every 20th)",
    fill = "Score"
  )

p1 + p2 +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') + 
  plot_layout(axis_titles = "collect")


# Filter for rows with score > 0, then tally counts per gene
varRisEst_gene_scored_positive <- varRisEst_gene_scored %>%
  filter(score > 0) %>%
  group_by(genename) %>%
  summarise(score_positive_total = sum(count), .groups = "drop") %>%
  arrange(desc(score_positive_total))
# slice_head(n = 15)

# Filter for rows with score > 0, then tally counts per gene
top_ranks <- varRisEst_gene_scored_positive %>%
  arrange(desc(score_positive_total)) %>%
  slice_head(n = 15)

# Filter the data to include only the top 10 genes and set the factor levels
top_ranks <- varRisEst_gene_scored %>%
  semi_join(top_ranks, by = "genename") %>%
  mutate(genename = factor(genename, levels = top_ranks$genename))

# Plot horizontal stacked bar chart for the top 10 genes
p3 <- ggplot(top_ranks, aes(x = count, y = genename, fill = score, group = score_factor)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_gradient2(low = "navy", mid = "lightblue", high = "red", midpoint = 0) +
  labs(
    # x = "Classifications per gene (count)",
    y = "Gene name\n(Top 15 pathogenic count)",
    fill = "Score"
  )

p3

# Complex layouts can be created with the `design` argument
design <- "
  123
  144
"
patch1 <- 
  p_score + p1 + p2 + p3 + plot_layout(design = design)  +
  plot_annotation(tag_levels = 'A',) +
  plot_layout(guides = 'collect', axis_titles = "collect")

ggsave(patch1, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PANEL, "_p_varRisEst_summary_scores.pdf"), height = 6, width = 10)

names(varRisEst_gene_scored)

gene_summary <- varRisEst_gene_scored %>%
  group_by(genename) %>%
  summarise(
    # 0-1,2-3,4,5,
    score5 = sum(if_else(score == 5, count, 0L)),
    score4 = sum(if_else(score >=3 & score <= 4, count, 0L)),
    score2 = sum(if_else(score >= 0 & score <= 2, count, 0L)),
    score0 = sum(if_else(score >= -5 & score <= -1, count, 0L)),
  ) %>%
  ungroup() %>%
  mutate(VariantCounts = paste(score5, score4, score2,score0, sep = " / "))


max_pathogenic <- max(gene_summary$Pathogenic, na.rm = TRUE)
max_total <- max(gene_summary$Pathogenic + gene_summary$Other, na.rm = TRUE)

varRisEst_summary <- merge(gene_summary, varRisEst_gene_scored_positive)
names(varRisEst_summary)
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'genename'] <- 'Genetic defect'
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'score5'] <- 'score5.VRE'
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'score4'] <- 'score4.VRE'
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'score2'] <- 'score2.VRE'
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'score0'] <- 'score0.VRE'
colnames(varRisEst_summary)[colnames(varRisEst_summary) == 'VariantCounts'] <- 'VariantCounts.VRE'

# Calculate a score for each row
varRisEst_var_scored <- varRisEst_var %>%
  mutate(score = map_dbl(clinvar_clnsig, get_score))

varRisEst_var_path <- varRisEst_var_scored |> filter(score >= 4)

# Calculate min, max, Q1, Q3 for probabilities
varRisEst_var_path_summary <- varRisEst_var_path %>%
  ungroup() %>%
  group_by(genename) %>%
  summarise(probabilities = list(occurrence_prob), .groups = "drop") %>%
  mutate(
    min_prob = sapply(probabilities, min, na.rm = TRUE),
    q1_prob = sapply(probabilities, function(x) quantile(x, 0.25, na.rm = TRUE)),
    median_prob = sapply(probabilities, median, na.rm = TRUE),
    q3_prob = sapply(probabilities, function(x) quantile(x, 0.75, na.rm = TRUE)),
    max_prob = sapply(probabilities, max, na.rm = TRUE)
  )

colnames(varRisEst_var_path_summary)[colnames(varRisEst_var_path_summary) == 'genename'] <- 'Genetic defect'# 

# Merge summary into the main dataframe and handle NA values
# df <- merge(df, varRisEst_var_path_summary, by = "Genetic defect", all.x = TRUE)
df <- merge(df_core, varRisEst_var_path_summary, by= "Genetic defect", all.x = TRUE )


# df1 <- merge(df, varRisEst_var_path_summary)



names(df)

library(tidyr)

df <- df %>%
  mutate(
    probabilities = map(probabilities, ~ if(all(is.na(.))) rep(0, 5) else replace_na(., 0)),
    min_prob = replace_na(min_prob, 0),
    q1_prob = replace_na(q1_prob, 0),
    median_prob = replace_na(median_prob, 0),
    q3_prob = replace_na(q3_prob, 0),
    max_prob = replace_na(max_prob, 0)
  )



