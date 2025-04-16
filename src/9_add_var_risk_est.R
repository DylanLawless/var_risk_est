library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(stringr)
library(ggplot2);theme_set(theme_bw())

varRisEst_gene <- readRDS(file= "../output/VarRiskEst_PanelAppRex_ID_398_gene_tally.Rds") # gene level 
varRisEst_var <- readRDS(file= "../output/VarRiskEst_PanelAppRex_ID_398_gene_variants.Rds") # all variants
# test2 <- readRDS(file= "../data/var_risk_est/panel_all_genes_df.Rds") # includes UK pop: "expected_cases", "prob_at_least_one"
# test3 <- readRDS(file= "../data/var_risk_est/panel_all_genes_df_tally.Rds") # includes UK pop: "expected_cases", "prob_at_least_one"

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
  scale_y_discrete(breaks = function(x) x[seq(1, length(x), 100)]) +
  labs(
    # x = "Classifications per gene (count)",
    y = "Gene name\n(every 100th)",
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
  scale_y_discrete(breaks = function(x) x[seq(1, length(x), 100)]) +
  labs(
    # x = " Classifications per gene (%)",
    y = "Gene name\n(every 100th)",
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

# p_score | ((p1 + p2) / p3) +
  # plot_annotation(tag_levels = 'A',) +
  # plot_layout(guides = 'collect',  widths = c(.3,3), axis_titles = "collect")

# Complex layouts can be created with the `design` argument
design <- "
  123
  144
"
patch1 <- 
  p_score + p1 + p2 + p3 + plot_layout(design = design)  +
  plot_annotation(tag_levels = 'A',) +
  plot_layout(guides = 'collect', axis_titles = "collect")



ggsave(patch1, file = "../output/p_varRisEst_summary_scores.pdf", height = 6, width = 10)
ggsave(patch1, file = "~/web/var_risk_est/images/p_varRisEst_summary_scores.pdf", height = 6, width = 10)


  
print("NEXT LOOK AT THE OVERLL_PROB score and start fresh - maybe the same plots")
print("THEN LOOK AT CLUSTERING with string db and color by score - quantify which cluster are the biggest threat. ")


# Aggregate counts per genename

# gene_summary <- varRisEst_gene_scored %>%
#   group_by(genename) %>%
#   summarise(
#     Pathogenic = sum(if_else(score > 0, count, 0L)),
#     Other = sum(if_else(score <= 0, count, 0L))
#   ) %>%
#   ungroup() %>%
#   mutate(VariantCounts = paste(Pathogenic, Other, sep = " / "))
# 
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

# clinvar_summary <- final_results_gene_scored %>%
#   group_by(GeneSymbol) %>%
#   summarise(
#     score5 = sum(if_else(score == 5, total, 0L)),
#     score4 = sum(if_else(score >= 2 & score <= 4, total, 0L)),
#     score2 = sum(if_else(score >= 0 & score <= 2, total, 0L)),
#     score0 = sum(if_else(score >= -5 & score <= -1, total, 0L)),
#   ) %>%
#   ungroup() %>%
#   mutate(VariantCounts = paste(score5, score4, score2,score0, sep = " / "))

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

df <- merge(df, varRisEst_summary, by= "Genetic defect", all.x = TRUE )



# # variant probability sparkline ----# 
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
df <- merge(df, varRisEst_var_path_summary, by = "Genetic defect", all.x = TRUE)

df <- df %>%
  mutate(
    probabilities = map(probabilities, ~ if(all(is.na(.))) rep(0, 5) else replace_na(., 0)),
    min_prob = replace_na(min_prob, 0),
    q1_prob = replace_na(q1_prob, 0),
    median_prob = replace_na(median_prob, 0),
    q3_prob = replace_na(q3_prob, 0),
    max_prob = replace_na(max_prob, 0)
  )


