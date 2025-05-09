# HPC large import ----
library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2); theme_set(theme_bw())

path_data <- "../data/"

# Ensure the output directory exists
if (!dir.exists("../images/")) dir.create("../images/")

# source PanelAppRex ----
source("panelapprex_import.R")

PanelAppRex_ID <- as.numeric(398) # args[1] on HPC from launch script
print(paste("Now running panel ID:", PanelAppRex_ID))

# Select our panel data ----
df_par <- df_par |> dplyr::select(entity_name, panel_id, Inheritance, name)
df_par <- df_par |> filter(panel_id == PanelAppRex_ID) # IUIS PID
colnames(df_par)[colnames(df_par) == 'entity_name'] <- 'genename'

# UK population
population_size <- 69433632
# population_size <- 83702  # births in 2023

# Define chromosomes to loop over
chromosomes <- c("")

# Start import data ----
header_line <- readLines("../data/nfkb1_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

# Function to process each chunk of dbNSFP data.
# (Note: All cleaning steps below will only be performed once per chunk.)
process_db_chunk <- function(chunk, pos) {
  # Assign column names from header_fields
  colnames(chunk) <- header_fields
  
  # Select relevant columns
  chunk <- chunk %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, 
                            clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)
  
  # Replace missing values and convert allele frequency and number to numeric
  chunk$gnomAD_genomes_AF[chunk$gnomAD_genomes_AF == "."] <- 0
  chunk$gnomAD_genomes_AN[chunk$gnomAD_genomes_AN == "."] <- 0
  chunk$gnomAD_genomes_AF <- as.numeric(chunk$gnomAD_genomes_AF)
  
  # Remove rows with missing clinvar_clnsig
  chunk <- chunk %>% filter(clinvar_clnsig != ".")
  
  # Copy original columns to new ones with "_all" appended
  chunk$HGVSc_VEP_all <- chunk$HGVSc_VEP
  chunk$HGVSp_VEP_all <- chunk$HGVSp_VEP
  chunk$genename_all  <- chunk$genename
  
  # Keep only the first transcript allele for simplicity
  chunk$HGVSc_VEP <- sapply(strsplit(as.character(chunk$HGVSc_VEP), ";"), `[`, 1)
  chunk$HGVSp_VEP <- sapply(strsplit(as.character(chunk$HGVSp_VEP), ";"), `[`, 1)
  chunk$genename   <- sapply(strsplit(as.character(chunk$genename), ";"), `[`, 1)
  
  # Ensure character vector
  chunk$HGVSc_VEP <- as.character(chunk$HGVSc_VEP)
  chunk$HGVSp_VEP <- as.character(chunk$HGVSp_VEP)
  chunk$genename <- as.character(chunk$genename)
  
  return(chunk)
}

# List to accumulate results from all chromosomes
accumulated <- list()

# Define chunk callback that appends processed chunks to the global list
chunk_callback <- function(x, pos) {
  accumulated[[length(accumulated) + 1]] <<- process_db_chunk(x, pos)
}

print(paste("Running on all chromosomes for panel ID:", PanelAppRex_ID))
for(chr in chromosomes) {
  data_file <- sprintf("../data/nfkb1%s", chr) # mimics the chromosome-split version on HPC
  message("Processing chromosome: ", chr)
  
  read_delim_chunked(
    file = data_file,
    delim = "\t",
    col_names = FALSE,  # header will be set in process_db_chunk
    skip = 1,           # skip header in the main file
    chunk_size = 100000,
    callback = DataFrameCallback$new(chunk_callback)
  )
}

# Combine all processed chunks into one data frame
df <- bind_rows(accumulated)
print(paste("Binding all rows complete for panel ID:", PanelAppRex_ID))

# Save the original processed table for future reference
# saveRDS(df, file = paste0(large_data, "/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants_large.Rds"))

# simple version gene 1
# df <- read.table("../data/nfkb1", 
#                  sep = "\t",
#                  header = FALSE, 
#                  stringsAsFactors = FALSE, 
#                  fill = TRUE)
# colnames(df) <- header_fields
df_gene1 <- df 

# simple version gene 2
# df <- read.table("../data/cftr", 
#                  sep = "\t",
#                  header = FALSE, 
#                  stringsAsFactors = FALSE, 
#                  fill = TRUE)

print(paste("Running on all chromosomes for panel ID:", PanelAppRex_ID))
for(chr in chromosomes) {
  data_file <- sprintf("../data/cftr%s", chr) # mimics the chromosome-split version on HPC
  message("Processing chromosome: ", chr)
  
  read_delim_chunked(
    file = data_file,
    delim = "\t",
    col_names = FALSE,  # header will be set in process_db_chunk
    skip = 1,           # skip header in the main file
    chunk_size = 100000,
    callback = DataFrameCallback$new(chunk_callback)
  )
}

# Combine all processed chunks into one data frame
df <- bind_rows(accumulated)
print(paste("Binding all rows complete for panel ID:", PanelAppRex_ID))

# colnames(df) <- header_fields
df_gene2 <- df

df <- rbind(df_gene1, df_gene2)

# print("Debug check for CFTR size")
# df |> filter(genename == "CFTR") |> dim()


# 2594
# End import data ----

# # simple version ----
# 
# 
# # simple version gene 1
# df <- read.table("../data/nfkb1",
#                  sep = "\t",
#                  header = FALSE,
#                  stringsAsFactors = FALSE,
#                  fill = TRUE)
# colnames(df) <- header_fields
# df_gene1 <- df
# 
# # simple version gene 2
# df <- read.table("../data/cftr",
#                  sep = "\t",
#                  header = FALSE,
#                  stringsAsFactors = FALSE,
#                  fill = TRUE)
# 
# colnames(df) <- header_fields
# df_gene2 <- df
# 
# df <- rbind(df_gene1, df_gene2)

# names(df)

# df <- df %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)
# 
# # Data Preparation for Population-Level Calculations ----
# df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
# df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
# df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
# 
# # Remove rows with missing clinvar_clnsig
# df <- df %>% dplyr::filter(clinvar_clnsig != ".")
# df |> count(clinvar_clnsig)
# 
# # keep just one transcript allele for simplicity
# df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 1)
# df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 1)
# df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)
# # Remove gene filter to process all genes
# # df <- df |> filter(genename == "NFKB1")

# get Inheritance from PanelAppRex
dfx <- merge(df, df_par)
dfx <- merge(df, df_par, by = "genename")

# df$Inheritance <- "AD"  # setting inheritance to AD for demonstration
# head(df)

# If no known variants per clinsig, consider a minimal risk with 1 de novo ----
df <- df %>%
  group_by(genename) %>%
  mutate(
    gnomAD_genomes_AF = as.numeric(gnomAD_genomes_AF),
    gnomAD_genomes_AN = as.numeric(gnomAD_genomes_AN),
    max_an = max(gnomAD_genomes_AN, na.rm = TRUE),
    synth_flag = gnomAD_genomes_AF == 0,
    gnomAD_genomes_AF = ifelse(synth_flag, 1 / (max_an + 1), gnomAD_genomes_AF)
  ) %>%
  ungroup()

# Run this but skip it is special case if we want a variant that is not on gnomad
if (!exists("KEPP_ALL_FOR_VALIDATION_SEARCH") || !KEPP_ALL_FOR_VALIDATION_SEARCH) {
  message("Running de novo estimate grouping and filtering step.")
  df <- df %>%
    group_by(genename, clinvar_clnsig) %>%
    group_modify(~ {
      # Always keep this comment as it is key to a subtle step.
      # For each group defined by genename and clinvar_clnsig:
      #   - If there is at least one row that is not synthetic (synth_flag == FALSE),
      #     then keep only the non-synthetic rows.
      #   - Otherwise, if all rows in the group are synthetic, keep only one synthetic row.
      if (any(!.x$synth_flag)) {
        .x %>% filter(!synth_flag)
      } else {
        .x %>% slice(1)
      }
    }) %>%
    ungroup()
} else {
  message("Skipping de novo estimate grouping and filtering step due to KEPP_ALL_FOR_VALIDATION_SEARCH being TRUE.")
}

names(df)

# Inheritance and expected cases ----
# For AR calculations, compute the total allele frequency per gene
df <- df %>%
  group_by(genename) %>%
  mutate(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE)) %>%
  ungroup()

# # Always keep this comment as it is key to a subtle step.
# # Calculate occurrence probability based on Inheritance model:
# #   - For AD and X-linked: occurrence probability equals the allele frequency.
# #   - For AR: occurrence probability is the sum of the homozygous (p^2) and compound heterozygous (2 * p * (total_AF - p)) probabilities.
# df <- df %>%
#   mutate(
#     occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
#                              gnomAD_genomes_AF,
#                              pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * (total_AF - gnomAD_genomes_AF), 0)),
#     expected_cases = population_size * occurrence_prob,
#     prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
#   )
# 
# head(df)


# Inheritance and expected cases ----
# For AR (autosomal recessive) cases, we first calculate the overall gene-level pathogenic
# allele frequency: total_AF = sum of all gnomAD_genomes_AF values for that gene.
# Under Hardy–Weinberg equilibrium, the probability that an individual is affected (by any
# combination of pathogenic variants in that gene) is (total_AF)^2. This accounts for both:
#   - homozygous variants (p_i^2 for each individual variant with frequency p_i), and 
#   - compound heterozygotes (2 * p_i * p_j for two different variants).
#
# The original per-variant calculation computed for each variant:
#     p_i^2 + 2 * p_i * (total_AF - p_i)
# This expression sums to: 2 * p_i * total_AF - p_i^2.
# If you sum this over all variants, the compound heterozygous events are double counted.
#
# To attribute the overall gene-level AR risk (total_AF^2) to each variant without double counting,
# we allocate the risk in proportion to the variant's allele frequency.
#
# The proportional per-variant risk is calculated as:
#     occurrence_prob for variant i = p_i * total_AF
#
# With this strategy, when you sum occurrence_prob over all variants in a gene, you get:
#     sum(p_i * total_AF) = total_AF * (sum(p_i)) = total_AF^2,
# which exactly equals the overall gene-level probability.
#
# Thus, for AR diseases, each variant's occurrence_prob reflects its share of the overall risk 
# that an affected individual carries that particular variant (either in homozygosity or as one 
# of two variants in a compound heterozygous genotype).
#
# version 2: We update occurrence_prob per variant accordingly.
df <- df %>%
  mutate(
    # For AD (autosomal dominant) and X-linked, the probability remains the allele frequency.
    # For AR, we assign the per-variant risk as the product of its allele frequency and the total
    # allele frequency for that gene (i.e. p_i * total_AF), which partitions the gene-level (total_AF)^2
    # risk among the variants.
    occurrence_prob = ifelse(
      Inheritance %in% c("AD", "X-linked"),
      gnomAD_genomes_AF,
      gnomAD_genomes_AF * total_AF
    ),
    # Multiply the per-variant probability by the population size to estimate the number of expected
    # cases attributable to that variant.
    expected_cases = population_size * occurrence_prob,
    # The probability of observing at least one event (case) in the population, given the occurrence
    # probability per individual.
    prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
  )


# test
df |> filter(genename == "CFTR") |> filter(clinvar_clnsig  == "Pathogenic") |> filter(HGVSp_VEP == "p.Arg36His") |> select(occurrence_prob)
# A tibble: 1 × 1


# > # test recessive calc ----
# >  df_v1 |> filter(genename == "CFTR") |> filter(clinvar_clnsig  == "Pathogenic") |> filter(HGVSp_VEP == "p.Arg36His") |> select(occurrence_prob)
# # A tibble: 1 × 1
# occurrence_prob
# <dbl>
#   1         0.00195
# > 
#   > df_v2 |> filter(genename == "CFTR") |> filter(clinvar_clnsig  == "Pathogenic") |> filter(HGVSp_VEP == "p.Arg36His") |> select(occurrence_prob)
# # A tibble: 1 × 1
# occurrence_prob
# <dbl>
#   1        0.000975

# Tally by ClinVar Category and Gene
clinvar_levels <- unique(df$clinvar_clnsig)

df_tally <- df %>%
  group_by(genename, clinvar_clnsig) %>%
  summarise(
    total_expected_cases = sum(expected_cases),
    overall_prob = 1 - prod(1 - occurrence_prob),
    .groups = "drop"
  ) %>%
  tidyr::complete(genename, clinvar_clnsig = clinvar_levels,
                  fill = list(total_expected_cases = 0, overall_prob = 0))

print(df_tally)

# Bar plot: Count of ClinVar Clinical Significance (across all genes)
p_count <- ggplot(df, aes(
  x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  fill = clinvar_clnsig)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = "Count of ClinVar Clinical Significance (All Genes)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/all_genes_clinvar_count_mini.png", plot = p_count, width = 6, height = 3)


# Scatter Plots for one clinvar category example (Pathogenic) across all genes ----
threshold_AF <- min(df$gnomAD_genomes_AF[df$prob_at_least_one >= 0.999])
threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)

unique_labels <- df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>%
  distinct(gnomAD_genomes_AF, expected_cases)





# Scatter Plots: Expected Cases and Probability vs Allele Frequency ----

# 
# # Compute unique labels per gene for pathogenic entries
# unique_labels <- df %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   distinct(genename, gnomAD_genomes_AF, expected_cases)
# 
# p_scatter1_path <-
#   df %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
#   geom_point(size = 3) +
#   geom_line(aes(group = 1)) +
#   ggrepel::geom_text_repel(data = unique_labels,
#                            aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)),
#                            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
#   labs(x = "Allele Frequency (log scale)",
#        y = "Expected Cases",
#        title = "Expected Cases vs\nAllele Frequency CFTR",
#        subtitle = paste0("Condition: population size ", population_size)) +
#   facet_wrap(~ genename, scales = "free")
# 
# p_scatter1_path
# 
# # Compute threshold AF per gene for pathogenic entries
# thresholds <- df %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   group_by(genename) %>%
#   summarise(threshold_AF = min(gnomAD_genomes_AF[prob_at_least_one >= 0.999], na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(threshold_AF_label = formatC(threshold_AF, format = "f", digits = 6))
# 
# p_scatter2_path <-
#   df %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, colour = clinvar_clnsig)) +
#   geom_point(size = 3) +
#   geom_line(aes(group = 1)) +
#   labs(x = "Allele Frequency (log scale)",
#        y = "Probability of ≥1 Event",
#        title = "Probability of At Least One\nEvent vs Allele Frequency CFTR",
#        subtitle = paste0("Condition: population size ", population_size)) +
#   geom_vline(data = thresholds, aes(xintercept = threshold_AF),
#              linetype = "dotted", colour = "black") +
#   geom_text(data = thresholds,
#             aes(x = threshold_AF, y = 1, label = threshold_AF_label),
#             vjust = -1, hjust = 0.5, colour = "black", size = 3) +
#   facet_wrap(~ genename, scales = "free")
# 
# p_scatter2_path
# 
# # Display the scatter plots.
# p_scatter1_path
# p_scatter2_path
# # 
# # Combine and save scatter plots vertically.
# # p_scatter1_path <- p_scatter1_path + 
#   # theme(axis.title.x = element_blank(),
#         # axis.text.x = element_blank(),
#         # axis.ticks.x = element_blank())
# 
# p_scatter <- p_scatter1_path / p_scatter2_path + 
#   plot_layout(guides = 'collect', axis = "collect")  + 
#   plot_annotation(tag_levels = 'A')
# print(p_scatter)
# 
# # Density histograms for Expected Cases by ClinVar Clinical Significance.
# p_density <- ggplot(df, aes(x = expected_cases, fill = clinvar_clnsig)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~ clinvar_clnsig, scales = "free", ncol = 4) +
#   scale_x_continuous(labels = function(x) format(round(x, 0), big.mark = ",")) +
#   guides(fill = "none") +
#   labs(x = "Expected Cases", 
#        y = "Density",
#        title = "Density of Expected Cases by ClinVar Clinical Significance",
#        subtitle = paste0("Condition: population size ", population_size)) 
# 
# print(p_density)
# 
# p_scatter_dense <- (p_density / (p_scatter1_path + p_scatter2_path)) +
#   plot_layout(widths = c(1, 1), guides = 'collect', axis = "collect") +
#   plot_annotation(tag_levels = 'A')
# print(p_scatter_dense)
# 
# ggsave("../images/cftr_scatterdense_expected_prob.png", plot = p_scatter_dense, width = 10, height = 8)
# 







p_scatter1_path <- 
  df |> filter(clinvar_clnsig == "Pathogenic") |>
  ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, colour = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  ylim(0, 800) +
  geom_text(data = unique_labels, 
            aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF)",
       y = "Expected Cases",
       title = "Expected Cases vs Allele Frequency (Pathogenic)") +
  facet_wrap(~ genename)  # facet by gene

vline_label <- tibble(threshold_AF = threshold_AF, 
                      threshold_AF_label = threshold_AF_label)

p_scatter2_path <- 
  df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>% 
  ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, colour = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF)",
       y = "Probability of ≥1 Case",
       title = "Probability of At Least One Case vs Allele Frequency (Pathogenic)") +
  geom_vline(xintercept = threshold_AF, linetype = "dotted", colour = "black") +
  geom_text(data = vline_label, 
            aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  ylim(0, 1.2) +
  facet_wrap(~ genename)

p_scatter1_path
p_scatter2_path

p_scatter1_path_mod <- p_scatter1_path + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_scatter <- p_scatter1_path_mod / p_scatter2_path + 
  plot_layout(guides = 'collect', axis = "collect") + 
  plot_annotation(tag_levels = 'A')
print(p_scatter)
# ggsave("../images/scatter_expected_prob_all_genes.png", plot = p_scatter, width = 6, height = 8)

p_bar <- ggplot(df_tally, aes(
  x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  ggrepel::geom_text_repel(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")),
                           direction = "y",
                           nudge_x = 0,
                           vjust = -1,
                           size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Example, Total Expected Cases in UK population size (~ 69.4M)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .75)),
                     labels = scales::comma) +
  facet_wrap(~ genename, scales = "free_y") +
  theme(text = element_text(face = "bold"))


p_bar
# ggsave("../images/bar_expected_cases_all_genes.png", plot = p_bar, width = 8, height = 6)

# Bar plot: Overall Probability by ClinVar Category for each gene
p_prob <- ggplot(df_tally, aes(
  x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.01)),
            vjust = -0.5, size = 3.5) +
  # geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.01)), 
  #           hjust = -.2,size = 3.5, angle = 90) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall\nProbability",
       title = "Overall Probability of an Affected Birth by ClinVar Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .5)), labels = scales::percent) +
  facet_wrap(~ genename, scales = "free_y") +
  theme(text = element_text(face = "bold"))

p_prob
# ggsave("../images/bar_overall_probability_all_genes.png", plot = p_prob, width = 8, height = 5)

# Combine bar charts using patchwork and save
# p_bar_mod <- p_bar + theme( axis.title.x = element_blank(),
#                            axis.text.x = element_blank(),
#                            axis.ticks.x = element_blank())
p_prob_mod <- p_prob + theme( axis.title.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())

p_bars <- (p_prob_mod / p_bar) + 
  plot_layout(guides = 'collect', axis = "collect") + 
  plot_annotation(tag_levels = 'A', title = "Recessive and Dominant Disease Genes") +
  theme(text = element_text(face = "bold"))
print(p_bars)
ggsave("../images/all_genes_combined_bar_charts_mini.png", plot = p_bars, width = 12, height = 6)

