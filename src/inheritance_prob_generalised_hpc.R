# HPC large import ----
library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2); theme_set(theme_bw())

print("Version 2 run")
setwd("~/var_risk_est/src")
large_data <- "~/data/var_risk_est_large"

# Retrieve PanelAppRex_ID from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
PanelAppRex_ID <- as.numeric(args[1])
print(paste("Now running panel ID:", PanelAppRex_ID))

# Ensure the output directory exists
if (!dir.exists("../images/")) dir.create("../images/")

# Source PanelAppRex (panel data)
source("panelapprex_import.R")

# Subset panel data for the current panel and rename column
df_par <- df_par %>% 
  select(entity_name, panel_id, Inheritance, name) %>% 
  filter(panel_id == PanelAppRex_ID)
colnames(df_par)[colnames(df_par) == "entity_name"] <- "genename"

# Set population size (e.g., 100 million)
population_size <- 100000000

# Define chromosomes to loop over
chromosomes <- c(1:22, "X", "Y", "M")

# Path to header file (assumed same for all chromosomes)
header_file <- "~/data/db/dbnsfp/h"
header_line <- read_lines(header_file, n_max = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]

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
  
  # Ensure genename is a character vector
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
  data_file <- sprintf("~/data/db/dbnsfp/dbNSFP4.4a_variant.chr%s.gz", chr)
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
saveRDS(df, file = paste0(large_data, "/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants_large.Rds"))

# The df is already cleaned from the chunk processing; no need to repeat splitting/transformation.
# Data Preparation for Population-Level Calculations ----

print(paste("dim df:", dim(df)))

# Merge with panel data to get Inheritance information
df <- merge(df, df_par)
print(paste("dim df after PAR merge:", dim(df)))

# Apply minimal risk for variants with zero observed allele frequency.
# For each gene, if gnomAD_genomes_AF is 0, assign 1/(max_an + 1).
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

print(paste("dim df after synth:", dim(df)))

# For each (genename, clinvar_clnsig) group, if non-synthetic rows exist, keep only those;
# otherwise, keep one synthetic row.
df <- df %>%
  group_by(genename, clinvar_clnsig) %>%
  group_modify(~ {
    if (any(!.x$synth_flag)) {
      .x %>% filter(!synth_flag)
    } else {
      .x %>% slice(1)
    }
  }) %>%
  ungroup()

print(paste("dim df synth filt:", dim(df)))

# For AR inheritance, calculate total allele frequency per gene.
df <- df %>%
  group_by(genename) %>%
  mutate(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE)) %>%
  ungroup()

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
## Version 1 ----
# df <- df %>%
#   mutate(
#     occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
#                              gnomAD_genomes_AF,
#                              pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * (total_AF - gnomAD_genomes_AF), 0)),
#     expected_cases = population_size * occurrence_prob,
#     prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
#   )
#
## Version 2: We update occurrence_prob per variant accordingly. ----
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

# Tally counts and expected cases by gene and ClinVar classification
clinvar_levels <- unique(df$clinvar_clnsig)

df_tally <- df %>%
  group_by(genename, clinvar_clnsig) %>%
  summarise(
    count = n(),
    total_expected_cases = sum(expected_cases),
    overall_prob = 1 - prod(1 - occurrence_prob),
    .groups = "drop"
  ) %>%
  tidyr::complete(genename, clinvar_clnsig = clinvar_levels,
                  fill = list(count = 0, total_expected_cases = 0, overall_prob = 0))

gene_count <- length(unique(df$genename))

print(paste("Gene count:", gene_count))

# Plot: Count of ClinVar Clinical Significance across all genes
p_count <- df %>%
  filter(synth_flag == FALSE) %>%
  ggplot(aes(x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
             fill = clinvar_clnsig)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = paste("Count of ClinVar Clinical Significance (All Genes =", gene_count, ")")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave(paste0("../images/genome_all_genes_clinvar_count_ID_", PanelAppRex_ID, ".png"), 
       plot = p_count, width = 6, height = 3)

# Save TSV and RDS outputs
write.table(df, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(df_tally, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_tally.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

saveRDS(df, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants.Rds"))
saveRDS(df_tally, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_tally.Rds"))

print("Data saved")
print(paste("Complete in R for panel ID:", PanelAppRex_ID))

print(paste("df dim:", dim(df)))
print(paste("df_tally dim:", dim(df_tally)))
