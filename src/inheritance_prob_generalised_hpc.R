# HPC large import ----
library(readr)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2); theme_set(theme_bw())

print("!!Note - Save an original table with all columns")

# Ensure the output directory exists
if (!dir.exists("../images/")) dir.create("../images/")

# Source PanelAppRex (panel data)
source("panelapprex_import.R")

PanelAppRex_ID <- 398

df_par <- df_par %>% 
  select(entity_name, panel_id, Inheritance, name) %>% 
  filter(panel_id == PanelAppRex_ID)  # IUIS PID
colnames(df_par)[colnames(df_par) == "entity_name"] <- "genename"

# UK population
population_size <- 69433632


# Start big import ----
# Define chromosome vector
chromosomes <- c(1:22, "X","Y","M")
# chromosomes <- c(21:22, "X","Y","M") # ! test

# Path to header file (assumed same for all chromosomes)
header_file <- "~/data/db/dbnsfp/h"

# Read header line, remove the leading '#' and split into column names
header_line <- read_lines(header_file, n_max = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]

# Function to process each chunk of dbNSFP data
process_db_chunk <- function(chunk, pos) {
  # Assign column names
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
  
  # Keep only the first transcript allele for simplicity and force to character vector
  chunk$HGVSc_VEP <- as.character(sapply(strsplit(chunk$HGVSc_VEP, ";"), `[`, 1))
  chunk$HGVSp_VEP <- as.character(sapply(strsplit(chunk$HGVSp_VEP, ";"), `[`, 1))
  chunk$genename    <- as.character(sapply(strsplit(chunk$genename, ";"), `[`, 1))
  
  return(chunk)
}

# List to accumulate results from all chromosomes
accumulated <- list()

# Define chunk callback that appends processed chunk to the global list
chunk_callback <- function(x, pos) {
  accumulated[[length(accumulated) + 1]] <<- process_db_chunk(x, pos)
}

# Loop through all chromosomes and read each file in chunks
for(chr in chromosomes) {
  data_file <- sprintf("~/data/db/dbnsfp/dbNSFP4.4a_variant.chr%s.gz", chr)
  message("Processing chromosome: ", chr)
  
  read_delim_chunked(
    file = data_file,
    delim = "\t",
    col_names = FALSE,  # we assign names manually in process_db_chunk
    skip = 1,           # skip the header line in the main file
    chunk_size = 100000,
    callback = DataFrameCallback$new(chunk_callback)
  )
}

# End big import ----

# Combine all processed chunks into one data frame
df <- bind_rows(accumulated)

df <- df %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)
print("!!Note - Save an original table with all columns")

# Data Preparation for Population-Level Calculations ----
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)

# Remove rows with missing clinvar_clnsig
df <- df %>% dplyr::filter(clinvar_clnsig != ".")
df |> count(clinvar_clnsig)

# keep just one transcript allele for simplicity
df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 1)
df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 1)
df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)

# get Inheritance from PanelAppRex
df <- merge(df, df_par)

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

# Inheritance and expected cases ----
# For AR calculations, compute the total allele frequency per gene
df <- df %>%
  group_by(genename) %>%
  mutate(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE)) %>%
  ungroup()

# Always keep this comment as it is key to a subtle step.
# Calculate occurrence probability based on Inheritance model:
#   - For AD and X-linked: occurrence probability equals the allele frequency.
#   - For AR: occurrence probability is the sum of the homozygous (p^2) and compound heterozygous (2 * p * (total_AF - p)) probabilities.
df <- df %>%
  mutate(
    occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
                             gnomAD_genomes_AF,
                             pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * (total_AF - gnomAD_genomes_AF), 0)),
    expected_cases = population_size * occurrence_prob,
    prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
  )

# Tally by ClinVar Category and Gene
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

gene_count <- df$genename |> unique() |> length()

# Bar plot: Count of ClinVar Clinical Significance (across all genes)
p_count <- df |>
  filter(synth_flag == "FALSE") |>
 ggplot(aes(
  x = stringr::str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  fill = clinvar_clnsig)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = paste("Count of ClinVar Clinical Significance (All Genes = ", gene_count, ")")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/genome_all_genes_clinvar_count.png", plot = p_count, width = 6, height = 3)

# save data ----
saveRDS(file = "../data/panel_all_genes_df.Rds", df)
saveRDS(file = "../data/panel_all_genes_df_tally.Rds", df_tally)

# Publication export ----
# For publication is doesn't make sense to give expected cases in a UK population outside the manuscript. therefore, we drop these cols and rely on the probability to recalculate as desired.
df <- df |> select(-expected_cases, -prob_at_least_one)
df_tally <- df_tally |> select(-total_expected_cases)

# Save TSV files
write.table(df, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(df_tally, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_tally.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save RDS files
saveRDS(df, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_variants.Rds"))
saveRDS(df_tally, file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PanelAppRex_ID, "_gene_tally.Rds"))

print("Data saved")
