library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(stringr)
library(patchwork)
library(ggplot2);theme_set(theme_bw())

print("PID is in ~1 per 1,000 births")
print("SCID is in ~1 per 100,000 births")
print("For a patient with a PID-like phenotype we provide the probabilities of observing the causal genetic determinant")
print("In these SCID we should take the observed number and x100 to give the observed count if we had 1,000,000 cases, rather than all PID")

print("In a global populaiton of 1,000,000/1000 (1 PID per 1000) we have 10,000 PID.")
print("In a global populaiton of 1,000,000/100,000 (1 SCID per 100000) we have 100 SCID.")

print("In 1,000,000 PID we have 100,000 SCID.")
print("In 1,000,000 PID we have 10,000 per SCID gene. (We have around 10 genes)")
print("So for each gene we expect up to 1000.")

# New ----

# File to read (using a single file for testing).
data_file <- "../output/scid_datasets.tsv"
df_scid <- read.table(data_file, header = TRUE)

library(tidyr)

# Convert the dataset from wide to long by pivoting the genetic defect columns
df_scid <- pivot_longer(
  df_scid,
  cols = c( "IL2RG", "RAG1", "DCLRE1C"),
  names_to = "genetic_defect",
  values_to = "count"
)

# df_scid <- df_scid |> filter(!Country == "Israel")

head(df_scid)

library(dplyr)
library(ggplot2)
library(tidyr)

# Summarise counts per Country and genetic defect, and calculate prevalence
df_summary <- df_scid %>%
  group_by(Country, genetic_defect, Population) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(prevalence = total_count / Population)

p_comparison <- ggplot(df_summary, aes(x = Country, y = prevalence, fill = genetic_defect)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_text(aes(label = total_count, y = prevalence + 1e-6),
            position = position_dodge(width = 0.9),
            size = 4) +
  labs(x = "Country",
       y = "Prevalence (Count / Population)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

print(p_comparison)

ggsave("../images/validation_studies_scid_gene_comparison.pdf", plot = p_comparison, width = 6, height = 4)


df_summary

# test ----
library(dplyr)
library(ggplot2)
library(binom)
library(purrr)
library(tidyr)

library(dplyr)
library(ggplot2)
library(binom)


# prep the count for SCID only rather than all PID.
df_scid$count <- (df_scid$count*100)

df_il2rg_all <- df_scid %>% 
  filter(genetic_defect == "IL2RG") %>% 
  mutate(rate_per_million = (count / Population) * 1e6)

p_lines <- ggplot() +
  geom_vline(data = df_il2rg_all, 
             aes(xintercept = rate_per_million, colour = Country), 
             linetype = "dotted", size = 1) +
  geom_text(data = df_il2rg_all, 
            aes(x = rate_per_million, y = 0, label = paste(Country, ":", count)),
            angle = 90, vjust = -0.5, hjust = 0, size = 4, show.legend = FALSE) +
  labs(title = "IL2RG-related SCID Cases per 1,000,000 Population",
       x = "Cases per 1,000,000",
       y = "") 
  # theme_bw() #+
  # # theme(axis.text.y = element_blank(),
  #       axis.ticks.y = element_blank())


p_lines
# ggsave("../images/il2rg_scid_cases_per_million.png", plot = p_lines, width = 9, height = 3)

# varestrisk ----
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
varRisEst_gene <- readRDS(file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_tally.Rds"))
varRisEst_var  <- readRDS(file = paste0("../output/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_variants.Rds"))

# varRisEst_large  <- readRDS(file = paste0("~/web/var_risk_est/data/full_run/var_risk_est_large/VarRiskEst_PanelAppRex_ID_", PANEL, "_gene_variants_large.Rds"))


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



# 
# library(dplyr)
# library(ggplot2)
# 
df_scid$genetic_defect |> unique()
# [1] "IL2RG"   "RAG1"    "DCLRE1C"
# 
# gene <- "IL2RG"
# 

# figure ----
library(dplyr)
library(ggplot2)
library(patchwork)





plot_gene_scid <- function(gene) {
  # Filter dataset and calculate adjusted counts and rate per million using expected SCID incidence (1 per 100,000)
  df_gene <- df_scid %>% 
    filter(genetic_defect == gene) %>% 
    mutate(expected_SCID = Population / 100000,
           ratio = if_else(SCID > 0, expected_SCID / SCID, NA_real_),
           adjusted_count = if_else(SCID > 0, count * ratio, NA_real_),
           rate_per_million = (adjusted_count / Population) * 1e6)
  
  # Get overall predicted values from varRisEst_gene_scored for the two score groups
  overall_pathogenic <- varRisEst_gene_scored %>% 
    filter(genename == gene, score == 5) %>% 
    pull(overall_prob) * 1e6
  
  overall_likely <- varRisEst_gene_scored %>% 
    filter(genename == gene, score > 2) %>% 
    summarise(total_prob = sum(overall_prob)) %>% 
    pull(total_prob) * 1e6
  
  # Summarise to get one vline per country for the legend (avoiding duplicates)
  df_countries <- df_gene %>% 
    group_by(Country) %>% 
    summarise(rate_per_million = mean(rate_per_million), .groups = "drop")
  
  # Base plot with vlines for each country and overall predictions 
  p <- ggplot() +
    # Country-specific vlines (mapped to country names)
    geom_vline(data = df_countries,
               aes(xintercept = rate_per_million, colour = Country),
               linetype = "dotted", size = 1) +
    # Overall predicted vlines mapped to dummy groups to force them into the legend
    geom_vline(aes(xintercept = overall_pathogenic, colour = "Predicted median"),
               linetype = "solid", size = 1) +
    geom_vline(aes(xintercept = overall_likely, colour = "Predicted incl LP"),
               linetype = "solid", size = 1) +
    annotate("text", x = overall_pathogenic, y = Inf,
             label = paste("Predicted\ntotal:", round(overall_pathogenic, 2)),
             vjust = 2, colour = "darkgreen") +
    annotate("text", x = overall_likely, y = Inf,
             label = paste("Predicted including\nLP:", round(overall_likely, 2)),
             vjust = 4, colour = "navy") +
    labs(title = paste0(gene, " - related SCID"),
         x = "Cases SCID per 1,000,000 PID",
         y = "") 
    # theme_bw() #+
    # theme(axis.text.y = element_blank(),
          # axis.ticks.y = element_blank())
  
  # Observed distribution based on adjusted rate per million (normal approximation)
  nsamples <- 10000
  mean_rate <- mean(df_gene$rate_per_million, na.rm = TRUE)
  sd_rate <- sd(df_gene$rate_per_million, na.rm = TRUE)
  norm_samples <- rnorm(nsamples, mean = mean_rate, sd = sd_rate)
  median_norm <- median(norm_samples)
  ci_norm <- quantile(norm_samples, probs = c(0.025, 0.975))
  
  p <- p +
    geom_density(data = data.frame(norm_samples = norm_samples),
                 aes(x = norm_samples, y = ..scaled..),
                 fill = "orange", alpha = 0.5) +
    geom_vline(xintercept = median_norm, linetype = "dotted", colour = "black", size = 1) +
    geom_vline(xintercept = ci_norm[1], linetype = "dotted", colour = "black", size = 1) +
    geom_vline(xintercept = ci_norm[2], linetype = "dotted", colour = "black", size = 1) +
    annotate("text", x = median_norm, y = Inf,
             label = paste("normal median:", round(median_norm, 2)),
             colour = "red", vjust = 0)
  
  # Predicted distributions for each score group and combine them into a single set of samples
  set.seed(666)
  nsamples_exp <- 10000
  
  sd_pathogenic <- overall_pathogenic * 0.1
  expected_pathogenic_samples <- rnorm(nsamples_exp, mean = overall_pathogenic, sd = sd_pathogenic)
  
  sd_likely <- overall_likely * 0.1
  expected_likely_samples <- rnorm(nsamples_exp, mean = overall_likely, sd = sd_likely)
  
  overall_expected_samples <- c(expected_pathogenic_samples, expected_likely_samples)
  median_overall_expected <- median(overall_expected_samples)
  ci_overall_expected <- quantile(overall_expected_samples, probs = c(0.025, 0.975))
  
  p <- p +
    geom_density(data = data.frame(overall_expected_samples = overall_expected_samples),
                 aes(x = overall_expected_samples, y = ..scaled..),
                 fill = "lightgreen", alpha = 0.5)
  
  # Annotate the normal distribution median if it is greater than zero
  if(median_overall_expected > 0) {
    p <- p + 
      geom_vline(xintercept = median_overall_expected, linetype = "dotdash", colour = "darkgreen", size = 1) +
      annotate("text", x = median_overall_expected, y = Inf,
               label = paste("Median:", round(median_overall_expected, 2)),
               colour = "darkgreen", vjust = 8)
  }
  
  # Build a custom colour mapping that includes each unique country plus the two predicted groups
  country_levels <- sort(unique(df_countries$Country))
  legend_levels <- c(country_levels, "Predicted median", "Predicted incl LP")
  # Set all countries to orange and the two predicted groups to darkgreen and navy respectively
  legend_colors <- c(rep("orange", length(country_levels)), "darkgreen", "navy")
  
  p <- p + scale_colour_manual(name = "Legend",
                               breaks = legend_levels,
                               values = setNames(legend_colors, legend_levels))
  
  
  
  # Build a custom colour mapping that includes each unique country plus the two predicted groups
  country_levels <- sort(unique(df_countries$Country))
  legend_levels <- c(country_levels, "Predicted median", "Predicted incl LP")
  # Set all countries to orange and the two predicted groups to darkgreen and navy respectively
  legend_colors <- c(rep("orange", length(country_levels)), "darkgreen", "navy")
  
  p <- p + scale_colour_manual(name = "Legend",
                               breaks = legend_levels,
                               values = setNames(legend_colors, legend_levels),
                               guide = guide_legend(
                                 override.aes = list(linetype = "solid",
                                                     shape = NA,
                                                     size = 5,
                                                     keywidth = 2)
                               ))
  
  return(p)
}

# Generate plots for each gene
p_il2rg   <- plot_gene_scid("IL2RG")
p_rag1    <- plot_gene_scid("RAG1")
p_dclre1c <- plot_gene_scid("DCLRE1C")

# Combine plots using patchwork (stack vertically)
combined_plot <- p_il2rg / p_rag1 / p_dclre1c +
  plot_layout(guides = 'collect', axis = "collect") +
  plot_annotation(tag_levels = 'A')

combined_plot

# Save the combined plot
ggsave("../images/validation_studies_scid_combined_plot.pdf", plot = combined_plot, width = 8, height = 7)
