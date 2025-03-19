library(ggplot2); theme_set(theme_bw())
library(patchwork)

# Ensure the output directory exists
if(!dir.exists("../images/")) dir.create("../images/")

# UK population
population_size <- 69433632
# population_size <- 83702  # births in 2023

# DBNSFP Data Import and Filtering ----
library(dplyr)

header_line <- readLines("../data/nfkb1_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

df <- read.table("../data/nfkb1", 
                 sep = "\t",
                 header = FALSE, 
                 stringsAsFactors = FALSE, 
                 fill = TRUE)
colnames(df) <- header_fields
head(df, 1)

# Keep a copy and filter out rows with clinvar_clnsig == "."
df <- df %>% dplyr::filter(!clinvar_clnsig == ".")
df |> count(clinvar_clnsig)

# Bar plot: Count of ClinVar Clinical Significance
p_count <- ggplot(df, aes(x = clinvar_clnsig, fill = clinvar_clnsig)) +
  geom_bar( 
    color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = "Count of ClinVar Clinical Significance NFKB1",
       subtitle = paste0("Condition: population size " , population_size)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/nfkb1_clinvar_count.png", plot = p_count, width = 6, height = 5)


# Data Preparation for Population-Level Calculations ----
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df <- df %>% select(`pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)

# keep just one transcript allele for simplicity
df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 1)
df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 1)
df$Inheritance <- "AD"  # setting inheritance to AD for demonstration
head(df)

# Filter for a specific ClinVar category: Uncertain_significance (example)
# df <- subset(df, clinvar_clnsig == "Pathogenic")
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)

# zeros? ----

# Substitute zero allele frequency with a synthesized minimum allele frequency.
# Here, we assume that if gnomAD_genomes_AF is 0, the minimum detectable AF is estimated as 1/(AN + 1),
# which is a conservative proxy based on the allele number (AN). Otherwise we may have no pathogenic candidates and we are most interested in a safe conservative estimate.
df$gnomAD_genomes_AN <- as.numeric(df$gnomAD_genomes_AN)
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)

max_af <- max(df$gnomAD_genomes_AN) |> as.numeric()

df$gnomAD_genomes_AF <- ifelse(df$gnomAD_genomes_AF == 0,
                                          1 / (max_af + 1),
                                          df$gnomAD_genomes_AF)

# Recalculate disease probability using the synthesized allele frequency
df$disease_prob <- ifelse(df$Inheritance %in% c("AD", "X-linked"),
                               df$gnomAD_genomes_AF,
                               df$gnomAD_genomes_AF^2)

# Calculate expected cases and probability of at least one affected individual in the population
df$expected_cases <- population_size * df$disease_prob
df$prob_at_least_one <- 1 - (1 - df$disease_prob)^population_size

# View the recalculated key results
print(df[, c("gnomAD_genomes_AF", "Inheritance", "disease_prob", 
                  "expected_cases", "prob_at_least_one")]) |> head()

# populations ----
# population_size <- 83702  
# population_size <- 69433632 # population_UK
# 
# df$disease_prob <- ifelse(df$Inheritance %in% c("AD", "X-linked"),
#                                df$gnomAD_genomes_AF,
#                                df$gnomAD_genomes_AF^2)
# df$expected_cases <- population_size * df$disease_prob
# df$prob_at_least_one <- 1 - (1 - df$disease_prob)^population_size

# View key results
# print(df[, c("gnomAD_genomes_AF", "Inheritance", "disease_prob", 
                  # "expected_cases", "prob_at_least_one")])

# Scatter Plots: Expected Cases and Probability vs Allele Frequency ----

# Compute threshold allele frequency at which probability is (almost) 1 (e.g., >= 0.999)
threshold_AF <- min(df$gnomAD_genomes_AF[df$prob_at_least_one >= 0.999])
threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)

p_scatter1_path <- 
  df |> filter(clinvar_clnsig == "Pathogenic") |>
  ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  # scale_x_log10() +
  ylim(0, 800) + # somthing sensible if only one value
  geom_text(aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)), 
            vjust = -1, hjust = 0.5, color = "black", size = 3) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Expected Cases",
       title = "Expected Cases vs\nAllele Frequency NFKB1",
       subtitle = paste0("Condition: population size " , population_size))

p_scatter2_path <- 
  df |> filter(clinvar_clnsig == "Pathogenic") |>
  ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, color = clinvar_clnsig)) +
  # geom_jitter( height = 0, size = 3) +
  geom_point( size = 3) +
  geom_line(aes(group = 1)) +
  # scale_x_log10() +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Probability of ≥1 Case",
       title = "Probability of At Least One\nCase vs Allele Frequency NFKB1",
       subtitle = paste0("Condition: population size " , population_size)) +
  geom_vline(xintercept = threshold_AF, linetype = "dotted", color = "black") +
  geom_text(aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
            vjust = -1, hjust = 0.5, color = "black", size = 3) +
  ylim(0,1.2)

# To display the plots:
p_scatter1_path
p_scatter2_path

# Combine and save scatter plots vertically
p_scatter1_path <- p_scatter1_path + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_scatter <- p_scatter1_path / p_scatter2_path + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
print(p_scatter)
# ggsave("../images/nfkb1_scatter_expected_prob.png", plot = p_scatter, width = 6, height = 8)

# Combine bar charts using patchwork and save

# p_bars <- (p_bar / p_prob) 

# Density histograms for Expected Cases by ClinVar Clinical Significance
p_density <- ggplot(df, aes(x = expected_cases, fill = clinvar_clnsig)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ clinvar_clnsig, scales = "free", ncol = 4) +
  scale_x_continuous(labels = function(x) format(round(x, 0), big.mark = ",")) +
  guides(fill = "none") +
  labs(x = "Expected Cases", 
       y = "Density",
       title = "Density of Expected Cases by ClinVar Clinical Significance",
       subtitle = paste0("Condition: population size ", population_size)) 

print(p_density)
# ggsave("../images/nfkb1_density_expected_cases.png", plot = p_density, width = 12, height = 5)

p_scatter_dense <- (p_density / (p_scatter1_path + p_scatter2_path)) +
  plot_layout(widths = c(1, 1), guides = 'collect', axis = "collect") +
  plot_annotation(tag_levels = 'A')
print(p_scatter_dense)
ggsave("../images/nfkb1_scatterdense_expected_prob.png", plot = p_scatter_dense, width = 12, height = 6)

 # Tally by ClinVar Category ----
df_calc <- df %>%
  
# df_calc <- df %>%
  mutate(allele_freq = as.numeric(gnomAD_genomes_AF)) %>%
  filter(!is.na(allele_freq)) %>%
  mutate(disease_prob = ifelse(Inheritance %in% c("AD", "X-linked"), 
                               allele_freq,
                               allele_freq^2),
         expected_cases = population_size * disease_prob,
         prob_at_least_one = 1 - (1 - disease_prob)^population_size)

clinvar_levels <- unique(df$clinvar_clnsig)

df_tally <- df_calc %>%
  group_by(clinvar_clnsig) %>%
  summarise(total_expected_cases = sum(expected_cases),
            overall_prob = 1 - prod(1 - disease_prob),
            .groups = "drop") %>%
  tidyr::complete(clinvar_clnsig = clinvar_levels,
                  fill = list(total_expected_cases = 0, overall_prob = 0))
print(df_tally)

# Bar plot: Total Expected Cases by ClinVar Category
p_bar <- ggplot(df_tally, aes(x = clinvar_clnsig, y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Total Expected Cases by ClinVar\nClinical Significance in NFKB1",
       subtitle = paste0("Condition: population size " , population_size)) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
p_bar
# ggsave("../images/nfkb1_bar_expected_cases.png", plot = p_bar, width = 8, height = 6)

# Bar plot: Overall Probability by ClinVar Category
p_prob <- ggplot(df_tally, aes(x = clinvar_clnsig, y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.001)), vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall\nProbability",
       title = "Overall Probability of an Affected\nBirth by ClinVar Category in NFKB1",
       subtitle = paste0("Condition: population size " , population_size)) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::percent)
p_prob
# ggsave("../images/nfkb1_bar_overall_probability.png", plot = p_prob, width = 8, height = 5)

# Combine bar charts using patchwork and save
p_bar <- p_bar + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_bars <- (p_bar / p_prob) + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
print(p_bars)
ggsave("../images/nfkb1_combined_bar_charts.png", plot = p_bars, width = 6, height = 7)
