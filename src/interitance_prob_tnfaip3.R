library(ggplot2); theme_set(theme_bw())
# Ensure the output directory exists
if(!dir.exists("../images/")) dir.create("../images/")

# UK population
population_size <- 69433632
# population_size <- 83702  # births in 2023
# 
# -------------------------
# DBNSFP Data Import and Filtering
# -------------------------
library(ggplot2); theme_set(theme_bw())
library(dplyr)

header_line <- readLines("~/Desktop/dbnsfp/data/tnfaip3_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

df <- read.table("~/Desktop/dbnsfp/data/tnfaip3", 
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
     title = "Count of ClinVar Clinical Significance TNFAIP3"
     # subtitle = paste0("Condition: population size " , population_size)
     ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")
p_count
ggsave("../images/tnfaip3_clinvar_count.png", plot = p_count, width = 6, height = 5)

# -------------------------
# Data Preparation for Population-Level Calculations
# -------------------------
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df <- df %>% select(`pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig)
df$Inheritance <- "AD"  # setting inheritance to AD for demonstration
head(df)

# Filter for a specific ClinVar category: Uncertain_significance (example)
df_path <- subset(df, clinvar_clnsig == "Uncertain_significance")
# df_path <- subset(df, clinvar_clnsig == "Pathogenic")
df_path$gnomAD_genomes_AF <- as.numeric(df_path$gnomAD_genomes_AF)

df_path$disease_prob <- ifelse(df_path$Inheritance %in% c("AD", "X-linked"),
                               df_path$gnomAD_genomes_AF,
                               df_path$gnomAD_genomes_AF^2)
df_path$expected_cases <- population_size * df_path$disease_prob
df_path$prob_at_least_one <- 1 - (1 - df_path$disease_prob)^population_size

# View key results
print(df_path[, c("gnomAD_genomes_AF", "Inheritance", "disease_prob", 
                  "expected_cases", "prob_at_least_one")]) |> head()

# -------------------------
# Scatter Plots: Expected Cases and Probability vs Allele Frequency
# -------------------------
library(patchwork)
df_path_nonNA <- subset(df_path, !is.na(gnomAD_genomes_AF))

# Compute threshold allele frequency at which probability is (almost) 1 (e.g., >= 0.999)
threshold_AF <- min(df_path_nonNA$gnomAD_genomes_AF[df_path_nonNA$prob_at_least_one >= 0.999])
threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)

p_scatter1 <- ggplot(df_path_nonNA, aes(x = gnomAD_genomes_AF, y = expected_cases)) +
  geom_point(size = 3, color = "#FF3366") +
  geom_line(aes(group = 1, color = "#FF3366")) +
  scale_x_log10() +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Expected Cases",
       title = "Expected Cases vs Allele Frequency \nTNFAIP3: Uncertain significance",
       subtitle = paste0("Condition: population size " , population_size))  +
  guides(color = "none")

p_scatter2 <- ggplot(df_path_nonNA, aes(x = gnomAD_genomes_AF, y = prob_at_least_one)) +
  geom_point(size = 3, color = "#FF3366") +
  geom_line(aes(group = 1, color = "#FF3366")) +
  scale_x_log10() +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Probability of ≥1 Case",
       title = "Probability of At Least One Case vs Allele Frequency\nTNFAIP3: Uncertain significance",
       subtitle = paste0("Condition: population size " , population_size)) +
  geom_vline(xintercept = threshold_AF, linetype = "dotted", color = "black") +
  geom_text(aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
            vjust = -1, hjust = 0.5, color = "black", size = 3) +
scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  guides(color = "none")

# To display the plots:
p_scatter1
p_scatter2

# Combine and save scatter plots vertically
p_scatter <- p_scatter1 / p_scatter2 + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
print(p_scatter)
ggsave("../images/tnfaip3_scatter_expected_prob.png", plot = p_scatter, width = 6, height = 8)

# -------------------------
# Tally by ClinVar Category
# -------------------------
df_calc <- df %>%
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
  # geom_text(aes(label = round(total_expected_cases, 1)), vjust = -0.5, size = 3.5) +
  geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Total Expected Cases by\nClinVar Clinical Significance TNFAIP3",
       subtitle = paste0("Condition: population size " , population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
p_bar
# ggsave("../images/bar_expected_cases.png", plot = p_bar, width = 8, height = 6)

# Bar plot: Overall Probability by ClinVar Category
p_prob <- ggplot(df_tally, aes(x = clinvar_clnsig, y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.001)), vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall Probability",
       title = "Overall Probability of an\nAffected Birth by ClinVar Category TNFAIP3",
       subtitle = paste0("Condition: population size " , population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = scales::percent)
p_prob
# ggsave("../images/bar_overall_probability.png", plot = p_prob, width = 8, height = 5)

# Combine bar charts using patchwork and save
p_bars <- p_bar + p_prob  + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
print(p_bars)
ggsave("../images/tnfaip3_combined_bar_charts.png", plot = p_bars, width = 10, height = 6)
