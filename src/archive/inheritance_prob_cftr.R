
library(ggplot2); theme_set(theme_bw())
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)

if(!dir.exists("../images/")) dir.create("../images/")

population_size <- 69433632

# Data Import and Filtering ----
header_line <- readLines("../data/cftr_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

df <- read.table("../data/cftr", 
                 sep = "\t",
                 header = FALSE, 
                 stringsAsFactors = FALSE, 
                 fill = TRUE)
colnames(df) <- header_fields

# Use NFE allele frequency column as our AF and remove the original AF column
df <- df |> select(-gnomAD_genomes_AF)
colnames(df)[colnames(df) == 'gnomAD_genomes_NFE_AF'] <- 'gnomAD_genomes_AF'

df <- df %>% filter(!clinvar_clnsig == ".")
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df$gnomAD_genomes_AN <- as.numeric(df$gnomAD_genomes_AN)
df <- df %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, 
                    clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)
df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 2)
df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 2)
df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)
df$Inheritance <- "AR"  # autosomal recessive

# Substitute zero AF with a minimal detectable value.
max_af <- max(df$gnomAD_genomes_AN, na.rm = TRUE)
df$gnomAD_genomes_AF <- ifelse(df$gnomAD_genomes_AF == 0,
                               1 / (max_af + 1),
                               df$gnomAD_genomes_AF)

# --- Compute Total Allele Frequency per Gene (for all variants) ---
# Here we do NOT filter by clinvar_clnsig; we use all variants with nonzero AN.
gene_totals <- df %>% 
  filter(genename == "CFTR", gnomAD_genomes_AN > 0) %>%
  group_by(genename) %>%
  summarise(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE), .groups = "drop")

# Merge gene totals into the main dataframe.
df <- df %>% left_join(gene_totals, by = "genename")

# --- Calculate Occurrence Probability ---
# For AR conditions, an event (i.e. having two pathogenic alleles) can occur via:
#   - Homozygous: (p)^2
#   - Compound heterozygous: 2 * p * (total_AF - p)
df <- df %>% mutate(
  other_AF = ifelse(Inheritance == "AR", total_AF - gnomAD_genomes_AF, 0),
  occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
                           gnomAD_genomes_AF,
                           pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * other_AF, 0))
)

# Calculate expected cases and probability of at least one event.
df <- df %>% mutate(
  expected_cases = population_size * occurrence_prob,
  prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
)

# View key results.
df |> head()

# --- Tally by ClinVar Category ---
clinvar_levels <- unique(df$clinvar_clnsig)

# Keep in mind that if there is overlap (non-independence) between variants (for example, if the same individual might carry two different variants), then this approach may overestimate the probability. But as an approximation under the assumption of independence, the code is correct for recessive disease.
df_tally <- df %>%
  group_by(clinvar_clnsig) %>%
  summarise(total_expected_cases = sum(expected_cases),
            overall_prob = 1 - prod(1 - occurrence_prob),
            .groups = "drop") %>%
  complete(clinvar_clnsig = clinvar_levels,
           fill = list(total_expected_cases = 0, overall_prob = 0))
print(df_tally)

# --- Bar Plots ---
# Total Expected Cases by ClinVar Category.
p_bar <- ggplot(df_tally, aes(
  x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Total Expected Cases by ClinVar\nClinical Significance in CFTR",
       subtitle = paste0("Population Size: ", population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
p_bar

# Overall Probability by ClinVar Category.
p_prob <- ggplot(df_tally, aes(
  x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.01)), vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall\nProbability",
       title = "Overall Probability of an Affected\nBirth by ClinVar Category in CFTR",
       subtitle = paste0("Population Size: ", population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::percent)
p_prob

# Combine bar charts.
p_bar_combined <- (p_bar / p_prob) + 
  plot_layout(guides = 'collect', axis = "collect") + 
  plot_annotation(tag_levels = 'A', title = "Recessive disease gene")
print(p_bar_combined)
ggsave("../images/cftr_combined_bar_charts.png", plot = p_bar_combined, width = 7, height = 8)






het_rate <- 0.0083       # 0.83% heterozygous frequency (as a proportion)
p <- het_rate / 2        # allele frequency p
hom_rate <- p^2          # homozygous frequency for one allele
comphet_rate <- 2 * p^2  # assuming two similar alleles, frequency of compound heterozygotes

cat("Allele frequency p =", p, "\n")
cat("Homozygous frequency =", hom_rate * 100, "%\n")
cat("Compound heterozygous frequency =", comphet_rate * 100, "%\n")
cat("Total recessive frequency =", (hom_rate + comphet_rate) * 100, "%\n")











# Validation ----

# Expected values of pathogenic variants
# https://www.cysticfibrosis.org.uk/sites/default/files/2024-10/CFT_2023_Annual_Data_Report_v9.pdf

#Number of cases in 2023: 11318
# Sec 1.47 CFTR variant combinations in the UK population
# This tabulation shows the proportion(%) of patients with the most common CFTR variant combinations in their genotype. For example, 4.0% of the UK population have one copy of F508del and one copy of G551D.

# print("Expected:\n
# variant	proportion
# R117H_het	5.2
# G551D_het	4.4
# R117H_hom	0.1
# G551D_hom	0.2
# R117H_G551D	0.2")
# 
# number_cases <- 11318
# proportion_cases_hom_R117H <- 0.1
# proportion_cases_het_R117H <-	5.2
# number_hom_R117H <- number_cases * proportion_cases_hom_R117H
# number_het_R117H <- number_cases * proportion_cases_het_R117H
# 
# print(population_size)
# print(number_cases)
# print(proportion_cases_hom_R117H)
# print(number_hom_R117H)
# print(number_het_R117H)

# Appendix 3: Full list of CFTR variants in the UK CF population
# The table below shows the number of people with CF who carry at least one of each variant.The groups are not mutually exclusive, as people with heterozygous variants appear twice in the table. Most common SNV: c.350G->A, p.Arg117His, 714, 6.3%

# 
# # Ensure the output directory exists
# if(!dir.exists("../images/")) dir.create("../images/")
# 
# # UK population
# population_size <- 69433632
# # population_size <- 83702  # births in 2023
# 
# # DBNSFP Data Import and Filtering ----
# header_line <- readLines("../data/cftr_head", n = 1)
# header_line <- sub("^#", "", header_line)
# header_fields <- strsplit(header_line, "\t")[[1]]
# rm(header_line)
# 
# df <- read.table("../data/cftr", 
#                  sep = "\t",
#                  header = FALSE, 
#                  stringsAsFactors = FALSE, 
#                  fill = TRUE)
# colnames(df) <- header_fields
# head(df, 1)
# 
# # Project specific population ----
# df <- df |> select(-gnomAD_genomes_AF)
# colnames(df)[colnames(df) == 'gnomAD_genomes_NFE_AF'] <- 'gnomAD_genomes_AF'
# 
# # Keep a copy and filter out rows with clinvar_clnsig == "."
# df <- df %>% filter(!clinvar_clnsig == ".")
# df |> count(clinvar_clnsig)
# 
# # Bar plot: Count of ClinVar Clinical Significance
# p_count <- ggplot(df, aes(
#   x = str_wrap(paste0(gsub("_", " ", clinvar_clnsig)), width = 20), 
#   fill = clinvar_clnsig)) +
#   geom_bar(color = "black") +
#   geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
#   labs(x = "ClinVar Clinical Significance",
#        y = "Count",
#        title = "Count of ClinVar Clinical Significance CFTR") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(fill = "none")
# 
# p_count
# ggsave("../images/cftr_clinvar_count.png", plot = p_count, width = 5, height = 4)
# 
# # Data Preparation for Population-Level Calculations ----
# df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
# df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
# df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
# df$gnomAD_genomes_AN <- as.numeric(df$gnomAD_genomes_AN)
# df <- df %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, 
#                     clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)
# 
# # Keep just one transcript allele for simplicity
# df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 2)  # cftr registry transcript
# df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 2)  # cftr registry transcript
# df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)
# df$Inheritance <- "AR"  # setting inheritance to AR for autosomal recessive conditions
# 
# # Substitute zero allele frequency with a synthesized minimum allele frequency.
# max_af <- max(df$gnomAD_genomes_AN, na.rm = TRUE)
# df$gnomAD_genomes_AF <- ifelse(df$gnomAD_genomes_AF == 0,
#                                1 / (max_af + 1),
#                                df$gnomAD_genomes_AF)
# 
# # Update occurrence probability calculation ----
# # First, get the total pathogenic allele frequency per gene.
# # We filter for pathogenic variants (clinvar_clnsig == "Pathogenic") with nonzero allele numbers.
# df_patho <- df %>% filter(gnomAD_genomes_AN > 0, clinvar_clnsig == "Pathogenic")
# gene_patho <- df_patho %>%
#   group_by(genename) %>%
#   summarise(total_pathogenic_AF = sum(gnomAD_genomes_AF, na.rm = TRUE), .groups = "drop")
# 
# # Join total pathogenic allele frequency into the main dataframe.
# df <- df %>% left_join(gene_patho, by = "genename")
# 
# # For AR conditions, an event (carriage) occurs when any two alleles are present.
# # This can be due to homozygous variants: (AF)^2, or compound heterozygous variants:
# # 2 * AF * (total_pathogenic_AF - AF).
# # For AD and X-linked, we simply use AF.
# # We call this the "occurrence probability" (always ≥ 0).
# df <- df %>% mutate(
#   other_AF = ifelse(Inheritance == "AR", total_pathogenic_AF - gnomAD_genomes_AF, 0),
#   occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
#                            gnomAD_genomes_AF,
#                            pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * other_AF, 0))
# )
# 
# # Calculate expected cases and probability of at least one event in the population.
# df <- df %>% mutate(
#   expected_cases = population_size * occurrence_prob,
#   prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
# )
# 
# # View the recalculated key results.
# head(df)
# 
# # Scatter Plots: Expected Cases and Probability vs Allele Frequency ----
# 
# # Compute threshold allele frequency at which probability is (almost) 1 (e.g., ≥ 0.999)
# threshold_AF <- min(df$gnomAD_genomes_AF[df$prob_at_least_one >= 0.999])
# threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)
# 
# unique_labels <- df %>% 
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   distinct(gnomAD_genomes_AF, expected_cases)
# 
# p_scatter1_path <- 
#   df %>% filter(clinvar_clnsig == "Pathogenic") %>%
#   ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
#   geom_point(size = 3) +
#   geom_line(aes(group = 1)) +
#   # scale_x_log10() +
#   ylim(0, 2400) +  # something sensible if only one value
#   ggrepel::geom_text_repel(data = unique_labels, 
#             aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)), 
#             vjust = -1, hjust = 0.5, colour = "black", size = 3) +
#   labs(x = "Allele Frequency (log scale)",
#        y = "Expected Cases",
#        title = "Expected Cases vs\nAllele Frequency CFTR",
#        subtitle = paste0("Condition: population size ", population_size))
# 
# vline_label <- tibble(threshold_AF = threshold_AF, 
#                       threshold_AF_label = threshold_AF_label)
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
#   geom_vline(xintercept = threshold_AF, linetype = "dotted", colour = "black") +
#   geom_text(data = vline_label, 
#             aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
#             vjust = -1, hjust = 0.5, colour = "black", size = 3) +
#   ylim(0, 1.2)
# 
# # Display the scatter plots.
# p_scatter1_path
# p_scatter2_path
# 
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
# # Tally by ClinVar Category ----
# df_calc <- df %>%
#   mutate(allele_freq = as.numeric(gnomAD_genomes_AF)) %>%
#   filter(!is.na(allele_freq)) %>%
#   mutate(occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"), 
#                                   allele_freq,
#                                   pmax(allele_freq^2 + 2 * allele_freq * (total_pathogenic_AF - allele_freq), 0)),
#          expected_cases = population_size * occurrence_prob,
#          prob_at_least_one = 1 - (1 - occurrence_prob)^population_size)
# 
# clinvar_levels <- unique(df$clinvar_clnsig)
# 
# df_tally <- df_calc %>%
#   group_by(clinvar_clnsig) %>%
#   summarise(total_expected_cases = sum(expected_cases),
#             overall_prob = 1 - prod(1 - occurrence_prob),
#             .groups = "drop") %>%
#   complete(clinvar_clnsig = clinvar_levels,
#            fill = list(total_expected_cases = 0, overall_prob = 0))
# print(df_tally)
# 
# # Bar plot: Total Expected Cases by ClinVar Category.
# p_bar <- ggplot(df_tally, aes(
#   x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
#   y = total_expected_cases, fill = clinvar_clnsig)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
#             vjust = -0.5, size = 3.5) +
#   labs(x = "ClinVar Clinical Significance",
#        y = "Total Expected\nCases",
#        title = "Total Expected Cases by ClinVar\nClinical Significance in CFTR",
#        subtitle = paste0("Condition: population size ", population_size)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(fill = "none") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
# 
# p_bar
# 
# # Bar plot: Overall Probability by ClinVar Category.
# p_prob <- ggplot(df_tally, aes(
#   x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
#   y = overall_prob, fill = clinvar_clnsig)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.01)), vjust = -0.5, size = 3.5) +
#   labs(x = "ClinVar Clinical Significance",
#        y = "Overall\nProbability",
#        title = "Overall Probability of an Affected\nBirth by ClinVar Category in CFTR",
#        subtitle = paste0("Condition: population size ", population_size)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   guides(fill = "none") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::percent)
# 
# p_prob
# 
# # Combine bar charts using patchwork and save.
# p_bar <- p_bar + 
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())
# print(p_bars)
# ggsave("../images/cftr_combined_bar_charts.png", plot = p_bars, width = 6, height = 6)



# other ----
# For CFTR, restrict to known pathogenic variants using gnomAD_genomes_AF.
cftr_patho <- df %>% 
  filter(genename == "CFTR", clinvar_clnsig == "Pathogenic", gnomAD_genomes_AN > 0)

# Total pathogenic allele frequency for CFTR.
total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)

# Focus on p.Arg117His.
df_arg117his <- cftr_patho %>% 
  filter(HGVSp_VEP == "p.Arg117His")

# p = allele frequency for p.Arg117His; 
# q = sum of all other pathogenic allele frequencies in CFTR.
p <- df_arg117his$gnomAD_genomes_AF
q <- total_AF - p

# Expected genotype counts in the general population using Hardy–Weinberg:
expected_hom <- population_size * p^2             # homozygous for p.Arg117His
expected_comphet <- population_size * 2 * p * q      # compound heterozygotes carrying p.Arg117His with any other pathogenic allele
expected_total_genotypes <- expected_hom + expected_comphet

# Option 1: Simple mortality adjustment (subtracting 0.4%):
mortality_rate <- 0.004
expected_after_mortality_simple <- expected_total_genotypes * (1 - mortality_rate)

# Option 2: Exponential survival model:
median_age <- 22
# Survival factor = exp(-annual mortality rate * median_age)
survival_factor <- exp(-mortality_rate * median_age)
adjusted_expected_after_mortality <- expected_total_genotypes * survival_factor

# p.Arg117His is known to have reduced penetrance.
# Observed in the registry: 714 cases among 11318 total CF cases (~6.3%).
observed_R117 <- 714
source_text <- "Source: UK Cystic Fibrosis Registry 2023 Annual Data Report, October 2024\nObserved: p.Arg117His in cases 714/11318 (~6.3%)"

# Create a summary data frame that includes:
# - Expected Homozygous
# - Expected Compound Het
# - Expected Total Genotypes (unadjusted)
# - Expected Total Genotypes After Mortality (Exponential Survival Adjustment)
# - Observed Count

summary_df <- data.frame(
  Category = c("Expected\nHomozygous", 
               "Expected\nCompound Het", 
               "Expected\nTotal Genotypes", 
               "After Mortality\n(Exponential)", 
               "Observed*"),
  Count = c(expected_hom, 
            expected_comphet, 
            expected_total_genotypes, 
            adjusted_expected_after_mortality, 
            observed_R117)
)

# Preserve the specified order by setting Category as a factor with levels.
summary_df$Category <- factor(summary_df$Category, levels = c("Expected\nHomozygous", 
                                                              "Expected\nCompound Het", 
                                                              "Expected\nTotal Genotypes", 
                                                              "After Mortality\n(Exponential)", 
                                                              "Observed*"))

# Plot the summary as a bar chart with annotations and source.
p_var <- ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(round(Count, 0), format = "d", big.mark = ",")),
            vjust = -0.5, size = 4) +
  labs(title = "Expected Genotype Counts for p.Arg117His in CFTR",
       subtitle = paste0("Population Size: ", population_size, 
                         "\nSurvival Factor (to median age 22): ", 
                         formatC(survival_factor, format = "f", digits = 3)),
       caption = source_text,
       x = "Genotype Category",
       y = "Count") +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) 

p_var
ggsave("../images/cftr_validation_pArg117His.png", plot = p_var, width = 6, height = 6)


