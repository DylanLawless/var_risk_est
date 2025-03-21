library(ggplot2); theme_set(theme_bw())
library(patchwork)

# Ensure the output directory exists
if(!dir.exists("../images/")) dir.create("../images/")

# UK population
population_size <- 69433632
# population_size <- 83702  # births in 2023

# DBNSFP Data Import and Filtering ----
library(dplyr)

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
head(df, 1)

# project specific population ----
df <- df |> select(-gnomAD_genomes_AF)
colnames(df)[colnames(df) == 'gnomAD_genomes_NFE_AF'] <- 'gnomAD_genomes_AF'


# Keep a copy and filter out rows with clinvar_clnsig == "."
df <- df %>% dplyr::filter(!clinvar_clnsig == ".")
df |> count(clinvar_clnsig)

# Bar plot: Count of ClinVar Clinical Significance
p_count <- ggplot(df, aes(
  x = stringr::str_wrap(paste0(gsub("_", " ", clinvar_clnsig)), width = 20), 
  fill = clinvar_clnsig)) +
  geom_bar( 
    color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = "Count of ClinVar Clinical Significance CFTR",
       # subtitle = paste0("Condition: population size " , population_size)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/cftr_clinvar_count.png", plot = p_count, width = 5, height = 4)


# Data Preparation for Population-Level Calculations ----
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df <- df %>% select(genename,`pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)

# keep just one transcript allele for simplicity
df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 2) # cftr registry transcript
df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 2) # cftr registry transcript
df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)
df$Inheritance <- "AR"  # setting inheritance to AD for demonstration
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
# gnomad AF -----
df$gnomAD_genomes_AF <- ifelse(df$gnomAD_genomes_AF == 0,
                               1 / (max_af + 1),
                               df$gnomAD_genomes_AF)
# 
# ggplot(df, aes(x = gnomAD_genomes_AF)) +
#   geom_histogram(bins = 30)

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
# Scatter Plots: Expected Cases and Probability vs Allele Frequency ----

# Compute threshold allele frequency at which probability is (almost) 1 (e.g., >= 0.999)
threshold_AF <- min(df$gnomAD_genomes_AF[df$prob_at_least_one >= 0.999])
threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)

unique_labels <- df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>%
  distinct(gnomAD_genomes_AF, expected_cases)

p_scatter1_path <- 
  df |> filter(clinvar_clnsig == "Pathogenic") |>
  ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  # scale_x_log10() +
  ylim(0, 800) + # somthing sensible if only one value
  geom_text(data = unique_labels, 
            aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Expected Cases",
       title = "Expected Cases vs\nAllele Frequency CFTR",
       subtitle = paste0("Condition: population size " , population_size))



vline_label <- tibble(threshold_AF = threshold_AF, 
                      threshold_AF_label = threshold_AF_label)

p_scatter2_path <- 
  df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>% 
  ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, colour = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Probability of ≥1 Case",
       title = "Probability of At Least One\nCase vs Allele Frequency CFTR",
       subtitle = paste0("Condition: population size ", population_size)) +
  geom_vline(xintercept = threshold_AF, linetype = "dotted", colour = "black") +
  geom_text(data = vline_label, 
            aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  ylim(0, 1.2)

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

p_scatter_dense <- (p_density / (p_scatter1_path + p_scatter2_path)) +
  plot_layout(widths = c(1, 1), guides = 'collect', axis = "collect") +
  plot_annotation(tag_levels = 'A')
print(p_scatter_dense)
ggsave("../images/cftr_scatterdense_expected_prob.png", plot = p_scatter_dense, width = 12, height = 6)

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
p_bar <- ggplot(df_tally, aes(
  x = stringr::str_wrap(paste0(gsub("_", " ", clinvar_clnsig)), width = 20), , 
  y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Total Expected Cases by ClinVar\nClinical Significance in CFTR",
       subtitle = paste0("Condition: population size " , population_size)) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
p_bar

# Bar plot: Overall Probability by ClinVar Category
p_prob <- ggplot(df_tally, aes(
  x = stringr::str_wrap(paste0(gsub("_", " ", clinvar_clnsig)), width = 20), , 
  y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.001)), vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall\nProbability",
       title = "Overall Probability of an Affected\nBirth by ClinVar Category in CFTR",
       subtitle = paste0("Condition: population size " , population_size)) +
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::percent)
p_prob
# ggsave("../images/CFTR_bar_overall_probability.png", plot = p_prob, width = 8, height = 5)

# Combine bar charts using patchwork and save
p_bar <- p_bar + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_bars <- (p_bar / p_prob) + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
print(p_bars)
ggsave("../images/CFTR_combined_bar_charts.png", plot = p_bars, width = 6, height = 6)










# test 1 ----

# 
# # Expected values of pathogenic variants
# # https://www.cysticfibrosis.org.uk/sites/default/files/2024-10/CFT_2023_Annual_Data_Report_v9.pdf
# 
# #Number of cases in 2023: 11318
# # Sec 1.47 CFTR variant combinations in the UK population
# # This tabulation shows the proportion(%) of patients with the most common CFTR variant combinations in their genotype. For example, 4.0% of the UK population have one copy of F508del and one copy of G551D.
# 
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
# print(proportion_cases_homozygous_R117H)
# print(number_hom_R117H)
# print(number_het_R117H)
# 
# df |> filter(HGVSp_VEP == "p.Arg117His")
# 
# # Calculate expected homozygous and heterozygous counts.
# # For autosomal recessive inheritance:
# # - homozygous: population_size * (allele frequency)^2
# # - heterozygous: population_size * 2 * allele frequency * (1 - allele frequency)
# df_test <- df %>%
#   mutate(expected_hom = population_size * (gnomAD_genomes_AF)^2,
#          expected_het = population_size * 2 * gnomAD_genomes_AF * (1 - gnomAD_genomes_AF))
# 
# # Filter for the p.Arg117His variant
# df_arg117his <- df_test %>% filter(HGVSp_VEP == "p.Arg117His")
# 
# print(df_arg117his[, c("HGVSp_VEP", "gnomAD_genomes_AF", "expected_hom", "expected_het")])
# 

# # test 2 ----
# Appendix 3: Full list of CFTR variants in the UK CF population
# The table below shows the number of people with CF who carry at least one of each variant.
# The groups are not mutually exclusive, as people with heterozygous variants appear twice in the table.
# 
# c.350G->A, p.Arg117His, 714, 6.3%

# # Calculate expected counts for heterozygous and homozygous carriage in the general population.
# # For autosomal recessive variants:
# #   - Homozygous count = population_size * (AF)^2
# #   - Heterozygous count = population_size * 2 * AF * (1 - AF)
# df_test <- df %>%
#   mutate(expected_hom = population_size * (gnomAD_genomes_AF)^2,
#          expected_het = population_size * 2 * gnomAD_genomes_AF * (1 - gnomAD_genomes_AF),
#          expected_total = expected_hom + expected_het)
# 
# # Focus on the variant c.350G>A, p.Arg117His.
# df_arg117his <- df_test %>% filter(HGVSp_VEP == "p.Arg117His")
# 
# # Print the general population estimates for this variant.
# print(df_arg117his[, c("HGVSp_VEP", "gnomAD_genomes_AF", "expected_hom", "expected_het", "expected_total")])
# 
# 
# 
# 
# # Calculate expected counts for heterozygous and homozygous carriage using HWE.
# df_test <- df %>%
#   mutate(expected_hom = population_size * (gnomAD_genomes_AF)^2,
#          expected_het = population_size * 2 * gnomAD_genomes_AF * (1 - gnomAD_genomes_AF),
#          expected_total = expected_hom + expected_het)
# 
# # Focus on the variant c.350G>A, p.Arg117His.
# df_arg117his <- df_test %>% filter(HGVSp_VEP == "p.Arg117His")
# 
# # Print the general population estimates for this variant.
# print(df_arg117his[, c("HGVSp_VEP", "gnomAD_genomes_AF", "expected_hom", "expected_het", "expected_total")])
# 
# # To obtain an equivalent measure in the CF population, we can rescale the expected_total.
# # Assuming that the allele-derived expected count in the general population is proportional to its contribution among all CFTR variants,
# # we compute the overall expected CF cases in the general population from our data.
# total_expected_CF <- sum(df$expected_cases, na.rm = TRUE)
# 
# # For the specific variant, the proportion among all expected CF cases is:
# variant_fraction <- df_arg117his$expected_cases / total_expected_CF
# 
# # Multiply this fraction by the registered CF cases to estimate the number of CF patients carrying at least one copy.
# expected_variant_in_CF <- number_cases * variant_fraction
# 
# print(expected_variant_in_CF)
# 
# 
# 
# 
# 
# 
# # Known values from the registry:
# number_cases <- 11318
# registry_count_R117H <- 714  # CF patients with at least one p.Arg117His allele
# 
# # For p.Arg117His (from DF):
# # expected_hom: expected homozygous carriers (all assumed to be on registry)
# # expected_het: expected heterozygous carriers in the general population
# # These were previously computed:
# #   expected_hom = population_size * (gnomAD_genomes_AF)^2
# #   expected_het = population_size * 2 * gnomAD_genomes_AF * (1 - gnomAD_genomes_AF)
# #
# # For the heterozygotes, only a fraction will be compound heterozygotes (i.e. have a second pathogenic variant).
# # This fraction (compound_factor) can be estimated from the registry data:
# compound_factor <- (registry_count_R117H - df_arg117his$expected_hom) / df_arg117his$expected_het
# 
# # Now, compute the estimated registry count for p.Arg117His from our DF:
# expected_registry_count <- df_arg117his$expected_hom + df_arg117his$expected_het * compound_factor
# 
# print(expected_registry_count)
# 
# 
# 
# # test 3 ----
# 
# # Define pathogenic statuses (adjust as needed)
# pathogenic_status <- c("Pathogenic")
# 
# # Filter for pathogenic variants
# df_patho <- df %>% 
#   filter(gnomAD_genomes_AN > 0) %>%
#   filter(clinvar_clnsig %in% pathogenic_status)
# 
# # Sum allele frequencies per gene for all pathogenic variants
# gene_patho <- df_patho %>%
#   group_by(genename) %>%
#   summarise(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE), .groups = "drop")
# 
# # Focus on p.Arg117His
# df_arg117his <- df_patho %>% 
#   filter(HGVSp_VEP == "p.Arg117His") %>%
#   left_join(gene_patho, by = "genename") %>%
#   # Compute allele frequency for other pathogenic variants in the gene
#   mutate(other_AF = total_AF - gnomAD_genomes_AF,
#          # Expected recessive disease: homozygous (p^2) + compound heterozygous (2*p*q)
#          expected_recessive = population_size * (gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * other_AF))
# 
# print(df_arg117his[, c("genename", "HGVSp_VEP", "gnomAD_genomes_AF", "expected_recessive")])
# 
# # For CF, we expect disease only when two pathogenic alleles are present.
# # The overall expected CF prevalence in the general population is given by:
# #   CF_total_expected = population_size * (total_pathogenic_AF)^2
# #
# # For p.Arg117His carriers, the expected count (in the general population) is:
# #   expected_R117 = population_size * (p_R117^2 + 2 * p_R117 * (total_pathogenic_AF - p_R117))
# #
# # But since the CF registry reports numbers among confirmed CF patients (number_cases),
# # we need to compute the fraction of p.Arg117His cases among all expected CF cases:
# #
# #   fraction_R117 = expected_R117 / CF_total_expected
# #
# # Then, the estimated number in the CF registry carrying p.Arg117His is:
# #
# #   expected_R117_in_CF = number_cases * fraction_R117
# 
# # First, compute total pathogenic AF for CFTR (from our filtered data):
# total_pathogenic_AF <- df_patho %>%
#   filter(genename == "CFTR") %>%
#   summarise(total_AF = sum(gnomAD_genomes_AF, na.rm = TRUE)) %>%
#   pull(total_AF)
# 
# # Extract p.Arg117His allele frequency:
# p_R117 <- df_arg117his$gnomAD_genomes_AF
# 
# # Compute expected counts in the general population:
# expected_R117 <- population_size * (p_R117^2 + 2 * p_R117 * (total_pathogenic_AF - p_R117))
# CF_total_expected <- population_size * (total_pathogenic_AF)^2
# 
# # Fraction of p.Arg117His cases among all CF cases:
# fraction_R117 <- expected_R117 / CF_total_expected
# # fraction_R117 <- expected_R117 / number_cases
# 
# # Estimated count in the CF registry:
# number_cases <- 11318
# expected_R117_in_CF <- number_cases * fraction_R117
# 
# print(expected_R117_in_CF)
# 
# 
# # test 5 ----
# # For CFTR, restrict to known pathogenic variants with nonzero allele counts.
# cftr_patho <- df %>% 
#   filter(genename == "CFTR", clinvar_clnsig %in% "Pathogenic", gnomAD_genomes_AN > 0)
# 
# # Total pathogenic allele frequency (for CFTR)
# total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)
# 
# # Focus on p.Arg117His
# df_arg117his <- cftr_patho %>% 
#   filter(HGVSp_VEP == "p.Arg117His")
# 
# # p = frequency of p.Arg117His and q = frequency of all other pathogenic alleles
# p <- df_arg117his$gnomAD_genomes_AF
# q <- total_AF - p
# 
# # Expected counts in the general population using Hardy–Weinberg:
# expected_hom <- population_size * p^2
# expected_comphet <- population_size * 2 * p * q
# expected_R117 <- expected_hom + expected_comphet
# 
# # Our registry (CF patient) data:
# number_cases <- 11318
# # Known CF registry count for p.Arg117His carriers is 714, which is ~6.3% of CF patients.
# observed_R117 <- 714
# 
# # The ratio of observed to expected (from allele frequencies) can be used as a scaling factor.
# scale_factor <- observed_R117 / expected_R117
# 
# # Adjusted expected number of CF patients carrying p.Arg117His (homozygous or compound heterozygous)
# adjusted_expected_R117 <- expected_R117 * scale_factor
# 
# print(adjusted_expected_R117)
# 
# 
# 
# 
# 
# 
# # For CFTR, restrict to known pathogenic variants using gnomAD_genomes_AF.
# cftr_patho <- df %>% 
#   filter(genename == "CFTR", clinvar_clnsig %in% "Pathogenic", gnomAD_genomes_AN > 0)
# 
# # This line sums the allele frequencies of all pathogenic CFTR variants
# # (as filtered in cftr_patho) to get the total pathogenic allele frequency.
# total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)
# 
# # Focus on p.Arg117His.
# df_arg117his <- cftr_patho %>% 
#   filter(HGVSp_VEP == "p.Arg117His")
# 
# # p = allele frequency for p.Arg117His; q = allele frequency for all other pathogenic variants in CFTR.
# p <- df_arg117his$gnomAD_genomes_AF
# q <- total_AF - p
# 
# # Expected genotype frequency (per Hardy–Weinberg) for CF causing genotypes involving p.Arg117His:
# #   - Homozygous: p^2
# #   - Compound heterozygous (with another pathogenic allele): 2 * p * q
# expected_genotype_rate <- p^2 + 2 * p * q
# 
# # Expected number of CF-causing genotypes (if fully expressed) in the general population.
# expected_R117 <- population_size * expected_genotype_rate
# 
# # Incorporate survival prior to registry inclusion.
# # Annual mortality among CF patients is ~0.4%. Although the registry reports a median age of 22,
# # we can approximate cumulative mortality up to that age.
# annual_mortality_rate <- 0.004
# average_age <- 22
# # Cumulative mortality using an exponential survival model.
# cumulative_mortality <- 1 - exp(-annual_mortality_rate * average_age)
# survival_factor <- 1 - cumulative_mortality
# 
# # Adjust the expected number by survival (i.e. those surviving to be captured in the registry).
# adjusted_expected_R117 <- expected_R117 * survival_factor
# 
# print(data.frame(
#   HGVSp_VEP = df_arg117his$HGVSp_VEP,
#   gnomAD_genomes_AF = p,
#   total_AF = total_AF,
#   expected_hom = population_size * p^2,
#   expected_comphet = population_size * 2 * p * q,
#   expected_R117 = expected_R117,
#   cumulative_mortality = cumulative_mortality,
#   survival_factor = survival_factor,
#   adjusted_expected_R117 = adjusted_expected_R117
# ))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # For CFTR, restrict to known pathogenic variants with nonzero allele counts using gnomAD_genomes_AF.
# cftr_patho <- df %>% 
#   filter(genename == "CFTR", clinvar_clnsig %in% "Pathogenic", gnomAD_genomes_AN > 0)
# 
# # Total pathogenic allele frequency for CFTR (using gnomAD_genomes_AF).
# total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)
# 
# # Focus on p.Arg117His.
# df_arg117his <- cftr_patho %>% 
#   filter(HGVSp_VEP == "p.Arg117His")
# 
# # p = frequency of p.Arg117His; q = frequency of all other pathogenic alleles.
# p <- df_arg117his$gnomAD_genomes_AF
# q <- total_AF - p
# 
# # Expected genotype counts from Hardy–Weinberg in the general population:
# expected_hom <- population_size * p^2             # homozygous for p.Arg117His
# expected_comphet <- population_size * 2 * p * q      # compound heterozygotes with p.Arg117His
# 
# # Sum of expected genotype counts (assuming full penetrance):
# expected_genotypes <- expected_hom + expected_comphet
# 
# # Adjust for mortality prior to registry (e.g. ~0.4% loss).
# mortality_rate <- 0.004
# expected_after_mortality <- expected_genotypes * (1 - mortality_rate)
# 
# # p.Arg117His is known to have reduced penetrance.
# # Observed in the registry: 714 cases among 11318 total CF cases (~6.3%).
# observed_R117 <- 714
# # Estimate the penetrance factor as the ratio of observed to expected (after mortality).
# penetrance_R117 <- observed_R117 / expected_after_mortality
# 
# # Adjust the expected count to account for reduced penetrance.
# expected_R117_adjusted <- expected_after_mortality * penetrance_R117
# 
# print(data.frame(
#   HGVSp_VEP = df_arg117his$HGVSp_VEP,
#   gnomAD_genomes_AF = p,
#   expected_hom = expected_hom,
#   expected_comphet = expected_comphet,
#   expected_after_mortality = expected_after_mortality,
#   penetrance_factor = penetrance_R117,
#   expected_R117_adjusted = expected_R117_adjusted
# ))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 


# new ----

# library(ggplot2); theme_set(theme_bw())
# library(patchwork)
# library(dplyr)
# library(stringr)
# library(tidyr)



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


# Ensure the output directory exists
if(!dir.exists("../images/")) dir.create("../images/")

# UK population
population_size <- 69433632
# population_size <- 83702  # births in 2023

# DBNSFP Data Import and Filtering ----
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
head(df, 1)

# Project specific population ----
df <- df |> select(-gnomAD_genomes_AF)
colnames(df)[colnames(df) == 'gnomAD_genomes_NFE_AF'] <- 'gnomAD_genomes_AF'

# Keep a copy and filter out rows with clinvar_clnsig == "."
df <- df %>% filter(!clinvar_clnsig == ".")
df |> count(clinvar_clnsig)

# Bar plot: Count of ClinVar Clinical Significance
p_count <- ggplot(df, aes(
  x = str_wrap(paste0(gsub("_", " ", clinvar_clnsig)), width = 20), 
  fill = clinvar_clnsig)) +
  geom_bar(color = "black") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Count",
       title = "Count of ClinVar Clinical Significance CFTR") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

p_count
ggsave("../images/cftr_clinvar_count.png", plot = p_count, width = 5, height = 4)

# Data Preparation for Population-Level Calculations ----
df$gnomAD_genomes_AF[df$gnomAD_genomes_AF == "."] <- 0
df$gnomAD_genomes_AN[df$gnomAD_genomes_AN == "."] <- 0
df$gnomAD_genomes_AF <- as.numeric(df$gnomAD_genomes_AF)
df$gnomAD_genomes_AN <- as.numeric(df$gnomAD_genomes_AN)
df <- df %>% select(genename, `pos(1-based)`, gnomAD_genomes_AN, gnomAD_genomes_AF, 
                    clinvar_clnsig, HGVSc_VEP, HGVSp_VEP)

# Keep just one transcript allele for simplicity
df$HGVSc_VEP <- sapply(strsplit(df$HGVSc_VEP, ";"), `[`, 2)  # cftr registry transcript
df$HGVSp_VEP <- sapply(strsplit(df$HGVSp_VEP, ";"), `[`, 2)  # cftr registry transcript
df$genename <- sapply(strsplit(df$genename, ";"), `[`, 1)
df$Inheritance <- "AR"  # setting inheritance to AR for autosomal recessive conditions

# Substitute zero allele frequency with a synthesized minimum allele frequency.
max_af <- max(df$gnomAD_genomes_AN, na.rm = TRUE)
df$gnomAD_genomes_AF <- ifelse(df$gnomAD_genomes_AF == 0,
                               1 / (max_af + 1),
                               df$gnomAD_genomes_AF)

# Update occurrence probability calculation ----
# First, get the total pathogenic allele frequency per gene.
# We filter for pathogenic variants (clinvar_clnsig == "Pathogenic") with nonzero allele numbers.
df_patho <- df %>% filter(gnomAD_genomes_AN > 0, clinvar_clnsig == "Pathogenic")
gene_patho <- df_patho %>%
  group_by(genename) %>%
  summarise(total_pathogenic_AF = sum(gnomAD_genomes_AF, na.rm = TRUE), .groups = "drop")

# Join total pathogenic allele frequency into the main dataframe.
df <- df %>% left_join(gene_patho, by = "genename")

# For AR conditions, an event (carriage) occurs when any two alleles are present.
# This can be due to homozygous variants: (AF)^2, or compound heterozygous variants:
# 2 * AF * (total_pathogenic_AF - AF).
# For AD and X-linked, we simply use AF.
# We call this the "occurrence probability" (always ≥ 0).
df <- df %>% mutate(
  other_AF = ifelse(Inheritance == "AR", total_pathogenic_AF - gnomAD_genomes_AF, 0),
  occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"),
                           gnomAD_genomes_AF,
                           pmax(gnomAD_genomes_AF^2 + 2 * gnomAD_genomes_AF * other_AF, 0))
)

# Calculate expected cases and probability of at least one event in the population.
df <- df %>% mutate(
  expected_cases = population_size * occurrence_prob,
  prob_at_least_one = 1 - (1 - occurrence_prob)^population_size
)

# View the recalculated key results.
head(df)

# Scatter Plots: Expected Cases and Probability vs Allele Frequency ----

# Compute threshold allele frequency at which probability is (almost) 1 (e.g., ≥ 0.999)
threshold_AF <- min(df$gnomAD_genomes_AF[df$prob_at_least_one >= 0.999])
threshold_AF_label <- formatC(threshold_AF, format = "f", digits = 6)

unique_labels <- df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>%
  distinct(gnomAD_genomes_AF, expected_cases)

p_scatter1_path <- 
  df %>% filter(clinvar_clnsig == "Pathogenic") %>%
  ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  # scale_x_log10() +
  ylim(0, 2400) +  # something sensible if only one value
  geom_text(data = unique_labels, 
            aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Expected Cases",
       title = "Expected Cases vs\nAllele Frequency CFTR",
       subtitle = paste0("Condition: population size ", population_size))

vline_label <- tibble(threshold_AF = threshold_AF, 
                      threshold_AF_label = threshold_AF_label)

p_scatter2_path <- 
  df %>% 
  filter(clinvar_clnsig == "Pathogenic") %>% 
  ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, colour = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  labs(x = "Allele Frequency (gnomAD_genomes_AF, log scale)",
       y = "Probability of ≥1 Event",
       title = "Probability of At Least One\nEvent vs Allele Frequency CFTR",
       subtitle = paste0("Condition: population size ", population_size)) +
  geom_vline(xintercept = threshold_AF, linetype = "dotted", colour = "black") +
  geom_text(data = vline_label, 
            aes(x = threshold_AF, y = 1, label = threshold_AF_label), 
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  ylim(0, 1.2)

# Display the scatter plots.
p_scatter1_path
p_scatter2_path

# Combine and save scatter plots vertically.
p_scatter1_path <- p_scatter1_path + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_scatter <- p_scatter1_path / p_scatter2_path + 
  plot_layout(guides = 'collect', axis = "collect")  + 
  plot_annotation(tag_levels = 'A')
print(p_scatter)

# Density histograms for Expected Cases by ClinVar Clinical Significance.
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

p_scatter_dense <- (p_density / (p_scatter1_path + p_scatter2_path)) +
  plot_layout(widths = c(1, 1), guides = 'collect', axis = "collect") +
  plot_annotation(tag_levels = 'A')
print(p_scatter_dense)
ggsave("../images/cftr_scatterdense_expected_prob.png", plot = p_scatter_dense, width = 12, height = 6)

# Tally by ClinVar Category ----
df_calc <- df %>%
  mutate(allele_freq = as.numeric(gnomAD_genomes_AF)) %>%
  filter(!is.na(allele_freq)) %>%
  mutate(occurrence_prob = ifelse(Inheritance %in% c("AD", "X-linked"), 
                                  allele_freq,
                                  pmax(allele_freq^2 + 2 * allele_freq * (total_pathogenic_AF - allele_freq), 0)),
         expected_cases = population_size * occurrence_prob,
         prob_at_least_one = 1 - (1 - occurrence_prob)^population_size)

clinvar_levels <- unique(df$clinvar_clnsig)

df_tally <- df_calc %>%
  group_by(clinvar_clnsig) %>%
  summarise(total_expected_cases = sum(expected_cases),
            overall_prob = 1 - prod(1 - occurrence_prob),
            .groups = "drop") %>%
  complete(clinvar_clnsig = clinvar_levels,
           fill = list(total_expected_cases = 0, overall_prob = 0))
print(df_tally)

# Bar plot: Total Expected Cases by ClinVar Category.
p_bar <- ggplot(df_tally, aes(
  x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = total_expected_cases, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = formatC(round(total_expected_cases, 0), format = "d", big.mark = ",")), 
            vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Total Expected\nCases",
       title = "Total Expected Cases by ClinVar\nClinical Significance in CFTR",
       subtitle = paste0("Condition: population size ", population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

p_bar

# Bar plot: Overall Probability by ClinVar Category.
p_prob <- ggplot(df_tally, aes(
  x = str_wrap(gsub("_", " ", clinvar_clnsig), width = 20), 
  y = overall_prob, fill = clinvar_clnsig)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(overall_prob, accuracy = 0.001)), vjust = -0.5, size = 3.5) +
  labs(x = "ClinVar Clinical Significance",
       y = "Overall\nProbability",
       title = "Overall Probability of an Affected\nBirth by ClinVar Category in CFTR",
       subtitle = paste0("Condition: population size ", population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::percent)

p_prob

# Combine bar charts using patchwork and save.
p_bar <- p_bar + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
print(p_bars)
ggsave("../images/CFTR_combined_bar_charts.png", plot = p_bars, width = 6, height = 6)


# check expected disease ----
# 
# 
# # For CFTR, restrict to known pathogenic variants using gnomAD_genomes_AF.
# cftr_patho <- df %>% 
#   filter(genename == "CFTR", clinvar_clnsig == "Pathogenic", gnomAD_genomes_AN > 0)
# 
# # Total pathogenic allele frequency for CFTR.
# total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)
# 
# # Focus on p.Arg117His.
# df_arg117his <- cftr_patho %>% 
#   filter(HGVSp_VEP == "p.Arg117His")
# 
# # p = allele frequency for p.Arg117His; 
# # q = sum of all other pathogenic allele frequencies in CFTR.
# p <- df_arg117his$gnomAD_genomes_AF
# q <- total_AF - p
# 
# # Expected genotype counts in the general population using Hardy–Weinberg:
# expected_hom <- population_size * p^2             # homozygous for p.Arg117His
# expected_comphet <- population_size * 2 * p * q      # compound heterozygotes carrying p.Arg117His with any other pathogenic allele
# expected_total_genotypes <- expected_hom + expected_comphet
# 
# # Optionally, adjust for a simple mortality rate (e.g. 0.4% loss)
# mortality_rate <- 0.004
# expected_after_mortality <- expected_total_genotypes * (1 - mortality_rate)
# 
# print(data.frame(
#   HGVSp_VEP = df_arg117his$HGVSp_VEP,
#   gnomAD_genomes_AF = p,
#   expected_hom = expected_hom,
#   expected_comphet = expected_comphet,
#   expected_total_genotypes = expected_total_genotypes,
#   expected_after_mortality = expected_after_mortality
# ))
# 
# 
# # p.Arg117His is known to have reduced penetrance.
# # Observed in the registry: 714 cases among 11318 total CF cases (~6.3%).
# observed_R117 <- 714
# source_text <- "Source: UK Cystic Fibrosis Registry 2023 Annual Data Report October 2024\nObserved: 714/11318 (~6.3%)"
# 
# # Create a summary data frame with the calculated expected counts and the observed count.
# summary_df <- data.frame(
#   Category = c("Expected\nHomozygous", "Expected\nCompound Het", "Expected\nTotal Genotypes", "Expected\nAfter Mortality", "Observed *"),
#   Count = c(expected_hom, expected_comphet, expected_total_genotypes, expected_after_mortality, observed_R117)
# )
# 
# # Plot the summary as a bar chart with annotations and add source text.
# ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
#   geom_bar(stat = "identity", color = "black") +
#   geom_text(aes(label = formatC(round(Count, 0), format = "d", big.mark = ",")),
#             vjust = -0.5, size = 4) +
#   labs(title = "Expected Genotype Counts for p.Arg117His in CFTR",
#        subtitle = paste0("Population Size: ", population_size),
#        caption = source_text,
#        x = "Genotype Category",
#        y = "Count") +
#   guides(fill = "none") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
# 
# 
# 
# 
# 

# i think it is important that our number of cases also relates to Total deaths reported during annual review year (%): 49 (0.4%)
# Age at death in years; median (95% CI): 46 (37, 55)
# Age in years; median: 22. 

# 
# 
# # Simulate the effect of mortality on the observed genotype counts.
# # For example, we can assume an exponential survival model.
# # Annual mortality among CF patients is ~0.4% (0.004), and the median age of cases is 22.
# annual_mortality_rate <- 0.004
# median_age <- 22
# # Cumulative survival probability (i.e. fraction surviving to age 22):
# survival_factor <- exp(-annual_mortality_rate * median_age)
# 
# # Adjust the expected total genotypes by the survival factor.
# adjusted_expected_after_mortality <- expected_total_genotypes * survival_factor
# 
# # Create a summary data frame that includes both the unadjusted and adjusted expected counts,
# # along with the observed count.
# summary_df <- data.frame(
#   Category = c("Homozygous", "Compound Het", "Total Genotypes", 
#                "After Mortality (Unadjusted)", "After Mortality (Adjusted)", "Observed"),
#   Count = c(expected_hom, expected_comphet, expected_total_genotypes,
#             expected_total_genotypes * (1 - 0.004),  # previous simple mortality adjustment (0.4%)
#             adjusted_expected_after_mortality,
#             observed_R117)
# )
# 
# # Plot the summary as a bar chart with annotations and include the source.
# source_text <- "Source: UK Cystic Fibrosis Registry 2023 Annual Data Report October 2024\nObserved: 714/11318 (~6.3%)"
# ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
#   geom_bar(stat = "identity") +
#   geom_text(aes(label = formatC(round(Count, 0), format = "d", big.mark = ",")), 
#             vjust = -0.5, size = 4) +
#   labs(title = "Expected Genotype Counts for p.Arg117His in CFTR",
#        subtitle = paste0("Population Size: ", population_size, "\nSurvival Factor (to median age 22): ", 
#                          formatC(survival_factor, format = "f", digits = 3)),
#        caption = source_text,
#        x = "Genotype Category",
#        y = "Count") +
#   guides(fill = "none") +
#   theme_minimal()













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

ggsave("../images/cftr_validation_pArg117His.png", plot = p_var, width = 6, height = 6)


