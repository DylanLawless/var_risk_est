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




# --- Bayesian Estimation for CFTR (p.Arg117His) Recessive Disease Example ---
# source("./inheritance_prob_cftr.R")

print("pos 117530975 (GrCh38 7-117530975-G-A, ENST00000003084.11 MANE Select transcript for CFTR HGVSpp.Arg117His) == p.Arg117His in transcript 2, and transcript 1 (p.Arg36His), but week keep trans 1 for simplicity across whole genome." )

# Skip one of the confitions because it is special case where we want a variant that is not listed in the clinsig
# KEPP_ALL_FOR_VALIDATION_SEARCH <- TRUE
rm(KEPP_ALL_FOR_VALIDATION_SEARCH)
source("./inheritance_prob_generalised_mini.R")
# source("./inheritance_prob_generalised_mini_v2.R")

# For CFTR, restrict to known pathogenic variants using gnomAD_genomes_AF.
cftr_patho <- df %>%
  filter(genename == "CFTR", clinvar_clnsig == "Pathogenic", gnomAD_genomes_AN > 0)

# Total pathogenic allele frequency for CFTR.
total_AF <- sum(cftr_patho$gnomAD_genomes_AF, na.rm = TRUE)

# Focus on p.Arg117His.
df_arg117his <- cftr_patho %>%
  # filter(HGVSp_VEP == "p.Arg117His")  # use transcript 2 if you want this aa coordinate
  filter(HGVSp_VEP == "p.Arg36His")  # use transcript 1 if you want this aa coordinate

# p = allele frequency for p.Arg117His;
# q = sum of all other pathogenic allele frequencies in CFTR.
p <- df_arg117his$gnomAD_genomes_AF
q <- total_AF - p

# Sanity check for p.Arg117His genotype expectations:
#
# For p.Arg117His, let:
#   p = allele frequency for p.Arg117His, and
#   q = sum of allele frequencies for all other pathogenic variants (q = total_AF - p).
#
# Under Hardy–Weinberg equilibrium, the full genotype probability for individuals carrying
# at least one copy of p.Arg117His is computed as:
#
#   p^2      -> probability of being homozygous for p.Arg117His, and
#   2 * p * q -> probability of being compound heterozygous (one copy is p.Arg117His and the other
#               is one of the other pathogenic variants).
#
# Thus, the total expected genotype probability becomes:
#
#   p^2 + 2 * p * q
#
# Multiplying this probability by the population size yields the expected number of individuals
# with any genotype involving p.Arg117His.
#
# Note: This full Hardy–Weinberg calculation (p^2 + 2 * p * q) is distinct from the partitioning
# strategy used to assign each variant an occurrence probability of p * total_AF (which equals p^2 + p*q).
# The partitioned value avoids double counting compound heterozygotes when summing across variants.
# However, for comparing with known disease cohort counts, the full genotype expectation (p^2 + 2*p*q)
# provides the appropriate theoretical estimate.

# Expected genotype counts in the general population using Hardy–Weinberg:
expected_hom <- population_size * p^2             # homozygous for p.Arg117His
expected_comphet <- population_size * 2 * p * q      # compound heterozygotes carrying p.Arg117His with any other pathogenic allele
expected_total_genotypes <- expected_hom + expected_comphet

# Note: Our population_size here is the total UK population, but the occurrence_prob value was 
# computed conditional on a patient having a specific phenotype (PID-related features). In this case, 
# the occurrence_prob of 0.00097506 represents the risk for a variant in that context.
#
# The sanity test compares the final reported occurrence_prob (conditioned on phenotype) 
# with the full Hardy–Weinberg (HW) expectation. Remember that our full HW calculation (p^2 + 2*p*q)
# would be calculated assuming the entire population is at risk.
#
# For p.Arg117His:
#   - p: allele frequency for p.Arg117His.
#   - q: contribution from all other pathogenic variants (q = total_AF - p).
#
# Under full HW expectations, the risk would be:
#   p^2         -> risk of being homozygous for p.Arg117His,
#   2 * p * q   -> risk of being compound heterozygous, where one allele is p.Arg117His.
#
# Multiplying this by population_size gives a full expected count based on the entire UK population.
#
# However, because our occurrence_prob value (here ~0.00097506) already incorporates a phenotype 
# restriction (i.e. a given patient has PID-related features), its corresponding expected count 
# should not be directly compared to the full HW expected count calculated for the whole population.

# Mortality adjustment using an exponential survival model:
# Annual mortality among CF patients ~0.4% (0.004) and median age is 22.
annual_mortality_rate <- 0.004
median_age <- 22
survival_factor <- exp(-annual_mortality_rate * median_age)
adjusted_expected_after_mortality <- expected_total_genotypes * survival_factor

# Observed data from the registry:
# p.Arg117His is known to have reduced penetrance.
# Observed in the registry: 714 cases among 11318 total CF cases (~6.3%).
observed_R117 <- 714
source_text <- "Source: UK Cystic Fibrosis Registry 2023 Annual Data Report, October 2024\nObserved: p.Arg117His in 714/11318 (~6.3%)"

# plot ----

expected_after_mortality_simple <- expected_total_genotypes * (1 - annual_mortality_rate)

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


# --- Bayesian Uncertainty Simulation ---

# We simulate uncertainty in p using a beta distribution.
# Here we use an effective allele count (pseudo-count) to reflect uncertainty in p.
effective_AN <- 1e6  # (chosen as a large number for illustration)
alpha_param <- p * effective_AN + 1
beta_param <- effective_AN - p * effective_AN + 1
nsamples <- 10000
p_samples <- rbeta(nsamples, alpha_param, beta_param)

# For each sampled p, compute the expected genotype count:
# (using the same formula: population_size * (p^2 + 2*p*(total_AF - p)))
expected_samples <- population_size * (p_samples^2 + 2 * p_samples * (total_AF - p))
# Adjust for mortality using the exponential survival model:
adjusted_samples <- expected_samples * survival_factor

# Bayesian mixture: combine the observed registry count and the literature-based extrapolation.
# Here, we use a 50:50 weight (w_new = 0.5) to reflect uncertainty from both sources.
w_new <- 0.5
bayes_adjusted_mix_cftr <- w_new * observed_R117 + (1 - w_new) * adjusted_samples

# Create a data frame for plotting the density distributions.
df_samples_cftr <- data.frame(
  LiteratureExtrapolated = adjusted_samples,
  BayesianAdjusted = bayes_adjusted_mix_cftr
)

# Compute medians and 95% credible intervals.
median_lit <- median(adjusted_samples)
ci_lit <- quantile(adjusted_samples, probs = c(0.025, 0.975))
median_bayes <- median(bayes_adjusted_mix_cftr)
ci_bayes <- quantile(bayes_adjusted_mix_cftr, probs = c(0.025, 0.975))

# Plot density distributions with vertical lines and annotations.
p_cftr_bayes <- ggplot() +
  geom_density(data = df_samples_cftr, aes(x = LiteratureExtrapolated, fill = "Literature Extrapolated"), alpha = 0.5) +
  geom_density(data = df_samples_cftr, aes(x = BayesianAdjusted, fill = "Bayesian Mixture Adjusted"), alpha = 0.5) +
  scale_fill_manual(name = "Estimate Type",
                    values = c("Literature Extrapolated" = "skyblue",
                               "Bayesian Mixture Adjusted" = "orange")) +
  labs(title = "Autosomal Recessive: CFTR-related p.Arg117His CF Cases",
    #title = "Bayesian Adjusted Estimates for CF-related p.Arg117His\nGenotype Counts (CFTR)",
       x = "Estimated Number of Cases",
       y = "Density") +
  
  geom_vline(xintercept = median_lit, linetype = "dotted", color = "blue", size = 1) +
  annotate("text", x = median_lit*1.1, y = max(density(adjusted_samples)$y) * .5,
           label = paste("Estimate\nMedian =", formatC(median_lit, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_lit[1], format = "f", digits = 0), ",",
                         formatC(ci_lit[2], format = "f", digits = 0), "]"),
           color = "blue", size = 4, hjust = 0) +
  
  geom_vline(xintercept = median_bayes, linetype = "dotted", color = "orange", size = 1) +
  annotate("text", x = median_lit *1.1, y = max(density(bayes_adjusted_mix_cftr)$y) * 0.6,
           label = paste("Bayesian\nMedian =", formatC(median_bayes, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_bayes[1], format = "f", digits = 0), ",",
                         formatC(ci_bayes[2], format = "f", digits = 0), "]"),
           color = "orange", size = 4, hjust = 0) +


  geom_vline(xintercept = observed_R117, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = observed_R117, y = max(density(adjusted_samples)$y) * 1.2,
           label = paste0("Reported    \ncases =", observed_R117),
           color = "red", size = 4, hjust = 1.5) +

  geom_vline(xintercept = expected_comphet, linetype = "solid", color = "darkgreen", size = 1) +
  annotate("text", x = expected_comphet, y =max(density(bayes_adjusted_mix_cftr)$y) * 0.9,
           label = paste0("Predicted Biallele = ", round(expected_comphet),
                         "\n(excl. Hom = ", round(expected_hom), ")"),
           color = "darkgreen", size = 4, hjust = 1) +
  
  geom_vline(xintercept = expected_total_genotypes, linetype = "solid", color = "darkgreen", size = 1) +
  annotate("text", x = median_lit*1.1, y =max(density(bayes_adjusted_mix_cftr)$y) * 0.8,
           label = paste0("Predicted total = ", round(expected_total_genotypes)),
           color = "darkgreen", size = 4, hjust = 0) +

  xlim(300, 1200)

print(p_cftr_bayes)

ggsave("../images/cftr_bayesian_adjusted_estimates.png", plot = p_cftr_bayes, width = 8, height = 4)































# 
# 
# 
# 
# 
# 
# # Compute unique labels per gene for pathogenic entries
# unique_labels <- cftr_patho %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   distinct(genename, gnomAD_genomes_AF, expected_cases)
# 
# p_scatter1_path <-
#   cftr_patho %>%
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
# thresholds <- cftr_patho %>%
#   filter(clinvar_clnsig == "Pathogenic") %>%
#   group_by(genename) %>%
#   summarise(threshold_AF = min(gnomAD_genomes_AF[prob_at_least_one >= 0.999], na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(threshold_AF_label = formatC(threshold_AF, format = "f", digits = 6))
# 
# p_scatter2_path <-
#   cftr_patho %>%
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