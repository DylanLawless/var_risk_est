# --- Bayesian Estimation for CFTR (p.Arg117His) Recessive Disease Example ---
source("./inheritance_prob_cftr.R")

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
bayes_adjusted_mix <- w_new * observed_R117 + (1 - w_new) * adjusted_samples

# Create a data frame for plotting the density distributions.
df_samples <- data.frame(
  LiteratureExtrapolated = adjusted_samples,
  BayesianAdjusted = bayes_adjusted_mix
)

# Compute medians and 95% credible intervals.
median_lit <- median(adjusted_samples)
ci_lit <- quantile(adjusted_samples, probs = c(0.025, 0.975))
median_bayes <- median(bayes_adjusted_mix)
ci_bayes <- quantile(bayes_adjusted_mix, probs = c(0.025, 0.975))

# Plot density distributions with vertical lines and annotations.
p_bayes <- ggplot() +
  geom_density(data = df_samples, aes(x = LiteratureExtrapolated, fill = "Literature Extrapolated"), alpha = 0.5) +
  geom_density(data = df_samples, aes(x = BayesianAdjusted, fill = "Bayesian Mixture Adjusted"), alpha = 0.5) +
  scale_fill_manual(name = "Estimate Type",
                    values = c("Literature Extrapolated" = "skyblue",
                               "Bayesian Mixture Adjusted" = "orange")) +
  labs(title = "Bayesian Adjusted Estimates for CF-related p.Arg117His\nGenotype Counts (CFTR)",
       x = "Estimated Number of Cases",
       y = "Density") +
  geom_vline(xintercept = median_lit, linetype = "dotted", color = "blue", size = 1) +
  annotate("text", x = median_lit, y = max(density(adjusted_samples)$y) * 1.2,
           label = paste("Estimate\nMedian =", formatC(median_lit, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_lit[1], format = "f", digits = 0), ",",
                         formatC(ci_lit[2], format = "f", digits = 0), "]"),
           color = "black", size = 4, hjust = .5) +
  geom_vline(xintercept = median_bayes, linetype = "dotted", color = "orange", size = 1) +
  annotate("text", x = median_bayes, y = max(density(bayes_adjusted_mix)$y) * 0.8,
           label = paste("Bayesian\nMedian =", formatC(median_bayes, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_bayes[1], format = "f", digits = 0), ",",
                         formatC(ci_bayes[2], format = "f", digits = 0), "]"),
           color = "black", size = 4, hjust = 0.5) +
  geom_vline(xintercept = observed_R117, linetype = "dotted", color = "red", size = 1) +
  annotate("text", x = observed_R117, y = max(density(adjusted_samples)$y) * 1,
           label = paste("Reported\ncases =", observed_R117),
           color = "black", size = 4, hjust = 0.5) +
  xlim(200, 2200)  
  # 1641.548


print(p_bayes)
ggsave("../images/cftr_bayesian_adjusted_estimates.png", plot = p_bayes, width = 8, height = 4)
