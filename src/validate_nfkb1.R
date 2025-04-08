library(ggplot2); theme_set(theme_bw())

# Skip one of the confitions because it is special case where we want a variant that is not listed in the clinsig
# KEPP_ALL_FOR_VALIDATION_SEARCH <- TRUE
rm(KEPP_ALL_FOR_VALIDATION_SEARCH)
source("./inheritance_prob_generalised_mini.R")

# For CFTR, restrict to known pathogenic variants using gnomAD_genomes_AF.
nfkb1_patho <- df %>%
  filter(genename == "NFKB1", clinvar_clnsig == "Pathogenic")

expected_total_genotypes <- nfkb1_patho$expected_cases 

# Set expected_total_genotypes_true based on synth_flag
# tehcnically this has to be 0 but lets do it properly anyway.
if(nfkb1_patho$synth_flag == TRUE){
  expected_total_genotypes_true <- 0
} else {
  expected_total_genotypes_true <- expected_total_genotypes
}

# From the study of NFKB1: Given that PID is a heterogeneous disease, with overlap in phenotypes and genetic causes across different diagnostic categories, we decided to perform an unbiased genetic analysis of all 846 unrelated index cases. Whole-genome sequence data were additionally available for 63 affected and 345 unaffected relatives. Within a broad range of phenotypes, CVID is the most common disease category, comprising 46% of the NIHRBR-RD PID cohort (n = 390 index cases; range, 0-93 years of age).

# Baseline ----
# Revised validation calculations for NFKB1-related CVID cases

# # Cohort data from a specialized clinical setting (assumed to represent nearly all national PID cases)
# N_cohort <- 846       # Total number of PID patients in the cohort
# n_NFKB1  <- 390       # Number of CVID cases observed in the cohort
# 
# # Calculate the observed prevalence of NFKB1-related CVID among PID patients in the cohort
# prevalence_cohort <- n_NFKB1 / N_cohort
# cat("Observed cohort prevalence of NFKB1-related CVID:", round(prevalence_cohort, 4), "\n")
# # This value (≈0.461) represents the proportion of PID cases that are NFKB1-related within the cohort.
# 
# # Literature-based approach:
# # CVID prevalence in the general population is ~1 in 25,000.
# prevalence_CVID <- 1 / 25000
# 
# # For a general population (e.g., the UK), calculate the expected number of CVID cases.
# population_UK <- 69433632
# expected_CVID_UK <- prevalence_CVID * population_UK
# cat("Expected number of CVID cases in the UK (from literature):", round(expected_CVID_UK), "\n")
# 
# # If we assume that the proportion of CVID cases attributable to NFKB1 is similar to the cohort value,
# # then the estimated number of NFKB1-related CVID cases in the UK is:
# estimated_NFKB1_UK <- expected_CVID_UK * prevalence_cohort
# cat("Estimated number of NFKB1-related CVID cases in the UK (extrapolated):", round(estimated_NFKB1_UK), "\n")
# 
# # Sensible estimate discussion:
# # Because our cohort is derived from a specialized clinical setting that likely captures most PID cases,
# # the observed count (390 NFKB1-related cases) may be a more realistic estimate of the national burden
# # for NFKB1-related CVID than the extrapolation from general population prevalence.
# cat("Alternatively, if the cohort represents nearly all PID cases in the region, the national estimate would be close to:", n_NFKB1, "\n")

# CI ----

# Load required packages
library(binom)
library(dplyr)

# Given values from the cohort study
N_cohort <- 846      # Total number of PID patients in the cohort
n_cvid  <- 390     # Number of CVID cases observed in the cohort
n_NFKB1  <- 16     # Number of CVID cases observed in the cohort


# Calculate the observed prevalence and its 95% confidence interval using Wilson's method
conf_int <- binom.confint(n_NFKB1, n_cvid, methods = "wilson")
p_point   <- conf_int$mean      # Point estimate 
p_lower   <- conf_int$lower     # Lower bound
p_upper   <- conf_int$upper     # Upper bound

cat("95% CI for cohort prevalence: (", round(p_lower, 4), ",", round(p_upper, 4), ")\n")

# Literature reports that the prevalence of CVID in the general population is approximately 1 in 25,000.
prevalence_CVID <- 1 / 25000

# UK population (example value)
population_UK <- 69433632

# Calculate the expected total number of CVID cases in the UK based on literature prevalence.
expected_CVID_UK <- prevalence_CVID * population_UK
cat("Expected number of CVID cases in the UK (from literature):", round(expected_CVID_UK), "\n")
# Expected_CVID_UK is approximately 2777 cases.

# Estimate the number of NFKB1-related CVID cases in the UK using the cohort prevalence
estimated_NFKB1_UK <- expected_CVID_UK * p_point
cat("Estimated number of NFKB1-related CVID cases in the UK (extrapolated):", round(estimated_NFKB1_UK), "\n")
# This gives a point estimate of about 1280 cases.

# Calculate the confidence interval for the estimated number of NFKB1-related cases
lower_estimate <- expected_CVID_UK * p_lower
upper_estimate <- expected_CVID_UK * p_upper
cat("95% CI for estimated NFKB1-related cases in the UK:", 
    round(lower_estimate), "to", round(upper_estimate), "\n")
# The 95% CI ranges approximately from 1188 to 1374 cases.

# Bayesian adjustment with uncertainty:
# We assume that 90% of the information comes from the specialized cohort and 10% from the literature-based extrapolation.

w <- 0.9  # weight for the cohort data

# Extrapolated estimates based on literature:
estimated_NFKB1_UK <- expected_CVID_UK * p_point  #  (point estimate)
lower_estimate <- expected_CVID_UK * p_lower      
upper_estimate <- expected_CVID_UK * p_upper      

# Bayesian adjusted estimate: weighted average of cohort count (n_NFKB1) and literature extrapolation
adjusted_estimate <- w * n_NFKB1 + (1 - w) * estimated_NFKB1_UK
adjusted_lower <- w * n_NFKB1 + (1 - w) * lower_estimate
adjusted_upper <- w * n_NFKB1 + (1 - w) * upper_estimate

cat("Bayesian adjusted estimated number of NFKB1-related CVID cases in the UK:", round(adjusted_estimate), "\n")
cat("95% Bayesian adjusted CI: (", round(adjusted_lower), ",", round(adjusted_upper), ")\n")


# plot ----

# Bayesian posterior sampling for cohort prevalence
set.seed(123)
nsamples <- 10000
alpha <- n_NFKB1 + 1     
beta_param <- n_cvid - n_NFKB1 + 1  
p_samples <- rbeta(nsamples, alpha, beta_param)

# Literature extrapolated estimates: using expected total CVID cases in the UK (≈2777)
estimated_samples <- expected_CVID_UK * p_samples  # estimated NFKB1-related cases from literature extrapolation

# Bayesian adjusted estimates: weighted average of the cohort count and literature extrapolation
w_new <- 0.5  # weight for the cohort data. We generate a mixture distribution (bayes_adjusted_mix_nfkb1) by setting the weight to 0.5, thereby allowing the final estimate to be influenced equally by the cohort data (constant 390) and the literature-based extrapolation (which is variable/conditional).
# Adjust the weight to produce a mixture distribution that better reflects uncertainty 
# from both the cohort and the literature-based extrapolation.


# Recompute the Bayesian adjusted samples using the new weight:
bayes_adjusted_mix_nfkb1 <- w_new * n_NFKB1 + (1 - w_new) * estimated_samples

# Create a data frame for plotting
df_samples_nfkb1 <- data.frame(
  Estimated = estimated_samples,
  BayesianAdjusted = bayes_adjusted_mix_nfkb1
)

# Compute median and 95% CI for the Bayesian adjusted distribution
median_mix <- median(bayes_adjusted_mix_nfkb1)
ci_mix <- quantile(bayes_adjusted_mix_nfkb1, probs = c(0.025, 0.975))
cat("Bayesian mixture adjusted estimate (w = 0.5) median:", round(median_mix, 0), "\n")
cat("95% CI for the Bayesian mixture adjusted estimate:", round(ci_mix[1], 0), "to", round(ci_mix[2], 0), "\n")
density_mix <- density(bayes_adjusted_mix_nfkb1)
max_y_mix <- max(density_mix$y)

# Compute median and 95% CI for the literature extrapolated distribution (Estimated)
median_est <- median(estimated_samples)
ci_est <- quantile(estimated_samples, probs = c(0.025, 0.975))
cat("Literature extrapolated estimate median:", round(median_est, 0), "\n")
cat("95% CI for the literature extrapolated estimate:", round(ci_est[1], 0), "to", round(ci_est[2], 0), "\n")
density_est <- density(estimated_samples)
max_y_est <- max(density_est$y)

# Create a combined density plot with vertical dotted lines and annotations for medians and 95% CI
p_nfkb1_bayes <- ggplot() +
  geom_density(data = df_samples_nfkb1, aes(x = Estimated, fill = "Literature Extrapolated"), alpha = 0.5) +
  geom_density(aes(x = bayes_adjusted_mix_nfkb1, fill = "Bayesian Mixture Adjusted"), alpha = 0.5) +
  scale_fill_manual(name = "Estimate Type", 
                    values = c("Literature Extrapolated" = "skyblue", 
                               "Bayesian Mixture Adjusted" = "orange")) +
  # labs(title = "Bayesian Adjusted Estimates for NFKB1-related\nCVID Case Genotype Counts",
  labs(title = "Autosomal Dominant: NFKB1-related all CVID Cases",
       x = "Estimated Number of Cases",
       y = "Density") +
  # Extend y-axis by 1.2 times to allow space for annotations
  scale_y_continuous(limits = c(0, max_y_mix * 1.2)) +
  
  # Predicted total genotypes (or 0 if synth_flag is TRUE)
  geom_vline(xintercept = expected_total_genotypes_true, linetype = "solid", color = "darkgreen", size = 1) +
  annotate("label", x = expected_total_genotypes_true, max_y_est * 2.2, 
           label = paste0("Predicted total = ", round(expected_total_genotypes_true)), 
           fill = scales::alpha("white", 0.5), label.size = 0,
           color = "darkgreen", vjust = 0, size = 4, hjust = 0) +

  geom_vline(xintercept = n_NFKB1 , linetype = "dotted", color = "red", size = 1) +
  annotate("label", x = median_est *1 , y = max_y_est * 1.7, 
           label = paste("Reported\ncases =", formatC(n_NFKB1, format = "f", digits = 0)), 
           fill = scales::alpha("white", 0.5),  label.size = 0,
           color = "red", vjust = 0, size = 4,  hjust = 0) +
  
  geom_vline(xintercept = median_mix, linetype = "dotted", color = "orange", size = 1) +
  annotate("text", x = median_est *1.4, y = max_y_est * 1, 
           label = paste("Bayesian\nMedian =", formatC(median_mix, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_mix[1], format = "f", digits = 0), ",", 
                         formatC(ci_mix[2], format = "f", digits = 0), "]"), 
           color = "orange", vjust = 0, size = 4, hjust = 0) +
  
  geom_vline(xintercept = median_est, linetype = "dotted", color = "blue", size = 1) +
  annotate("text", x = median_est *1.6, y = max_y_est * .2, 
           label = paste("Estimate Max\nMedian =", formatC(median_est, format = "f", digits = 0),
                         "\n95% CI: [", formatC(ci_est[1], format = "f", digits = 0), ",", 
                         formatC(ci_est[2], format = "f", digits = 0), "]"), 
           color = "blue", vjust = 0, size = 4, hjust = 0) +
  
  geom_vline(xintercept = expected_total_genotypes, linetype = "solid", color = "darkgreen", size = 1) +
  annotate("label", x = expected_total_genotypes *0.75, y = max_y_est * .1,
           label = paste0("Predicted\ntotal\nsynthetic = ", round(expected_total_genotypes)),
           fill = scales::alpha("white", 0.5),  label.size = 0,
           color = "darkgreen", vjust = 0, size = 4, hjust = 0)
  #xlim(200, 2200)

# print(p_nfkb1_bayes)
ggsave("../images/nfkb1_case_est_distribution_combined_mixture.png", plot = p_nfkb1_bayes, width = 9, height = 3)

