library(patchwork)

# get validation result plots
# Be careful here as overlapping environment variables could get confused 
# Likewise, note that variables would be resued so see individual scripts for datasets

# dominant example
plot_nfkb1 <- local({
  source("validate_nfkb1.R")
  p_nfkb1_bayes
})

# recessive example
plot_cftr <- local({
  source("validate_cftr.R")
  p_cftr_bayes
})

patch1 <- plot_nfkb1 / plot_cftr
patch1 <- patch1  + plot_layout(guides = 'collect', axis = "collect")  + plot_annotation(tag_levels = 'A')
patch1

ggsave("../images/validation_studies_bayesian_adjusted_estimates.png", plot = patch1, width = 9, height = 6)




# An interpretation of results for these genes of interst ----

# Compute unique labels per gene for pathogenic entries
unique_labels <- df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  distinct(genename, gnomAD_genomes_AF, expected_cases)

p_scatter1_path <-
  df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  ggplot(aes(x = gnomAD_genomes_AF, y = expected_cases, color = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  ggrepel::geom_text_repel(data = unique_labels,
                           aes(x = gnomAD_genomes_AF, y = expected_cases, label = round(expected_cases)),
                           vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  labs(x = "log(Allele Frequency)",
       y = "Expected\nCases",
       # title = "Expected Cases vs\nAllele Frequency CFTR",
       # subtitle = paste0("Condition: population size ", population_size)
       ) +
  facet_wrap(~ genename, scales = "free") + 
  theme(legend.position = 'top') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_scatter1_path

# Compute threshold AF per gene for pathogenic entries
thresholds <- df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  group_by(genename) %>%
  summarise(threshold_AF = min(gnomAD_genomes_AF[prob_at_least_one >= 0.999], na.rm = TRUE),
            .groups = "drop") %>%
  mutate(threshold_AF_label = formatC(threshold_AF, format = "f", digits = 6))

df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  filter(genename == "CFTR") %>%
  summarise(threshold_AF = min(gnomAD_genomes_AF[prob_at_least_one >= 0.999], na.rm = TRUE),
            .groups = "drop")

df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  filter(genename == "NFKB1") %>%
  summarise(threshold_AF = min(gnomAD_genomes_AF[prob_at_least_one >= 0.999], na.rm = TRUE),
            .groups = "drop")

p_scatter2_path <-
  df %>%
  filter(clinvar_clnsig == "Pathogenic") %>%
  ggplot(aes(x = gnomAD_genomes_AF, y = prob_at_least_one, colour = clinvar_clnsig)) +
  geom_point(size = 3) +
  geom_line(aes(group = 1)) +
  labs(x = "log(Allele Frequency)",
       y = "AF for Probability\nof ≥1 Event",
       # title = "Probability of At Least One\nEvent vs Allele Frequency CFTR",
       # subtitle = paste0("Condition: population size ", population_size)
       ) +
  geom_vline(data = thresholds, aes(xintercept = threshold_AF),
             linetype = "dotted", colour = "black") +
  geom_text(data = thresholds,
            aes(x = threshold_AF, y = 1, label = threshold_AF_label),
            vjust = -1, hjust = 0.5, colour = "black", size = 3) +
  facet_wrap(~ genename, scales = "free") + 
  theme(legend.position = 'top') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_scatter2_path

p_density <- ggplot(df, aes(x = expected_cases, fill = clinvar_clnsig)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ clinvar_clnsig, scales = "free", ncol = 4,
             labeller = labeller(clinvar_clnsig = function(x) 
               stringr::str_wrap(gsub("_", " ", x), width = 20))) +
  scale_x_continuous(labels = scales::comma) +
  guides(fill = "none") +
  labs(x = "Expected Cases",
       y = "Density") +
       #        # title = "Density of Expected Cases by ClinVar Clinical Significance",
       #        # subtitle = paste0("Condition: population size ", population_size)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_density)

# Combine and save scatter plots vertically.
p_scatter1_path <- p_scatter1_path +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_scatter <- p_scatter1_path / p_scatter2_path +
  plot_layout(guides = 'collect', axis = "collect")  +
  plot_annotation(tag_levels = 'A') & theme(legend.position = 'top') 

print(p_scatter)

p_scatter_dense <- (p_density | p_scatter) +
  plot_layout(widths = c(2, 1), heights = c(1, 1), guides = 'collect', axis = "collect")  +
  plot_annotation(tag_levels = 'A', subtitle = paste0("Condition: population size ", population_size, ", phenotype PID-related, genes CFTR and NFKB1.")) 
print(p_scatter_dense)

ggsave("../images/validation_studies_scatterdense_expected_prob.png", plot = p_scatter_dense, width = 10, height = 6)


