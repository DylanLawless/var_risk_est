library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(forcats)
library(tidytext)
library(knitr)
library(ggplot2)
library(ggridges)
library(patchwork)
theme_set(theme_bw())

set.seed(666)

# paths and ids ---------------------------------------------------------------
PANEL          <- 398
genetic_defect <- "NFKB1"
path_data      <- "~/web/PanelAppRex/data"
df_core        <- readRDS(file.path(path_data, "path_PanelAppData_genes_combined_Rds")) |>
  rename(`Genetic defect` = entity_name) |>
  filter(panel_id == PANEL)

varRisEst_gene <- readRDS(sprintf("../output/VarRiskEst_PanelAppRex_ID_%d_gene_tally.Rds",  PANEL))
varRisEst_var  <- readRDS(sprintf("../output/VarRiskEst_PanelAppRex_ID_%d_gene_variants.Rds", PANEL))

# score map -------------------------------------------------------------------
score_map <- c(
  "Pathogenic"                                  = 5,
  "Likely pathogenic"                           = 4,
  "Pathogenic, low penetrance"                  = 3,
  "likely pathogenic, low penetrance"           = 3,
  "Conflicting classifications of pathogenicity" = 2,
  "risk factor"                                 = 1,
  "association"                                 = 1,
  "likely risk allele"                          = 1,
  "drug response"                               = 0,
  "Uncertain significance"                      = 0,
  "no classification for the single variant"    = 0,
  "no classifications from unflagged records"   = 0,
  "Affects"                                     = 0,
  "other"                                       = 0,
  "not provided"                                = 0,
  "uncertain risk allele"                       = 0,
  "protective"                                  = -3,
  "Likely benign"                               = -4,
  "Benign"                                      = -5
)

get_score <- function(clinvar) {
  terms  <- str_split(str_replace_all(clinvar, "_", " "), "[/|]", simplify = TRUE) |> str_trim()
  scores <- ifelse(terms %in% names(score_map), score_map[terms], 0)
  mean(scores)
}

# helper functions ------------------------------------------------------------
prior_weight   <- function(s) ifelse(s > 0, s + 1, 1)
posterior_mean <- function(a, b) a / (a + b)

# score and evidence ----------------------------------------------------------
varRisEst_gene_scored <- varRisEst_gene |>
  filter(genename == genetic_defect) |>
  mutate(score = map_dbl(clinvar_clnsig, get_score))

varRisEst_var_scored <- varRisEst_var |>
  filter(genename == genetic_defect) |>
  mutate(score          = map_dbl(clinvar_clnsig, get_score),
         evidence_score = score)

# patient data ----------------------------------------------------------------
patient_data <- varRisEst_var_scored |>
  mutate(
    known_var_observed = case_when(
      HGVSp_VEP == "p.Ser237Ter" ~  1L,
      # HGVSc_VEP == "c.159+1G>A"  ~ -9L,
      HGVSp_VEP == "p.Arg231His" ~  1L,
      HGVSp_VEP == "p.Gly650Arg" ~  1L,
      HGVSp_VEP == "p.Val236Ile" ~ -9L,
      HGVSp_VEP == "p.Thr567Ile" ~ -9L,
      TRUE                       ~  0L
    ),
    pathogenic_weight   = rescale(score, to = c(0, 1), from = c(-5, 5)),
    adj_occurrence_prob = occurrence_prob * pathogenic_weight,
    prior_w             = prior_weight(score)
  )

# cohort data - all REF alleles with 10% missing ----
cohort_size <- 200

patient_data <- patient_data %>%
  mutate(
    cohort_missing = map_int(
      known_var_observed,
      ~ {
        base_missing <- if (.x == -9L) 1L else 0L
        base_missing + rbinom(1, cohort_size - 1L, 0.10)
      }
    ),
    cohort_AN = 2L * (cohort_size - cohort_missing)
  ) %>%
  select(-cohort_missing)

# variants to model
risk_variants <- patient_data |>
  dplyr::filter((known_var_observed == 1 & score >= 0) | known_var_observed == -9) |>
  dplyr::mutate(
    flag        = dplyr::if_else(known_var_observed == -9, "missing", "present"),
    class       = dplyr::if_else(evidence_score > 3, "causal", "other"),
    variant_id  = dplyr::if_else(HGVSp_VEP == ".", HGVSc_VEP, HGVSp_VEP),
    # variant_lab = paste(flag, variant_id, sep = " "),
    variant_lab = paste(variant_id),
    alpha       = round(adj_occurrence_prob * cohort_AN) + prior_w,
    beta        = cohort_AN - round(adj_occurrence_prob * cohort_AN) + 1,
    group       = factor(paste(flag, class),
                         levels = c("present causal", "present other",
                                    "missing causal", "missing other"))
  )

# posterior sampling
n_sim   <- 10000
sim_mat <- sapply(seq_len(nrow(risk_variants)),
                  \(i) stats::rbeta(n_sim,
                                    risk_variants$alpha[i],
                                    risk_variants$beta[i]))
colnames(sim_mat) <- risk_variants$variant_lab

top_hit        <- apply(sim_mat, 1, which.max)
prob_is_topvec <- prop.table(table(factor(top_hit, levels = seq_len(nrow(risk_variants)))))

# posterior expectations & causal damaging share
posterior_share <- colMeans(sim_mat / rowSums(sim_mat))

risk_variants <- risk_variants |>
  dplyr::mutate(
    prob_is_top          = as.numeric(prob_is_topvec),
    posterior_share      = posterior_share[variant_lab],
    prob_causal_damaging = dplyr::if_else(class == "causal", posterior_share, 0)
  )

# variant‑level CI
ci_present <- risk_variants |>
  dplyr::filter(flag == "present") |>
  dplyr::mutate(
    lower  = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) stats::quantile(x, 0.025)),
    median = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) stats::quantile(x, 0.500)),
    upper  = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) stats::quantile(x, 0.975))
  ) |>
  dplyr::select(variant_lab, lower, median, upper)

plot_df <- risk_variants |>
  dplyr::left_join(ci_present, by = "variant_lab") |>
  dplyr::mutate(
    lower    = dplyr::if_else(flag == "present", lower, prob_is_top),
    median   = dplyr::if_else(flag == "present", median, prob_is_top),
    upper    = dplyr::if_else(flag == "present", upper, prob_is_top),
    var_plot = reorder_within(variant_lab, median, group)
  )

# plots -----------------------------------------------------------------------
p_points <- ggplot(plot_df,
                            aes(x = median, xmin = lower, xmax = upper,
                                         y = var_plot, fill = score)) +
  geom_pointrange(shape = 21, colour = "black", size = 0.8, stroke = 0.3) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
                      labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks  = c(-5, 0, 5),
                       limits = c(-5,5),
                       labels  = c("benign", "unknown", "pathogenic"),
                       na.value = "grey80") +
  scale_y_reordered() +
  labs(x = "probability", y = NULL, #fill = NULL
       ) +
  scale_x_continuous(limits = c(0, .5), breaks = c(0, 0.2, 0.4))

share_mat <- sim_mat / rowSums(sim_mat)

share_df <- tibble::as_tibble(share_mat) |>
  dplyr::mutate(draw = dplyr::row_number()) |>
  tidyr::pivot_longer(-draw, names_to = "variant_lab", values_to = "share") |>
  dplyr::left_join(risk_variants |>
                     dplyr::select(variant_lab, group, score), by = "variant_lab") |>
  dplyr::mutate(var_plot = reorder_within(variant_lab, share, group))

p_dist <- ggplot(share_df,
                          aes(x = share, y = var_plot, fill = score)) +
  ggridges::geom_density_ridges(scale = 1, rel_min_height = 0.01,
                                colour = "grey30", alpha = 0.7) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
                      labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                                breaks  = c(-5, 0, 5),
                                limits = c(-5,5),
                                labels  = c("benign", "unknown", "pathogenic"),
                                na.value = "grey80") +
  scale_y_reordered() +
  labs(x = "posterior probability", y = NULL) +
  # theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, .5), breaks = c(0, 0.2, 0.4))
p_dist

# new plot – per‑variant damaging‑causal probabilities
p_causal <- ggplot(
  risk_variants |>
    dplyr::mutate(var_plot = reorder_within(variant_lab,
                                                     prob_causal_damaging, group)),
  aes(x = prob_causal_damaging, y = var_plot, fill = score)
) +
  geom_point(shape = 21, colour = "black", size = 3) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
                      labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks  = c(-5, 0, 5),
                       limits = c(-5,5),
                       labels  = c("benign", "unknown", "pathogenic"),
                       na.value = "grey80") +
  scale_y_reordered() +
  labs(x = "p(causal & damaging)", y = NULL) +
  scale_x_continuous(limits = c(0, .5), breaks = c(0, 0.2, 0.4))
p_causal

# overall distribution
damaging_share_vec <- rowSums(
  share_mat[, risk_variants$class == "causal", drop = FALSE]
)

ci_overall <- stats::quantile(damaging_share_vec, c(0.025, 0.5, 0.975))

p_overall <- ggplot(tibble::tibble(share = damaging_share_vec),
                             aes(x = share)) +
  geom_histogram(bins = 50, fill = "red", alpha = 0.6, colour = "black") +
  geom_vline(xintercept = ci_overall[c(1, 3)], linetype = "dashed") +
  geom_vline(xintercept = ci_overall[2], colour = "purple") +
  annotate("text", x = ci_overall[2], 
           y = Inf,
           vjust = 1.2,
                    label = sprintf("median = %.3f\n95%% CI = [%.3f, %.3f]",
                                    ci_overall[2], ci_overall[1], ci_overall[3]),
                    hjust = 0.5, size = 4) +
  labs(x = "total probability damaging causal", y = "count") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.6)))

p_overall
# (p_points + p_dist + p_causal) / p_overall +
#   patchwork::plot_layout(guides = "collect", axis = "collect") +
#   patchwork::plot_annotation(tag_levels = "A")

# ────────────────────────────────────────────────────────────────────────────
# tables
summary_causal_damaging <- tibble::tibble(
  metric = "p(causal_damaging)",
  lower  = ci_overall[1],
  median = ci_overall[2],
  upper  = ci_overall[3]
)

variant_table_causal <- risk_variants |>
  dplyr::left_join(ci_present, by = "variant_lab") |>
  dplyr::mutate(
    lower  = dplyr::coalesce(lower, prob_is_top),
    median = dplyr::coalesce(median, prob_is_top),
    upper  = dplyr::coalesce(upper, prob_is_top)
  ) |>
  dplyr::arrange(group, dplyr::desc(prob_causal_damaging)) |>
  dplyr::select(group,
                variant              = variant_lab,
                evidence_score,
                lower, median, upper,
                posterior_share,
                prob_causal_damaging)

# console view
knitr::kable(summary_causal_damaging, digits = 3,
             caption = "overall probability of a damaging causal variant (95 % CI)")

knitr::kable(variant_table_causal, digits = 3,
             caption = "per‑variant probabilities with 95 % credible intervals")

# LaTeX view
latex_summary <- knitr::kable(summary_causal_damaging, digits = 3,
                              caption = "overall probability of a damaging causal variant (95\\% CI)",
                              format = "latex", booktabs = TRUE)

latex_variants <- knitr::kable(variant_table_causal, digits = 3,
                               caption = "per‑variant probabilities with 95\\% credible intervals",
                               format = "latex", booktabs = TRUE)

cat(latex_summary, "\n\n")
cat(latex_variants, "\n")




















# ---------------------------------------------------------------------------
# per‑variant damaging‑causal distributions (plot C)
share_df_causal <- share_df %>% 
  filter(risk_variants$class[match(variant_lab, risk_variants$variant_lab)] == "causal")

p_causal_dist <- ggplot(share_df_causal,
                   aes(x = share,
                       y = var_plot,
                       # fill = score
                       fill = group
                       )) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01,
                      colour = "grey30", alpha = 0.7) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(group = label_wrap_gen(width = 9))) +
  theme(strip.text.y = element_text(size = 7)) +
  # scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
  #                      breaks  = c(-5, 0, 5),
  #                      labels  = c("benign", "unknown", "pathogenic"),
  #                      na.value = "grey80") +
  scale_fill_manual(values = c(`present causal` = "forestgreen", `missing causal` = "orange")) +
  scale_y_reordered() +
  labs(x = "posterior p(causal & damaging)", y = NULL) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1))
# p_causal_dist

# ---------------------------------------------------------------------------
# contribution of present vs missing causal variants (plot D)
# source_df <- risk_variants %>% 
#   group_by(flag) %>% 
#   summarise(p_damaging = sum(prob_causal_damaging), .groups = "drop") %>% 
#   mutate(prop = p_damaging / sum(p_damaging))
# 
# p_source <- ggplot(source_df, aes(x = flag, y = p_damaging, fill = flag)) +
#   geom_col(width = 0.6, colour = "black") +
#   geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
#             vjust = -0.3, size = 3.5) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   scale_fill_manual(values = c(present = "forestgreen", missing = "orange")) +
#   labs(x = NULL,
#        y = "total\ncontribution\nprobability",
#        fill = NULL) +
#   theme_bw()



# ────────────────────────────────────────────────────────────────────────────
# stacked contribution of individual variants (plot D‑stack)

# keep every variant; an epsilon is added only for plotting so that ‘0’ variants
# still get a sliver and appear in the legend
eps <- 1e-2
variant_source_df <- risk_variants %>% 
  select(flag, score, prob_causal_damaging) %>% 
  mutate(prob_plot = prob_causal_damaging + eps)

source_totals <- variant_source_df %>% 
  group_by(flag) %>% 
  summarise(p_damaging = sum(prob_causal_damaging), .groups = "drop") %>% 
  mutate(prop = p_damaging / sum(p_damaging))

p_source_stack <- ggplot(variant_source_df,
                         aes(x = flag, y = prob_plot, fill = score)) +
  geom_col(width = 0.6, colour = "black", position = position_stack()) +
  geom_text(data = source_totals,
            aes(x = flag, y = p_damaging,
                label = scales::percent(prop, accuracy = 0.1)),
            vjust = -2, size = 3.5, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_fill_gradientn(colours = c("navy", "lightblue", "red"),
                       breaks  = c(-5, 0, 5),
                       limits  = c(-5, 5),
                       labels  = c("benign", "unknown", "pathogenic"),
                       na.value = "grey80") +
  labs(x = NULL,
       y = "total\np(causal & damaging)",
       # fill = NULL
       ) +
  theme_bw()
p_source_stack


# ---------------------------------------------------------------------------
# overall damaging‑causal distribution (plot E, unchanged)
# (code for p_overall already defined earlier)

# combined display: A = p_points, B = p_dist, C = p_causal, D = p_source, E = p_overall
# (p_points + p_dist + p_causal + p_source) / p_overall +
  # patchwork::plot_layout(guides = "collect", axis = "collect") +
  # patchwork::plot_annotation(tag_levels = "A")
  # 

p_quant <- (p_points + p_dist + p_causal) /  (p_source_stack + p_causal_dist) / p_overall +
  patchwork::plot_layout(guides = "collect", axis = "collect") +
  patchwork::plot_annotation(tag_levels = "A")  + 
  plot_layout(heights = c(2, 1, 1))

p_quant
ggsave("../images/plot_quant_uncert_ci.pdf", plot = p_quant, width = 9, height = 7)


# The results (Figure A–F) show that, of the total posterior probability that a damaging variant in NFKB1 underlies this proband’s phenotype, roughly half is carried by the single observed nonsense variant p.Ser237Ter and the other half is still attributable to known but unobserved pathogenic alleles at this locus.
# In panel A we plot each variant’s posterior share, with 95 % credible intervals: p.Ser237Ter stands out at about 0.24 (CI ~0.11–0.42), while the next highest is an unobserved splice donor variant c.159+1G>A at roughly 0.26 (no interval because it’s fixed by definition). Panel B’s ridge plots trace the full posterior distribution for every variant, illustrating the overlap and skewness in their shares. Panel C then isolates only the damaging calls (score > 3), confirming that p.Ser237Ter and the missing splice donor together account for almost all of the damaging‑causal probability.
# Panel D (the stacked bar) quantifies this directly: present variants contribute about 51.7 % of the damaging probability, and missing variants about 48.3 %. Panel E’s causal‑only density ridges further break down how each causal group (present vs missing) distributes its share. Finally, panel F shows the histogram of the total damaging‑causal probability (summing across all variants) with a median of 0.472 and a 95 % credible interval from 0.300 to 0.655. Together, these plots demonstrate that although p.Ser237Ter is the single most likely cause, nearly half of the aetiological probability remains tied to other known alleles not observed in this patient, underscoring the residual uncertainty and arguing for deeper or orthogonal validation.

# ---------------------------------------------------------------------------
# tabulating contributions
summary_source <- source_df %>% 
  select(flag, p_damaging, prop)

kable(summary_source, digits = 4,
      caption = "proportion of the total damaging‑causal probability attributable to present vs missing variants")

latex_source <- kable(summary_source, digits = 4,
                      caption = "proportion of the total damaging--causal probability attributable to present vs\\ missing variants",
                      format = "latex", booktabs = TRUE)

cat(latex_source, "\n")


# If you need the exact per‑variant numbers again
risk_variants %>%
  select(flag, variant_lab, evidence_score, posterior_share) %>%
  filter(evidence_score > 3) %>%           # damaging only
  arrange(desc(posterior_share))
# 
# That table will show the ~0.25 for the present variant and ~0.23 summed over the missing ones (each of those missing variants individually has a tiny share, but they add up).
