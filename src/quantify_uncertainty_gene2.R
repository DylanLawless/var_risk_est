library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(stringr)
library(patchwork)
library(ggplot2);theme_set(theme_bw())

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

# set example gene ----
genetic_defect <- "NFKB1"
varRisEst_gene <- varRisEst_gene |> filter(genename == genetic_defect)
varRisEst_var <- varRisEst_var |> filter(genename == genetic_defect)

# score classifications ----
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

varRisEst_var_scored <- varRisEst_var %>%
  mutate(score = map_dbl(clinvar_clnsig, get_score))

# make example patient data ----
# known variant from ClinVar was observed in our dataset: 0 - Reference allele in our data, 1 - heterozygous  ALT (matching clinvar) in our data. 2 - homozygous  ALT (matching clinvar) in our data. -9 is missing from our data meaning that we do not know the nucleotide genotype and it could be present for this patient. 

# single causal
patient_data <- varRisEst_var_scored |>
  filter(genename == "NFKB1") |>
  mutate(
    known_var_observed = case_when(
      HGVSp_VEP == "p.Ser237Ter" ~ 1L, # heterozygous known pathogenic
      HGVSp_VEP == "p.Val236Ile" ~ -9L, # missing
      HGVSp_VEP == "p.Thr567Ile" ~ -9L, # missing
      TRUE                         ~ 0L # all present
    )
  )



# multiple causal
patient_data <- varRisEst_var_scored |>
  filter(genename == "NFKB1") |>
  mutate(
    known_var_observed = case_when(
      HGVSp_VEP == "p.Ser237Ter" ~ 1L, # heterozygous known pathogenic
      HGVSc_VEP == "c.159+1G>A" ~ 1L, # het nown likely pathogenic
      HGVSp_VEP == "p.Val236Ile" ~ -9L, # missing
      HGVSp_VEP == "p.Thr567Ile" ~ -9L, # missing
      TRUE                         ~ 0L # all present
    )
  )


names(patient_data)

# reorder ----
patient_data <- patient_data |> select(
"genename",
"Inheritance",
"clinvar_clnsig",
"HGVSc_VEP",
"HGVSp_VEP",
"occurrence_prob",
"score",
"known_var_observed",
"gnomAD_genomes_AF",
"panel_id",
"name",
"max_an",
"synth_flag",
"total_AF")
patient_data


# calc ----
# ---- incorporate pathogenicity score into probabilities ----
patient_data <- patient_data %>% 
  mutate(
    pathogenic_weight   = scales::rescale(score, to = c(0, 1), from = c(-5, 5)),
    adj_occurrence_prob = occurrence_prob * pathogenic_weight
  )

# ---- recalculate probabilities with weighting ----
# probability any missing weighted variant is present
p_missing_any_w <- patient_data %>%
  filter(known_var_observed == -9) %>%
  pull(adj_occurrence_prob) %>%
  { 1 - prod(1 - .) }

# weighted candidate probability
p_candidate_w <- patient_data %>%
  filter(known_var_observed == 1) %>%
  pull(adj_occurrence_prob)

# ---- per‑score summary of missing variants ----
missing_score_summary <- patient_data %>% 
  filter(known_var_observed == -9) %>% 
  group_by(score) %>% 
  summarise(
    n_variants = n(),
    prob_any   = 1 - prod(1 - adj_occurrence_prob),
    .groups    = "drop"
  ) %>% 
  arrange(desc(score))

# ---- summary table (candidate vs missing, split by score) ----
total_prob <- p_candidate_w + p_missing_any_w

summary_df_w2 <- bind_rows(
  tibble(
    category    = "candidate variant (weighted)",
    probability = p_candidate_w,
    proportion  = p_candidate_w / total_prob
  ),
  missing_score_summary %>% 
    transmute(
      category    = paste0("missing variants score ", score, " (weighted)"),
      probability = prob_any,
      proportion  = prob_any / total_prob
    )
)

# ---- detailed variant‑level table ----
variant_table <- patient_data %>% 
  filter(known_var_observed == -9) %>% 
  arrange(desc(score), desc(adj_occurrence_prob)) %>% 
  select(
    variant          = HGVSp_VEP,
    occurrence_prob  = adj_occurrence_prob,
    clinvar          = clinvar_clnsig,
    score
  )

# ---- output ----
kable(summary_df_w2, digits = 9, caption = "candidate vs missing variants, grouped by score")
kable(variant_table,   digits = 9, caption = "weighted missing variants (detailed)")




# ---- detailed variant‑level table (candidate causal + missing only) ----
variant_table_all <- patient_data %>% 
  mutate(flag = if_else(known_var_observed == -9, "missing", "present")) %>% 
  filter(flag == "missing" | (flag == "present" & score > 3)) %>% 
  arrange(flag, desc(score), desc(adj_occurrence_prob)) %>% 
  select(
    flag,
    variant          = HGVSp_VEP,
    occurrence_prob  = adj_occurrence_prob,
    clinvar          = clinvar_clnsig,
    score
  )

# ---- output ----
kable(summary_df_w2, digits = 9,
      caption = "candidate vs missing variants, grouped by score")

kable(variant_table_all, digits = 9,
      caption = "candidate causal and missing variants (detailed)")






















# wilson method ----

# ---- 95 % credible interval for candidate dominance ----
set.seed(666)

# ---- refine priors with score ----
prior_weight <- function(s) ifelse(s > 0, s + 1, 1)

patient_data <- patient_data %>% 
  mutate(
    prior_w = prior_weight(score)
  )

# candidate prior
cand <- patient_data %>% filter(known_var_observed == 1)
N_c   <- cand$max_an
k_c   <- round(cand$adj_occurrence_prob * N_c)
alpha_c <- k_c + cand$prior_w
beta_c  <- N_c - k_c + 1

# missing priors (only score ≥ 0 variants carry risk)
miss <- patient_data %>% 
  filter(known_var_observed == -9, score >= 0)

alpha_m <- round(miss$adj_occurrence_prob * miss$max_an) + miss$prior_w
beta_m  <- miss$max_an - round(miss$adj_occurrence_prob * miss$max_an) + 1

# ---- posterior simulation ----
n_sim <- 10000
cand_sim <- rbeta(n_sim, alpha_c, beta_c)
miss_sim <- if (nrow(miss)) {
  sapply(seq_along(alpha_m), \(i) rbeta(n_sim, alpha_m[i], beta_m[i]))
} else {
  matrix(0, nrow = n_sim, ncol = 1)
}

miss_any_sim <- 1 - apply(1 - miss_sim, 1, prod)
rel_conf_sim <- cand_sim / (cand_sim + miss_any_sim)
ci_rel_conf  <- quantile(rel_conf_sim, c(0.025, 0.5, 0.975))

# ---- output credible interval ----
kable(
  tibble(
    quantile = c("2.5 %", "median", "97.5 %"),
    rel_conf = ci_rel_conf
  ),
  digits = 4,
  caption = "credible interval for candidate variant being top causal (95 %)"
)













# ---- identify present (observed) alternate variants (>0 score) ----
present_alt <- patient_data %>% 
  filter(known_var_observed %in% c(1, 2), HGVSp_VEP != "p.Ser237Ter", score > 0)

# ---- posterior parameters for all observed variants ----
obs_variants <- bind_rows(
  candidate = cand,
  present_alt
) %>% 
  mutate(var_id = row_number())

alpha_obs <- round(obs_variants$adj_occurrence_prob * obs_variants$max_an) + obs_variants$prior_w
beta_obs  <- obs_variants$max_an - round(obs_variants$adj_occurrence_prob * obs_variants$max_an) + 1

# ---- combine missing (score ≥0) ----
miss <- patient_data %>% 
  filter(known_var_observed == -9, score >= 0)

alpha_m <- round(miss$adj_occurrence_prob * miss$max_an) + miss$prior_w
beta_m  <- miss$max_an - round(miss$adj_occurrence_prob * miss$max_an) + 1

# ---- posterior simulation ----
set.seed(777)
n_sim <- 10000

obs_sim  <- sapply(seq_along(alpha_obs), \(i) rbeta(n_sim, alpha_obs[i],  beta_obs[i]))
miss_sim <- if (nrow(miss)) {
  sapply(seq_along(alpha_m),   \(i) rbeta(n_sim, alpha_m[i],    beta_m[i]))
} else {
  matrix(0, nrow = n_sim, ncol = 1)
}

miss_any_sim <- 1 - apply(1 - miss_sim, 1, prod)

# ---- probability candidate has highest posterior ----
cand_sim <- obs_sim[, 1]
alt_sim  <- if (ncol(obs_sim) > 1) obs_sim[, -1, drop = FALSE] else matrix(0, nrow = n_sim, ncol = 1)
max_other <- pmax(miss_any_sim, apply(alt_sim, 1, max))

prob_cand_top <- mean(cand_sim > max_other)
ci_cand_top   <- quantile(cand_sim / (cand_sim + max_other), c(0.025, 0.5, 0.975))

# ---- summary output ----
kable(
  tibble(
    metric    = c("probability candidate highest", "CI2.5", "CImedian", "CI97.5"),
    value     = c(prob_cand_top, ci_cand_top)
  ),
  digits = 4,
  caption = "confidence that candidate variant is top causal (observed + missing)"
)


# plot ----

# ---- credible interval data ----
ci_df <- tibble(
  variant = "candidate",
  lower   = ci_cand_top[1],
  median  = ci_cand_top[2],
  upper   = ci_cand_top[3]
)

gg_ci <- ggplot(ci_df, aes(x = variant, y = median, ymin = lower, ymax = upper)) +
  geom_pointrange(colour = "steelblue", linewidth = 0.8) +
  labs(y = "relative confidence", x = NULL) +
  coord_flip() +
  theme_bw()

# ---- posterior mean probability for each variant ----
posterior_mean <- function(a, b) a / (a + b)

var_probs <- bind_rows(
  cand  %>% mutate(prob = posterior_mean(alpha_c, beta_c)) %>% select(variant = HGVSp_VEP, prob),
  present_alt %>% mutate(prob = posterior_mean(alpha_obs[-1], beta_obs[-1])) %>% select(variant = HGVSp_VEP, prob)
)

if (nrow(miss)) {
  var_probs <- bind_rows(
    var_probs,
    miss %>% 
      mutate(variant = paste0("missing_", HGVSp_VEP),
             prob    = posterior_mean(alpha_m, beta_m)) %>% 
      select(variant, prob)
  )
}

gg_var <- ggplot(var_probs, aes(x = reorder(variant, prob), y = prob)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(x = NULL, y = "posterior mean probability") +
  theme_bw()

# ---- combined plot ----
(gg_ci | gg_var) + plot_annotation(tag_levels = "a")




















# ---- assemble combined dataset ----
ci_all <- bind_rows(
  # candidate CI
  tibble(
    variant = "candidate",
    lower   = ci_cand_top[1],
    median  = ci_cand_top[2],
    upper   = ci_cand_top[3],
    type    = "candidate"
  ),
  # other variants (point only)
  var_probs %>% 
    mutate(
      lower = prob,
      median = prob,
      upper = prob,
      type = "other"
    )
)

# order variants by median
ci_all <- ci_all %>% mutate(variant = reorder(variant, median))

# ---- single plot ----
ggplot(ci_all, aes(x = variant, y = median, ymin = lower, ymax = upper,
                   colour = type, shape = type)) +
  geom_pointrange(size = 1) +
  scale_colour_manual(values = c(candidate = "steelblue", other = "grey40")) +
  scale_shape_manual(values  = c(candidate = 17,          other = 16)) +
  coord_flip() +
  labs(x = NULL, y = "probability", colour = NULL, shape = NULL) +
  theme_bw()






















library(tidytext)   # for reorder_within / scale_y_reordered

# ---- candidate CI row ----
cand_label <- paste0("present_", cand$HGVSp_VEP)
ci_candidate <- tibble(
  variant = cand_label,
  lower   = ci_cand_top[1],
  median  = ci_cand_top[2],
  upper   = ci_cand_top[3],
  type    = "present"
)

# ---- other variants (point‑only rows) ----
other_variants <- var_probs %>% 
  filter(variant != cand_label) %>% 
  mutate(
    lower = prob,
    median = prob,
    upper = prob,
    type  = if_else(str_detect(variant, "^missing_"), "missing", "present")
  )

# ---- combined dataset ----
ci_all <- bind_rows(ci_candidate, other_variants) %>% 
  mutate(
    var_plot = reorder_within(variant, median, type)
  )

# ---- plot ----
ggplot(ci_all,
       aes(x = median, xmin = lower, xmax = upper,
           y = var_plot, colour = type, shape = type)) +
  geom_pointrange(size = 0.8) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(present = "steelblue", missing = "grey40")) +
  scale_shape_manual(values  = c(present = 17,         missing = 16)) +
  scale_y_reordered() +
  labs(x = "probability", y = NULL, colour = NULL, shape = NULL) +
  theme_bw()



















library(dplyr)
library(ggplot2)
library(tidytext)    # for reorder_within / scale_y_reordered
library(forcats)

# helper
posterior_mean <- function(a, b) a / (a + b)
prior_weight   <- function(s) ifelse(s > 0, s + 1, 1)

# ----- candidate row (CI) ----------------------------------------------------
cand_label <- paste0("present_", cand$HGVSp_VEP)

ci_candidate <- tibble(
  variant = cand_label,
  lower   = ci_cand_top[1],
  median  = ci_cand_top[2],
  upper   = ci_cand_top[3],
  type    = "present"
)

# ----- other present variants (point only) -----------------------------------
present_alt_rows <- present_alt %>% 
  mutate(
    variant = paste0("present_", HGVSp_VEP),
    prob    = posterior_mean(alpha_obs[-1], beta_obs[-1]),
    lower   = prob,
    median  = prob,
    upper   = prob,
    type    = "present"
  ) %>% 
  select(variant, lower, median, upper, type)

# ----- all missing variants (include benign) ---------------------------------
miss_all <- patient_data %>% 
  filter(known_var_observed == -9) %>% 
  mutate(
    variant = paste0("missing_", HGVSp_VEP),
    prior_w = prior_weight(score),
    alpha   = round(adj_occurrence_prob * max_an) + prior_w,
    beta    = max_an - round(adj_occurrence_prob * max_an) + 1,
    prob    = posterior_mean(alpha, beta),
    lower   = prob,
    median  = prob,
    upper   = prob,
    type    = "missing"
  ) %>% 
  select(variant, lower, median, upper, type)

# ----- combine & order -------------------------------------------------------
ci_all <- bind_rows(ci_candidate, present_alt_rows, miss_all) %>% 
  mutate(var_plot = reorder_within(variant, median, type))

# ----- plot ------------------------------------------------------------------
ggplot(ci_all,
       aes(x = median, xmin = lower, xmax = upper,
           y = var_plot, colour = type, shape = type)) +
  geom_pointrange(size = 0.8) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(present = "steelblue", missing = "grey40")) +
  scale_shape_manual(values  = c(present = 17,         missing = 16)) +
  scale_y_reordered() +
  labs(x = "probability", y = NULL, colour = NULL, shape = NULL) +
  theme_bw()





























# single causal variant ----
# ----- libraries -------------------------------------------------------------
library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(ggplot2); theme_set(theme_bw())
library(knitr)
library(tidytext)      # for reorder_within / scale_y_reordered
library(forcats)

# ----------------------------------------------------------------------------- 
#                               PRE‑PROCESSING
# ----------------------------------------------------------------------------- 
#  Assumes patient_data already created as in prior steps, with `score`,
#  `occurrence_prob`, `max_an`, and `known_var_observed` columns.

# pathogenicity‑weighted occurrence probability
patient_data <- patient_data %>% 
  mutate(
    pathogenic_weight   = scales::rescale(score, to = c(0, 1), from = c(-5, 5)),
    adj_occurrence_prob = occurrence_prob * pathogenic_weight
  )

# ----------------------------------------------------------------------------- 
#                        SUMMARY TABLES  (weighted probs)
# ----------------------------------------------------------------------------- 
# probability any missing variant present (weighted)
p_missing_any_w <- patient_data %>% 
  filter(known_var_observed == -9) %>% 
  pull(adj_occurrence_prob) %>% 
  { 1 - prod(1 - .) }

# candidate variant (observed pathogenic, score > 3)
cand <- patient_data %>% filter(known_var_observed == 1, score > 3)
p_candidate_w <- cand$adj_occurrence_prob

# per‑score summary of missing variants
missing_score_summary <- patient_data %>% 
  filter(known_var_observed == -9) %>% 
  group_by(score) %>% 
  summarise(
    n_variants = n(),
    prob_any   = 1 - prod(1 - adj_occurrence_prob),
    .groups    = "drop"
  ) %>% 
  arrange(desc(score))

total_prob <- p_candidate_w + p_missing_any_w

summary_df_w2 <- bind_rows(
  tibble(
    category    = "candidate variant (weighted)",
    probability = p_candidate_w,
    proportion  = p_candidate_w / total_prob
  ),
  missing_score_summary %>% 
    transmute(
      category    = paste0("missing variants score ", score, " (weighted)"),
      probability = prob_any,
      proportion  = prob_any / total_prob
    )
)

# detailed variant‑level table (candidate causal & missing)
variant_table_all <- patient_data %>% 
  mutate(flag = if_else(known_var_observed == -9, "missing", "present")) %>% 
  filter(flag == "missing" | (flag == "present" & score > 3)) %>% 
  arrange(flag, desc(score), desc(adj_occurrence_prob)) %>% 
  select(
    flag,
    variant          = HGVSp_VEP,
    occurrence_prob  = adj_occurrence_prob,
    clinvar          = clinvar_clnsig,
    score
  )

# display summary tables
kable(summary_df_w2, digits = 9,
      caption = "candidate vs missing variants, grouped by score")
kable(variant_table_all, digits = 9,
      caption = "candidate causal and missing variants (detailed)")

# ----------------------------------------------------------------------------- 
#         CREDIBLE INTERVAL & PROBABILITY CANDIDATE IS TOP CAUSAL
# ----------------------------------------------------------------------------- 
set.seed(666)

prior_weight <- function(s) ifelse(s > 0, s + 1, 1)
posterior_mean <- function(a, b) a / (a + b)

patient_data <- patient_data %>% mutate(prior_w = prior_weight(score))

# candidate posterior parameters
N_c   <- cand$max_an
k_c   <- round(cand$adj_occurrence_prob * N_c)
alpha_c <- k_c + cand$prior_w
beta_c  <- N_c - k_c + 1

# other present variants with score > 3 (exclude candidate)
present_alt <- patient_data %>% 
  filter(known_var_observed %in% c(1, 2),
         HGVSp_VEP != cand$HGVSp_VEP,
         score > 3)

# posterior parameters for all observed variants
obs_variants <- bind_rows(cand, present_alt)
alpha_obs <- round(obs_variants$adj_occurrence_prob * obs_variants$max_an) + 
  obs_variants$prior_w
beta_obs  <- obs_variants$max_an - 
  round(obs_variants$adj_occurrence_prob * obs_variants$max_an) + 1

# all missing variants (include benign to show but risk only if score ≥0)
miss_all <- patient_data %>% filter(known_var_observed == -9)
miss_risk <- miss_all %>% filter(score >= 0)
alpha_m <- round(miss_risk$adj_occurrence_prob * miss_risk$max_an) + 
  miss_risk$prior_w
beta_m  <- miss_risk$max_an - 
  round(miss_risk$adj_occurrence_prob * miss_risk$max_an) + 1

# posterior simulation
set.seed(777)
n_sim <- 10000
obs_sim  <- sapply(seq_along(alpha_obs), \(i) rbeta(n_sim, alpha_obs[i], beta_obs[i]))
miss_sim <- if (nrow(miss_risk)) {
  sapply(seq_along(alpha_m), \(i) rbeta(n_sim, alpha_m[i], beta_m[i]))
} else {
  matrix(0, nrow = n_sim, ncol = 1)
}

miss_any_sim <- 1 - apply(1 - miss_sim, 1, prod)
cand_sim     <- obs_sim[, 1]
alt_sim      <- if (ncol(obs_sim) > 1) obs_sim[, -1, drop = FALSE] else 
  matrix(0, nrow = n_sim, ncol = 1)
max_other    <- pmax(miss_any_sim, apply(alt_sim, 1, max))

prob_cand_top <- mean(cand_sim > max_other)
ci_cand_top   <- quantile(cand_sim / (cand_sim + max_other),
                          c(0.025, 0.5, 0.975))

kable(
  tibble(
    metric = c("probability candidate highest", "CI2.5", "CImedian", "CI97.5"),
    value  = c(prob_cand_top, ci_cand_top)
  ),
  digits = 4,
  caption = "confidence that candidate variant is top causal (observed + missing)"
)

# ----------------------------------------------------------------------------- 
#                           PLOT (CI + points)
# ----------------------------------------------------------------------------- 
# candidate CI row
cand_label <- paste0("present_", cand$HGVSp_VEP)
ci_candidate <- tibble(
  variant = cand_label,
  lower   = ci_cand_top[1],
  median  = ci_cand_top[2],
  upper   = ci_cand_top[3],
  type    = "present"
)

# other present variants (point only)
present_alt_rows <- present_alt %>% 
  mutate(
    variant = paste0("present_", HGVSp_VEP),
    prob    = posterior_mean(alpha_obs[-1], beta_obs[-1]),
    lower   = prob,
    median  = prob,
    upper   = prob,
    type    = "present"
  ) %>% 
  select(variant, lower, median, upper, type)

# missing variants (point only, include benign)
miss_rows <- miss_all %>% 
  mutate(
    variant = paste0("missing_", HGVSp_VEP),
    prob    = posterior_mean(
      round(adj_occurrence_prob * max_an) + prior_w,
      max_an - round(adj_occurrence_prob * max_an) + 1),
    lower   = prob,
    median  = prob,
    upper   = prob,
    type    = "missing"
  ) %>% 
  select(variant, lower, median, upper, type)

ci_all <- bind_rows(ci_candidate, present_alt_rows, miss_rows) %>% 
  mutate(var_plot = reorder_within(variant, median, type))

ggplot(ci_all,
       aes(x = median, xmin = lower, xmax = upper,
           y = var_plot, colour = type, shape = type)) +
  geom_pointrange(size = 0.8) +
  facet_grid(type ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(present = "steelblue", missing = "grey40")) +
  scale_shape_manual(values  = c(present = 17,         missing = 16)) +
  scale_y_reordered() +
  labs(x = "probability", y = NULL, colour = NULL, shape = NULL) +
  theme_bw()














# multiple causal -----
# ----------------------------------------------------------------------------- 
# COMPLETE SCRIPT — handles multiple present causal variants (score > 3)
# ----------------------------------------------------------------------------- 
library(dplyr)
library(stringr)
library(scales)
library(ggplot2); theme_set(theme_bw())
library(knitr)
library(tidytext)
library(forcats)

set.seed(666)

prior_weight   <- function(s) ifelse(s > 0, s + 1, 1)
posterior_mean <- function(a, b) a / (a + b)

# ----------------------------------------------------------------- patient_data
patient_data <- varRisEst_var_scored |>
  filter(genename == "NFKB1") |>
  mutate(
    known_var_observed = case_when(
      HGVSp_VEP == "p.Ser237Ter"   ~ 1L,  # pathogenic
      HGVSc_VEP == "c.159+1G>A"    ~ 1L,  # likely pathogenic
      HGVSp_VEP == "p.Val236Ile"   ~ -9L, # missing
      HGVSp_VEP == "p.Thr567Ile"   ~ -9L, # missing
      TRUE                         ~ 0L
    ),
    pathogenic_weight   = rescale(score, to = c(0, 1), from = c(-5, 5)),
    adj_occurrence_prob = occurrence_prob * pathogenic_weight,
    prior_w             = prior_weight(score)
  )


# ----------------------------------------------------------------- patient_data
patient_data <- varRisEst_var_scored |>
  filter(genename == "NFKB1") |>
  mutate(
    known_var_observed = case_when(
      HGVSp_VEP == "p.Ser237Ter"   ~ 1L,  # pathogenic
      HGVSc_VEP == "c.159+1G>A"    ~ -9L,  # likely pathogenic missing
      HGVSp_VEP == "p.Val236Ile"   ~ -9L, # benign missing
      HGVSp_VEP == "p.Thr567Ile"   ~ -9L, # benign missing
      TRUE                         ~ 0L
    ),
    pathogenic_weight   = rescale(score, to = c(0, 1), from = c(-5, 5)),
    adj_occurrence_prob = occurrence_prob * pathogenic_weight,
    prior_w             = prior_weight(score)
  )

# # --------------------------------------------------------------- risk_variants
# risk_variants <- patient_data %>% 
#   filter((known_var_observed == 1 & score > 3) | known_var_observed == -9) %>% 
#   mutate(
#     flag        = if_else(known_var_observed == -9, "missing", "present"),
#     variant_lab = paste0(flag, "_", coalesce(HGVSp_VEP, HGVSc_VEP)),
#     alpha       = round(adj_occurrence_prob * max_an) + prior_w,
#     beta        = max_an - round(adj_occurrence_prob * max_an) + 1
#   )

# --------------------------------------------------------------- risk_variants
risk_variants <- patient_data %>% 
  filter((known_var_observed == 1 & score > 3) | known_var_observed == -9) %>% 
  mutate(
    flag = if_else(known_var_observed == -9, "missing", "present"),
    variant_id = if_else(HGVSp_VEP == ".", HGVSc_VEP, HGVSp_VEP),
    variant_lab = paste0(flag, "_", variant_id),
    alpha = round(adj_occurrence_prob * max_an) + prior_w,
    beta  = max_an - round(adj_occurrence_prob * max_an) + 1
  )

# ---------------------------------------------------------------- probabilities
n_sim <- 10000
sim_mat <- sapply(seq_len(nrow(risk_variants)),
                  \(i) rbeta(n_sim, risk_variants$alpha[i], risk_variants$beta[i]))
colnames(sim_mat) <- risk_variants$variant_lab

top_hit        <- apply(sim_mat, 1, which.max)
prob_is_topvec <- prop.table(table(factor(top_hit, levels = seq_len(nrow(risk_variants)))))

risk_variants <- risk_variants %>% 
  mutate(prob_is_top = as.numeric(prob_is_topvec))

# credible intervals for each present variant’s posterior share
ci_present <- risk_variants %>% 
  filter(flag == "present") %>% 
  mutate(
    lower  = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) quantile(x, 0.025)),
    median = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) quantile(x, 0.500)),
    upper  = apply(sim_mat[, variant_lab, drop = FALSE] /
                     rowSums(sim_mat), 2, \(x) quantile(x, 0.975))
  ) %>% 
  select(variant_lab, lower, median, upper)

# attach CI / point values
plot_df <- risk_variants %>% 
  left_join(ci_present, by = c("variant_lab")) %>% 
  mutate(
    lower  = if_else(flag == "present", lower, prob_is_top),
    median = if_else(flag == "present", median, prob_is_top),
    upper  = if_else(flag == "present", upper, prob_is_top),
    var_plot = reorder_within(variant_lab, median, flag)
  )

# ------------------------------------------------------------- summary tables
summary_df <- risk_variants %>% 
  group_by(flag) %>% 
  summarise(probability = sum(adj_occurrence_prob),
            .groups = "drop") %>% 
  mutate(proportion = probability / sum(probability))

variant_table <- risk_variants %>% 
  arrange(flag, desc(score), desc(adj_occurrence_prob)) %>% 
  select(flag,
         variant     = variant_lab,
         occurrence_prob = adj_occurrence_prob,
         clinvar     = clinvar_clnsig,
         score)

kable(summary_df, digits = 9,
      caption = "candidate vs missing variant groups (weighted)")
kable(variant_table, digits = 9,
      caption = "candidate causal and missing variants (detailed)")

# --------------------------------------------------------------------------- plot
ggplot(plot_df,
       aes(x = median, xmin = lower, xmax = upper,
           y = var_plot, colour = flag, shape = flag)) +
  geom_pointrange(size = 0.8) +
  facet_grid(flag ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(present = "steelblue", missing = "grey40")) +
  scale_shape_manual(values  = c(present = 17,          missing = 16)) +
  scale_y_reordered() +
  labs(x = "probability", y = NULL, colour = NULL, shape = NULL) +
  theme_bw()

