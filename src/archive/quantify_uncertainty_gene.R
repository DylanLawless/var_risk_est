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
