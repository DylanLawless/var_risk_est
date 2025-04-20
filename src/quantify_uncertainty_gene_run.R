library(dplyr)
library(stringr)
library(purrr)
library(scales)
library(forcats)
library(tidytext)
library(ggplot2)
library(ggridges)
library(patchwork)
library(tidyr)
theme_set(theme_bw())

set.seed(666)

# — data loading & preprocessing ——————————————————————————————————————
PANEL          <- 398
path_data      <- "~/web/PanelAppRex/data"

varRisEst_var  <- readRDS(sprintf(
  "../output/VarRiskEst_PanelAppRex_ID_%d_gene_variants.Rds", PANEL
))

# ClinVar score map ---------------------------------------------------------
score_map <- c(
  Pathogenic = 5, `Likely pathogenic` = 4,
  `Pathogenic, low penetrance` = 3, `likely pathogenic, low penetrance` = 3,
  `Conflicting classifications of pathogenicity` = 2,
  `risk factor` = 1, association = 1, `likely risk allele` = 1,
  `drug response` = 0, `Uncertain significance` = 0,
  protective = -3, `Likely benign` = -4, Benign = -5
)
get_score <- function(cln) {
  terms <- str_split(str_replace_all(cln, "_", " "), "[/|]", simplify = TRUE) |> str_trim()
  mean(ifelse(terms %in% names(score_map), score_map[terms], 0))
}
prior_weight <- function(s) ifelse(s > 0, s + 1, 1)

# build patient data --------------------------------------------------------
build_patient_data <- function(include_missing_variant = FALSE) {
  varRisEst_var %>%
    filter(genename == genetic_defect) %>%
    mutate(
      score               = map_dbl(clinvar_clnsig, get_score),
      evidence_score      = score,
      known_var_observed  = case_when(
        HGVSp_VEP == "p.Ser237Ter" ~ 1L,
        HGVSp_VEP == "p.Arg231His" ~ 1L,
        HGVSp_VEP == "p.Gly650Arg" ~ 1L,
        HGVSp_VEP == "p.Val236Ile" ~ -9L,
        HGVSp_VEP == "p.Thr567Ile" ~ -9L,
        include_missing_variant & HGVSc_VEP == "c.159+1G>A" ~ -9L,
        TRUE ~ 0L
      ),
      pathogenic_weight   = rescale(score, to = c(0, 1), from = c(-5, 5)),
      adj_occurrence_prob = occurrence_prob * pathogenic_weight,
      prior_w             = prior_weight(score)
    )
}



# scenario 1
scenario     <- "1"
genetic_defect <- "NFKB1"
patient_data <- build_patient_data(include_missing_variant = FALSE)
source("quantify_uncertainty_gene_functions.R")

# scenario 2
scenario     <- "2"
genetic_defect <- "NFKB1"
patient_data <- build_patient_data(include_missing_variant = TRUE)
source("quantify_uncertainty_gene_functions.R")

# scenario 3
# set one variant from every classification as missing 
build_patient_data <- function(include_missing_variant = FALSE) {
  var <- varRisEst_var %>%
    filter(genename == genetic_defect) %>%
    mutate(score = map_dbl(clinvar_clnsig, get_score))
  
  # select one representative variant per classification to mark as -9
  mark_as_uncertain <- var %>%
    group_by(clinvar_clnsig) %>%
    slice(1) %>%
    ungroup() %>%
    pull(HGVSp_VEP)
  
  var %>%
    mutate(
      evidence_score      = score,
      known_var_observed  = case_when(
        HGVSp_VEP %in% c("p.Ser237Ter", "p.Arg231His", "p.Gly650Arg") ~ 1L,
        HGVSp_VEP %in% mark_as_uncertain ~ -9L,
        include_missing_variant & HGVSc_VEP == "c.159+1G>A" ~ -9L,
        TRUE ~ 0L
      ),
      pathogenic_weight   = rescale(score, to = c(0, 1), from = c(-5, 5)),
      adj_occurrence_prob = occurrence_prob * pathogenic_weight,
      prior_w             = prior_weight(score)
    )
}

scenario     <- "3"
genetic_defect <- "TNFAIP3"
patient_data <- build_patient_data(include_missing_variant = TRUE)
source("quantify_uncertainty_gene_functions.R")

# risk_variants  |> select(HGVSp_VEP, clinvar_clnsig, score) |> arrange(score)
