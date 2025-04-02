# Read input data
df <- read_tsv("../output/VarRiskEst_PanelAppRex_ID_398_gene_variants.tsv", show_col_types = FALSE)

# map scores ----
library(scales)   # for rescale()
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2);theme_set(theme_bw())
library(biomaRt)
library(rentrez)
library(tibble)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(viridis)

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
varRisEst_gene_scored <- df %>%
  mutate(score = map_dbl(clinvar_clnsig, get_score))

# Filter for CFTR
cftr_df <- varRisEst_gene_scored %>%
  dplyr::filter(genename == "CFTR") %>%
  mutate(
    chrom = "7",  # CFTR is on chr7
    pos = as.numeric(`pos(1-based)`),
    yvar = as.numeric(occurrence_prob)
  )

# Special condition for displaying CFTR variants as common practice -----
# the transcript used on CFTR2 database and typically reported is transcript 2 in the raw data. We can source these labels from the version prepared in the standalone CFTR tests. By default the full HPC run only keeps transcript 1 for simplicity.
header_line <- readLines("../data/cftr_head", n = 1)
header_line <- sub("^#", "", header_line)
header_fields <- strsplit(header_line, "\t")[[1]]
rm(header_line)

df_hgvs2 <- read.table("../data/cftr", 
                 sep = "\t",
                 header = FALSE, 
                 stringsAsFactors = FALSE, 
                 fill = TRUE)
colnames(df_hgvs2) <- header_fields
df_hgvs2 <- df_hgvs2 %>% dplyr::select(genename, `pos(1-based)`, HGVSc_VEP, HGVSp_VEP)
# keep just one transcript allele for simplicity
df_hgvs2$HGVSc_VEP_1 <- sapply(strsplit(df_hgvs2$HGVSc_VEP, ";"), `[`, 1)
df_hgvs2$HGVSc_VEP_2 <- sapply(strsplit(df_hgvs2$HGVSc_VEP, ";"), `[`, 2)
df_hgvs2$HGVSp_VEP_1 <- sapply(strsplit(df_hgvs2$HGVSp_VEP, ";"), `[`, 1)
df_hgvs2$HGVSp_VEP_2 <- sapply(strsplit(df_hgvs2$HGVSp_VEP, ";"), `[`, 2)
df_hgvs2$genename <- sapply(strsplit(df_hgvs2$genename, ";"), `[`, 1)
# to replace cftr_df HGVSc_VEP HGVSp_VEP which match to HGVSc_VEP_1 HGVSc_VEP_1
df_hgvs2_dedup <- df_hgvs2 %>%
  distinct(HGVSc_VEP_1, HGVSp_VEP_1, .keep_all = TRUE)

cftr_df2 <- cftr_df %>%
  left_join(
    df_hgvs2_dedup %>%
      dplyr::select(HGVSc_VEP_1, HGVSp_VEP_1, HGVSc_VEP_2, HGVSp_VEP_2),
    by = c("HGVSc_VEP" = "HGVSc_VEP_1", "HGVSp_VEP" = "HGVSp_VEP_1")
  ) %>%
  mutate(
    HGVSc_VEP = ifelse(!is.na(HGVSc_VEP_2), HGVSc_VEP_2, HGVSc_VEP),
    HGVSp_VEP = ifelse(!is.na(HGVSp_VEP_2), HGVSp_VEP_2, HGVSp_VEP)
  ) %>%
  dplyr::select(-HGVSc_VEP_2, -HGVSp_VEP_2)

cftr_df <- cftr_df2
rm(cftr_df2)
# End CFTR special transcript handling -----

# vector of "chr7:position" strings
cftr_coords <- paste0("7[CHR] AND ", cftr_df$pos, "[CHRPOS] AND Homo sapiens[Organism]")

get_rs_from_pos <- function(query) {
  res <- entrez_search(db = "snp", term = query, retmax = 1)
  if (length(res$ids) > 0) {
    return(res$ids[[1]])
  } else {
    return(NA)
  }
}

# batch with progress
# rsids <- map_chr(cftr_coords, possibly(get_rs_from_pos, NA_character_))
# saveRDS(rsids, file="../output/VarRiskEst_PanelAppRex_ID_398_cftr_rsids.Rds")
rsids <- readRDS( file="../output/VarRiskEst_PanelAppRex_ID_398_cftr_rsids.Rds")

cftr_df$rsid <- rsids
cftr_df <- cftr_df %>% dplyr::filter(!is.na(rsid))  # keep only matched
cftr_df <- cftr_df %>% dplyr::filter(score > 0)
cftr_df <- cftr_df |> dplyr::select(HGVSc_VEP, HGVSp_VEP, chrom, pos, rsid, occurrence_prob, score)

# assign fill colour based on raw probability
# loc_cftr$data$bg <- viridis::viridis(n = nrow(loc_cftr$data), option = "inferno")[rank(loc_cftr$data$occurrence_prob)]
cftr_df$occurrence_prob_log <- log10(cftr_df$occurrence_prob)

# linearly scale occurrence_prob_log to [1, N] for colour mapping
cftr_df <- cftr_df %>%
  mutate(
    occurrence_prob_log = log10(as.numeric(occurrence_prob)),
    bg = viridis::inferno(100)[
      rescale(occurrence_prob_log, to = c(1, 100)) |> round()
    ]
  )

# regenerate locus object
loc_cftr <- locus(
  data = cftr_df,
  gene = "CFTR",
  flank = 1e5,
  yvar = "occurrence_prob_log",
  ens_db = EnsDb.Hsapiens.v86
)
#
# add recombination rate
# loc_cftr <- link_recomb(loc_cftr)
# saveRDS(loc_cftr, file="../output/VarRiskEst_PanelAppRex_ID_398_cftr_loc.Rds")
loc_cftr <- readRDS( file="../output/VarRiskEst_PanelAppRex_ID_398_cftr_loc.Rds")



cftr_df <- cftr_df %>%
  mutate(
    size = log10(score + 1) * 2,
    colour_val = score,  # for two-sided gradient
    bg = scales::col_numeric(
      palette = scales::gradient_n_pal(c("navy", "lightblue", "red"))(seq(0, 1, length.out = 100)),
      domain = c(-5, 5)
    )(colour_val)
  )

yrng <- range(cftr_df$occurrence_prob_log, na.rm = TRUE)

top_variants <- cftr_df |>
  arrange(desc(occurrence_prob_log)) |>
  slice_head(n = 5) |>
  mutate(label = ifelse(HGVSp_VEP != ".", HGVSp_VEP, HGVSc_VEP))

offsets <- rep(c(0, 1, -1, 2, -2), length.out = nrow(top_variants)) * diff(yrng) * 0.07

legend_cols <- scales::gradient_n_pal(c( "lightblue", "red"))(seq(0, 1, length.out = 5))
legend_vals <- signif(seq(0, 5, length.out = 5), digits = 1)

pl <- quote({
  with(top_variants, {
    text(
      pos,
      occurrence_prob_log + offsets,
      labels = label,
      pos = 3, cex = 1, 
      col = "black", offset = 0.3
    )
  })
  legend("topleft",
         legend = legend_vals,
         pt.bg = legend_cols,
         pch = 21,
         title = "Variant score",
         bty = "n",
         pt.cex = 1.5,
         y.intersp = 1.1,
         cex = 0.8
  )
})

pdf("../images/locusplot_cftr.pdf", width = 5, height = 4) 
locus_plot(
  loc_cftr,
  bg = cftr_df$bg,
  ylab = "Prob (-log10)",
  pch = 21,
  col = "black",
  cex = cftr_df$size,
  pcutoff = NULL,
  ylim = yrng + c(-0.5, 1),
  panel.last = pl,
  cex.axis = .7,
  cex.lab = 1.1,
  cex.main = 1,
  cex.text = 1.1
)
dev.off()

