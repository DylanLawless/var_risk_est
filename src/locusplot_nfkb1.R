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

# Filter for nfkb1
nfkb1_df <- varRisEst_gene_scored %>%
  dplyr::filter(genename == "NFKB1") %>%
  mutate(
    chrom = "4",  # nfkb1 is on chr7
    pos = as.numeric(`pos(1-based)`),
    yvar = as.numeric(occurrence_prob)
  )

# vector of "chr7:position" strings
nfkb1_coords <- paste0("4[CHR] AND ", nfkb1_df$pos, "[CHRPOS] AND Homo sapiens[Organism]")

get_rs_from_pos <- function(query) {
  res <- entrez_search(db = "snp", term = query, retmax = 1)
  if (length(res$ids) > 0) {
    return(res$ids[[1]])
  } else {
    return(NA)
  }
}

# batch with progress
# rsids <- map_chr(nfkb1_coords, possibly(get_rs_from_pos, NA_character_))
# saveRDS(rsids, file="../output/VarRiskEst_PanelAppRex_ID_398_nfkb1_rsids.Rds")
rsids <- readRDS( file="../output/VarRiskEst_PanelAppRex_ID_398_nfkb1_rsids.Rds")

nfkb1_df$rsid <- rsids
nfkb1_df <- nfkb1_df %>% dplyr::filter(!is.na(rsid))  # keep only matched
nfkb1_df <- nfkb1_df %>% dplyr::filter(score >= 0)
nfkb1_df <- nfkb1_df |> dplyr::select(HGVSc_VEP, HGVSp_VEP, chrom, pos, rsid, occurrence_prob, score)

# assign fill colour based on raw probability
# loc_nfkb1$data$bg <- viridis::viridis(n = nrow(loc_nfkb1$data), option = "plasma")[rank(loc_nfkb1$data$occurrence_prob)]
nfkb1_df$occurrence_prob_log <- log10(nfkb1_df$occurrence_prob)

# linearly scale occurrence_prob_log to [1, N] for colour mapping
nfkb1_df <- nfkb1_df %>%
  mutate(
    occurrence_prob_log = log10(as.numeric(occurrence_prob)),
    bg = viridis::plasma(100)[
      rescale(occurrence_prob_log, to = c(1, 100)) |> round()
    ]
  )

# regenerate locus object
loc_nfkb1 <- locus(
  data = nfkb1_df,
  gene = "NFKB1",
  flank = 1e5,
  yvar = "occurrence_prob_log",
  ens_db = EnsDb.Hsapiens.v86
)
#
# add recombination rate
# loc_nfkb1 <- link_recomb(loc_nfkb1)
# saveRDS(loc_nfkb1, file="../output/VarRiskEst_PanelAppRex_ID_398_nfkb1_loc.Rds")
loc_nfkb1 <- readRDS( file="../output/VarRiskEst_PanelAppRex_ID_398_nfkb1_loc.Rds")

# add labels ----
# Get top 5 variants to label with HGVS
top_variants <- nfkb1_df %>%
  arrange(desc(occurrence_prob_log)) %>%
  slice_head(n = 5) %>%
  mutate(label = ifelse(HGVSp_VEP != ".", HGVSp_VEP, HGVSc_VEP))

# Compute offsets proportional to range (e.g. 10% of y range)
yrng <- range(nfkb1_df$occurrence_prob_log, na.rm = TRUE)
offset_unit <- diff(yrng) * 0.07  # adjustable spread factor

# Use staggered pattern of y-offsets
offsets <- rep(c(0, 1, -1, 2, -2), length.out = nrow(top_variants)) * offset_unit

library(scales)

nfkb1_df <- nfkb1_df %>%
  mutate(
    size = log10(score + 1) + 1,
    colour_val = score,
    bg = scales::col_numeric(
      palette = scales::gradient_n_pal(c("lightblue", "red"))(seq(0, 1, length.out = 100)),
      domain = c(0, 5)
    )(colour_val)
  )

yrng <- range(nfkb1_df$occurrence_prob_log, na.rm = TRUE)

top_variants <- nfkb1_df |>
  arrange(desc(occurrence_prob_log)) |>
  slice_head(n = 5) |>
  mutate(label = ifelse(HGVSp_VEP != ".", HGVSp_VEP, HGVSc_VEP))

n <- nrow(top_variants)

xrange <- range(top_variants$pos, na.rm = TRUE)
x_spread_factor <- 0.02
x_offsets <- seq(-1, 1, length.out = n) * diff(xrange) * x_spread_factor
top_variants$pos_offset <- top_variants$pos + x_offsets

y_spread_factor <- 0.4
y_offsets <- seq(-1, 1, length.out = n) * diff(yrng) * y_spread_factor
top_variants$y_offset <- top_variants$occurrence_prob_log + y_offsets

legend_cols <- scales::gradient_n_pal(c("lightblue", "red"))(seq(0, 1, length.out = 5))
legend_vals <- signif(seq(0, 5, length.out = 5), digits = 1)

pl <- quote({
  with(top_variants, {
    text(
      pos_offset,
      y_offset,
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

pdf("../images/locusplot_nfkb1.pdf", width = 5, height = 4)
locus_plot(
  loc_nfkb1,
  bg = nfkb1_df$bg,
  ylab = "Prob (-log10)",
  pch = 21,
  col = "black",
  cex = nfkb1_df$size,
  pcutoff = NULL,
  ylim = yrng + c(-0.5, 1),
  panel.last = pl,
  cex.axis = .7,
  cex.lab = 1.1,
  cex.main = 1,
  cex.text = 1.1
)
dev.off()
 
